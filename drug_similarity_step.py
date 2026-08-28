from __future__ import annotations

import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import BulkTanimotoSimilarity

from resource_config import DRUGBANK_SDF, DRUGBANK_XML, DRUGCENTRAL_TSV


NULL_STRINGS = {"", "nan", "na", "n/a", "null", "<na>", "-", "none"}

DRUGBANK_SHEET_COLUMNS = [
    "HuMANet_ID",
    "Query_Name",
    "Matched_Name",
    "Study_Folder",
    "Rank",
    "Drug",
    "DrugBank_ID",
    "Drug_Groups",
    "Similarity",
    "ATC_Codes",
    "ATC_Classes",
    "Indication",
    "Mechanism",
    "Therapeutic_Area",
    "Drug_InChIKey",
    "Drug_SMILES",
]

DRUGCENTRAL_SHEET_COLUMNS = [
    "HuMANet_ID",
    "Query_Name",
    "Matched_Name",
    "Study_Folder",
    "Rank",
    "Drug",
    "DrugCentral_ID",
    "CAS_RN",
    "Similarity",
    "Drug_InChIKey",
    "Drug_SMILES",
]

DRUGBANK_XML_CACHE = Path(DRUGBANK_XML).with_suffix(".metadata.pkl")


def _clean_scalar(value: object) -> Optional[str]:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text.lower() in NULL_STRINGS:
        return None
    return text


def _therapeutic_area_from_atc_classes(atc_classes: object) -> str:
    text = _clean_scalar(atc_classes)
    if not text:
        return "No ATC"
    text = text.upper()
    if "ALIMENTARY TRACT AND METABOLISM" in text:
        return "Metabolism"
    if "CARDIOVASCULAR SYSTEM" in text:
        return "Cardiovascular"
    if "NERVOUS SYSTEM" in text:
        return "Neurological"
    if "ANTINEOPLASTIC" in text:
        return "Cancer"
    if "ANTIINFECTIVES" in text or "ANTI-INFECTIVE" in text:
        return "Infectious disease"
    if "RESPIRATORY" in text:
        return "Respiratory"
    if "MUSCULOSKELETAL" in text:
        return "Musculoskeletal"
    if "BLOOD" in text:
        return "Blood"
    if "DERMATOLOGICAL" in text:
        return "Dermatology"
    if "GENITOURINARY" in text:
        return "Genitourinary"
    if "SENSORY ORGANS" in text:
        return "Sensory organs"
    return "No ATC"


def _first_prop(mol: Chem.Mol, names: Sequence[str]) -> Optional[str]:
    for name in names:
        if mol.HasProp(name):
            value = _clean_scalar(mol.GetProp(name))
            if value:
                return value
    return None


def _join_unique(values: Iterable[object], sep: str = "; ") -> Optional[str]:
    seen = set()
    out: List[str] = []
    for value in values:
        text = _clean_scalar(value)
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return sep.join(out) if out else None


def _mol_to_smiles(mol: Chem.Mol) -> Optional[str]:
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def _smiles_to_mol(value: object) -> Optional[Chem.Mol]:
    text = _clean_scalar(value)
    if not text:
        return None
    candidates = [part.strip() for part in text.split("|") if part.strip()]
    if not candidates:
        candidates = [text]
    for candidate in candidates:
        try:
            mol = Chem.MolFromSmiles(candidate)
        except Exception:
            mol = None
        if mol is not None:
            return mol
    return None


@lru_cache(maxsize=1)
def _fingerprint_generator():
    return rdFingerprintGenerator.GetMorganGenerator(radius=2, includeChirality=True)


def _count_fingerprint(mol: Chem.Mol):
    return _fingerprint_generator().GetCountFingerprint(mol)


@lru_cache(maxsize=1)
def _load_drugbank_xml_metadata() -> pd.DataFrame:
    xml_path = Path(DRUGBANK_XML)
    if not xml_path.exists():
        return pd.DataFrame(
            columns=[
                "DrugBank_ID",
                "ATC_Codes",
                "ATC_Classes",
                "Indication",
                "Mechanism",
            ]
        )

    if DRUGBANK_XML_CACHE.exists():
        return pd.read_pickle(DRUGBANK_XML_CACHE)

    ns = {"db": "http://www.drugbank.ca"}
    rows = []

    try:
        for _, elem in ET.iterparse(str(xml_path), events=("end",)):
            if elem.tag != f"{{{ns['db']}}}drug":
                continue

            drugbank_id = None
            for node in elem.findall("db:drugbank-id", ns):
                if node.attrib.get("primary") == "true":
                    drugbank_id = _clean_scalar(node.text)
                    break
            if not drugbank_id:
                elem.clear()
                continue

            atc_codes = []
            atc_classes = []
            for atc_node in elem.findall("db:atc-codes/db:atc-code", ns):
                code = _clean_scalar(atc_node.attrib.get("code"))
                if code:
                    atc_codes.append(code)
                for level_node in atc_node.findall("db:level", ns):
                    level_text = _clean_scalar(level_node.text)
                    if level_text:
                        atc_classes.append(level_text)

            indication = _clean_scalar(elem.findtext("db:indication", default=None, namespaces=ns))
            mechanism = _clean_scalar(elem.findtext("db:mechanism-of-action", default=None, namespaces=ns))

            rows.append(
                {
                    "DrugBank_ID": drugbank_id,
                    "ATC_Codes": _join_unique(atc_codes),
                    "ATC_Classes": _join_unique(atc_classes),
                    "Indication": indication,
                    "Mechanism": mechanism,
                }
            )
            elem.clear()
    except ET.ParseError as exc:
        print(f"[WARN] DrugBank XML could not be fully parsed ({xml_path}): {exc}. Falling back to SDF-only DrugBank metadata.")
        return pd.DataFrame(
            columns=[
                "DrugBank_ID",
                "ATC_Codes",
                "ATC_Classes",
                "Indication",
                "Mechanism",
            ]
        )

    out = pd.DataFrame(rows).drop_duplicates(subset=["DrugBank_ID"], keep="first")
    try:
        DRUGBANK_XML_CACHE.parent.mkdir(parents=True, exist_ok=True)
        out.to_pickle(DRUGBANK_XML_CACHE)
    except Exception:
        pass
    return out


@lru_cache(maxsize=1)
def _load_drugbank_hits_source() -> pd.DataFrame:
    sdf_paths = [Path(DRUGBANK_SDF)]
    rows = []
    for sdf_path in sdf_paths:
        if not sdf_path.exists():
            continue
        supplier = Chem.SDMolSupplier(
            str(sdf_path),
            removeHs=False,
            sanitize=True,
            strictParsing=False,
        )
        for mol in supplier:
            if mol is None:
                continue
            smiles = _mol_to_smiles(mol)
            if not smiles:
                continue
            rows.append(
                {
                    "DrugBank_ID": _first_prop(mol, ["DRUGBANK_ID", "DATABASE_ID"]),
                    "Drug": _first_prop(mol, ["GENERIC_NAME", "NAME"]),
                    "Drug_Groups": _first_prop(mol, ["DRUG_GROUPS"]),
                    "ATC_Codes": _first_prop(mol, ["ATC_CODES", "ATC_CODE", "ATC"]),
                    "ATC_Classes": _first_prop(mol, ["ATC_CLASSES", "ATC_CLASS", "ATC_DESCRIPTION"]),
                    "Indication": _first_prop(mol, ["INDICATION"]),
                    "Mechanism": _first_prop(mol, ["MECHANISM-OF-ACTION", "MECHANISM_OF_ACTION"]),
                    "Drug_InChIKey": _first_prop(mol, ["INCHI_KEY"]),
                    "Drug_SMILES": smiles,
                    "Mol": mol,
                }
            )

    if not rows:
        return pd.DataFrame(columns=DRUGBANK_SHEET_COLUMNS + ["Mol", "FP"])

    df = pd.DataFrame(rows)
    xml_meta = _load_drugbank_xml_metadata()
    if not xml_meta.empty:
        df = df.merge(xml_meta, on="DrugBank_ID", how="left", suffixes=("", "_xml"))
        for col in ["ATC_Codes", "ATC_Classes", "Indication", "Mechanism"]:
            xml_col = f"{col}_xml"
            if xml_col in df.columns:
                if col in df.columns:
                    df[col] = df[col].astype("string").fillna(df[xml_col].astype("string"))
                    df = df.drop(columns=[xml_col])
                else:
                    df = df.rename(columns={xml_col: col})
    df["Drug"] = df["Drug"].astype("string").fillna(df["DrugBank_ID"].astype("string"))
    df["Therapeutic_Area"] = df["ATC_Classes"].map(_therapeutic_area_from_atc_classes).astype("string")
    df["FP"] = df["Mol"].apply(_count_fingerprint)
    return df


@lru_cache(maxsize=1)
def _load_drugcentral_hits_source() -> pd.DataFrame:
    path = Path(DRUGCENTRAL_TSV)
    if not path.exists():
        return pd.DataFrame(columns=DRUGCENTRAL_SHEET_COLUMNS + ["Mol", "FP"])

    df = pd.read_csv(path, sep="\t", dtype="string", low_memory=False)
    required = {"SMILES", "InChIKey", "ID", "INN"}
    if not required.issubset(df.columns):
        return pd.DataFrame(columns=DRUGCENTRAL_SHEET_COLUMNS + ["Mol", "FP"])

    out = df.rename(
        columns={
            "SMILES": "Drug_SMILES",
            "InChIKey": "Drug_InChIKey",
            "ID": "DrugCentral_ID",
            "INN": "Drug",
        }
    ).copy()
    if "CAS_RN" not in out.columns:
        out["CAS_RN"] = pd.NA
    out["Mol"] = out["Drug_SMILES"].map(_smiles_to_mol)
    out = out[out["Mol"].notna()].copy()
    out["FP"] = out["Mol"].apply(_count_fingerprint)
    return out


def _prepare_input_metabolites(mapped_df: pd.DataFrame) -> pd.DataFrame:
    if mapped_df.empty:
        return pd.DataFrame(columns=["HuMANet_ID", "Query_Name", "Matched_Name", "Study_Folder", "SMILES", "Mol", "FP"])

    out = mapped_df.copy()
    if "HuMANet_ID" not in out.columns:
        if "Unnamed: 0" in out.columns:
            out["HuMANet_ID"] = out["Unnamed: 0"].astype("string")
        else:
            out["HuMANet_ID"] = pd.Series(
                [f"LIBRARY_{i:07d}" for i in range(len(out))],
                index=out.index,
                dtype="string",
            )
    if "Query_Name" not in out.columns:
        if "Standardized_Name" in out.columns:
            out["Query_Name"] = out["Standardized_Name"].astype("string")
        else:
            out["Query_Name"] = out["HuMANet_ID"].astype("string")
    if "matched_name" in out.columns:
        out["Matched_Name"] = out["matched_name"].astype("string")
    else:
        out["Matched_Name"] = pd.NA
    out["Matched_Name"] = out["Matched_Name"].fillna(out["Query_Name"].astype("string"))
    if "Study_Folder" not in out.columns:
        out["Study_Folder"] = pd.NA
    if "SMILES" not in out.columns:
        out["SMILES"] = pd.NA
    out["Mol"] = out["SMILES"].map(_smiles_to_mol)
    out = out[out["Mol"].notna()].copy()
    out["FP"] = out["Mol"].apply(_count_fingerprint)
    return out[["HuMANet_ID", "Query_Name", "Matched_Name", "Study_Folder", "SMILES", "Mol", "FP"]].copy()


def _top_hits(input_df: pd.DataFrame, ref_df: pd.DataFrame, top_n: int = 5, min_similarity: float = 0.0) -> list[tuple[int, int, float]]:
    if input_df.empty or ref_df.empty:
        return []
    ref_fps = ref_df["FP"].tolist()
    hits: list[tuple[int, int, float]] = []
    for left_idx, row in input_df.iterrows():
        sims = BulkTanimotoSimilarity(row["FP"], ref_fps)
        order = sorted(range(len(sims)), key=lambda i: sims[i], reverse=True)
        for rank, ref_pos in enumerate(order[:top_n], start=1):
            sim = float(sims[ref_pos])
            if sim < min_similarity:
                continue
            hits.append((left_idx, ref_pos, sim, rank))
    return hits


def build_drugbank_hits_sheet(mapped_df: pd.DataFrame, top_n: int = 5, min_similarity: float = 0.0) -> pd.DataFrame:
    inputs = _prepare_input_metabolites(mapped_df)
    drugbank = _load_drugbank_hits_source()
    if inputs.empty or drugbank.empty:
        return pd.DataFrame(columns=DRUGBANK_SHEET_COLUMNS)

    rows = []
    hits = _top_hits(inputs, drugbank, top_n=top_n, min_similarity=min_similarity)
    for left_idx, ref_pos, sim, rank in hits:
        left = inputs.loc[left_idx]
        right = drugbank.iloc[ref_pos]
        rows.append(
            {
                "HuMANet_ID": left["HuMANet_ID"],
                "Query_Name": left["Query_Name"],
                "Matched_Name": left["Matched_Name"],
                "Study_Folder": left["Study_Folder"],
                "Rank": rank,
                "Drug": right.get("Drug", pd.NA),
                "DrugBank_ID": right.get("DrugBank_ID", pd.NA),
                "Drug_Groups": right.get("Drug_Groups", pd.NA),
                "Similarity": round(sim, 4),
                "ATC_Codes": right.get("ATC_Codes", pd.NA),
                "ATC_Classes": right.get("ATC_Classes", pd.NA),
                "Indication": right.get("Indication", pd.NA),
                "Mechanism": right.get("Mechanism", pd.NA),
                "Therapeutic_Area": right.get("Therapeutic_Area", pd.NA),
                "Drug_InChIKey": right.get("Drug_InChIKey", pd.NA),
                "Drug_SMILES": right.get("Drug_SMILES", pd.NA),
            }
        )
    return pd.DataFrame(rows, columns=DRUGBANK_SHEET_COLUMNS).drop_duplicates().reset_index(drop=True)


def build_drugcentral_hits_sheet(mapped_df: pd.DataFrame, top_n: int = 5, min_similarity: float = 0.0) -> pd.DataFrame:
    inputs = _prepare_input_metabolites(mapped_df)
    drugcentral = _load_drugcentral_hits_source()
    if inputs.empty or drugcentral.empty:
        return pd.DataFrame(columns=DRUGCENTRAL_SHEET_COLUMNS)

    rows = []
    hits = _top_hits(inputs, drugcentral, top_n=top_n, min_similarity=min_similarity)
    for left_idx, ref_pos, sim, rank in hits:
        left = inputs.loc[left_idx]
        right = drugcentral.iloc[ref_pos]
        rows.append(
            {
                "HuMANet_ID": left["HuMANet_ID"],
                "Query_Name": left["Query_Name"],
                "Matched_Name": left["Matched_Name"],
                "Study_Folder": left["Study_Folder"],
                "Rank": rank,
                "Drug": right.get("Drug", pd.NA),
                "DrugCentral_ID": right.get("DrugCentral_ID", pd.NA),
                "CAS_RN": right.get("CAS_RN", pd.NA),
                "Similarity": round(sim, 4),
                "Drug_InChIKey": right.get("Drug_InChIKey", pd.NA),
                "Drug_SMILES": right.get("Drug_SMILES", pd.NA),
            }
        )
    return pd.DataFrame(rows, columns=DRUGCENTRAL_SHEET_COLUMNS).drop_duplicates().reset_index(drop=True)


def build_drug_annotation_sheets(mapped_df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    return {
        "drugbank_hits": build_drugbank_hits_sheet(mapped_df),
        "drugcentral_hits": build_drugcentral_hits_sheet(mapped_df),
    }
