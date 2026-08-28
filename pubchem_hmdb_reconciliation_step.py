import argparse
import re
import sqlite3
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
import requests

from knowledge_annotation_step import build_stage5_annotation_sheets
from pipeline_utils import (
    HAS_RPY2,
    _clean_cid,
    _clean_hmdb_id,
    _clean_inchikey,
    _clean_scalar,
    _first_existing_column,
    _first_hmdb_id,
    _normalize_hmdb_field,
    _split_ids,
    _split_synonyms,
    _sqlite_table_exists,
    annotate_classyfire,
    extract_ids_from_pugview,
    fetch_cids_by_name_sqlite,
    fetch_pugview_ids_for_cids,
    load_classyfire_bridge,
    timed_step,
    write_excel_workbook,
)
from resource_config import CLASSYFIRE_BRIDGE_R, HMDB_LITE_CSV, PUBCHEM_OFFLINE_SQLITE

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except Exception:
    HAS_RDKIT = False

def _normalize_chebi_field(val: object) -> Optional[str]:
    ids: List[str] = []
    for part in _split_ids(val):
        m = re.search(r"(?:CHEBI:)?(\d+)", str(part), flags=re.IGNORECASE)
        if m:
            ids.append(m.group(1))
    if not ids:
        return None
    seen = set()
    out = []
    for x in ids:
        if x in seen:
            continue
        seen.add(x)
        out.append(x)
    return ";".join(out)


def _naify_string_series(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .replace({"": pd.NA, "NA": pd.NA, "Na": pd.NA, "null": pd.NA, "NA_character_": pd.NA, "<NA>": pd.NA})
    )


def fetch_smiles_sqlite(db_path: str, cids: Set[str]) -> pd.DataFrame:
    if not cids:
        return pd.DataFrame(columns=["PubChem_CID", "SMILES"], dtype="string")

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows: List[Tuple[str, str]] = []
    cid_list = list(cids)
    for i in range(0, len(cid_list), 900):
        chunk = cid_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = f"SELECT cid, smiles FROM cid_smiles WHERE cid IN ({placeholders})"
        rows.extend(cur.execute(q, chunk).fetchall())
    con.close()

    return pd.DataFrame(rows, columns=["PubChem_CID", "SMILES"]).astype("string").drop_duplicates("PubChem_CID")


def fetch_inchikey_sqlite(db_path: str, cids: Set[str]) -> pd.DataFrame:
    if not cids:
        return pd.DataFrame(columns=["PubChem_CID", "InChIKey"], dtype="string")

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows: List[Tuple[str, str]] = []
    cid_list = list(cids)
    for i in range(0, len(cid_list), 900):
        chunk = cid_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = f"SELECT cid, inchikey FROM cid_inchi WHERE cid IN ({placeholders})"
        rows.extend(cur.execute(q, chunk).fetchall())
    con.close()

    return pd.DataFrame(rows, columns=["PubChem_CID", "InChIKey"]).astype("string").drop_duplicates("PubChem_CID")


def fetch_synonyms_sqlite(db_path: str, cids: Set[str], max_synonyms_per_cid: int = 200) -> pd.DataFrame:
    if not cids:
        return pd.DataFrame(columns=["PubChem_CID", "Synonyms"], dtype="string")

    con = sqlite3.connect(db_path)
    cur = con.cursor()

    syn_map: Dict[str, List[str]] = {cid: [] for cid in cids}
    seen: Dict[str, Set[str]] = {cid: set() for cid in cids}

    cid_list = list(cids)
    for i in range(0, len(cid_list), 900):
        chunk = cid_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = f"SELECT cid, syn FROM cid_syn WHERE cid IN ({placeholders})"
        for cid, syn in cur.execute(q, chunk):
            if len(syn_map[cid]) >= max_synonyms_per_cid:
                continue
            if syn in seen[cid]:
                continue
            seen[cid].add(syn)
            syn_map[cid].append(syn)

    con.close()

    rows = [(cid, "; ".join(vals) if vals else pd.NA) for cid, vals in syn_map.items()]
    return pd.DataFrame(rows, columns=["PubChem_CID", "Synonyms"]).astype("string")


def _sqlite_table_columns(db_path: str, table_name: str) -> List[str]:
    if not _sqlite_table_exists(db_path, table_name):
        return []
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows = cur.execute(f"PRAGMA table_info({table_name})").fetchall()
    con.close()
    return [r[1] for r in rows]


def fetch_pubchem_title_sqlite(db_path: str, cids: Set[str]) -> pd.DataFrame:
    if not cids or not _sqlite_table_exists(db_path, "cid_title"):
        return pd.DataFrame(columns=["PubChem_CID", "Standardized_Name"], dtype="string")

    cols = _sqlite_table_columns(db_path, "cid_title")
    cid_col = _first_existing_column(cols, ["cid"])
    title_col = _first_existing_column(cols, ["title", "standardized_name", "name"])
    if not cid_col or not title_col:
        return pd.DataFrame(columns=["PubChem_CID", "Standardized_Name"], dtype="string")

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows: List[Tuple[str, str]] = []
    cid_list = list(cids)
    for i in range(0, len(cid_list), 900):
        chunk = cid_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = f"SELECT {cid_col}, {title_col} FROM cid_title WHERE {cid_col} IN ({placeholders})"
        rows.extend(cur.execute(q, chunk).fetchall())
    con.close()

    return pd.DataFrame(rows, columns=["PubChem_CID", "Standardized_Name"]).astype("string").drop_duplicates("PubChem_CID")

def fetch_pubchem_iupac_sqlite(db_path: str, cids: Set[str]) -> pd.DataFrame:
    if not cids or not _sqlite_table_exists(db_path, "cid_iupac"):
        return pd.DataFrame(columns=["PubChem_CID", "IUPAC_Name"], dtype="string")

    cols = _sqlite_table_columns(db_path, "cid_iupac")
    cid_col = _first_existing_column(cols, ["cid"])
    iupac_col = _first_existing_column(cols, ["iupac_name", "iupac", "name"])
    if not cid_col or not iupac_col:
        return pd.DataFrame(columns=["PubChem_CID", "IUPAC_Name"], dtype="string")

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows: List[Tuple[str, str]] = []
    cid_list = list(cids)
    for i in range(0, len(cid_list), 900):
        chunk = cid_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = f"SELECT {cid_col}, {iupac_col} FROM cid_iupac WHERE {cid_col} IN ({placeholders})"
        rows.extend(cur.execute(q, chunk).fetchall())
    con.close()

    return pd.DataFrame(rows, columns=["PubChem_CID", "IUPAC_Name"]).astype("string").drop_duplicates("PubChem_CID")


def fetch_pubchem_mass_sqlite(db_path: str, cids: Set[str]) -> pd.DataFrame:
    if not cids or not _sqlite_table_exists(db_path, "cid_mass"):
        return pd.DataFrame(columns=["PubChem_CID", "Molecular_Formula", "Exact_Mass"], dtype="string")

    cols = _sqlite_table_columns(db_path, "cid_mass")
    cid_col = _first_existing_column(cols, ["cid"])
    formula_col = _first_existing_column(cols, ["molecular_formula", "formula"])
    mass_col = _first_existing_column(cols, ["exact_mass", "monoisotopic_mass", "mass"])
    if not cid_col:
        return pd.DataFrame(columns=["PubChem_CID", "Molecular_Formula", "Exact_Mass"], dtype="string")

    select_cols = [cid_col]
    out_cols = ["PubChem_CID"]
    if formula_col:
        select_cols.append(formula_col)
        out_cols.append("Molecular_Formula")
    if mass_col:
        select_cols.append(mass_col)
        out_cols.append("Exact_Mass")

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows = []
    cid_list = list(cids)
    for i in range(0, len(cid_list), 900):
        chunk = cid_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = f"SELECT {', '.join(select_cols)} FROM cid_mass WHERE {cid_col} IN ({placeholders})"
        rows.extend(cur.execute(q, chunk).fetchall())
    con.close()

    df = pd.DataFrame(rows, columns=out_cols).astype("string").drop_duplicates("PubChem_CID")
    for c in ["Molecular_Formula", "Exact_Mass"]:
        if c not in df.columns:
            df[c] = pd.NA
    return df[["PubChem_CID", "Molecular_Formula", "Exact_Mass"]]


def fetch_pubchem_class_sqlite(db_path: str, cids: Set[str]) -> pd.DataFrame:
    if not cids or not _sqlite_table_exists(db_path, "cid_classification"):
        return pd.DataFrame(columns=["PubChem_CID", "Super_Class", "Main_Class", "Sub_Class"], dtype="string")

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows: List[Tuple[str, str, str, str]] = []
    cid_list = list(cids)
    for i in range(0, len(cid_list), 900):
        chunk = cid_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = (
            f"SELECT cid, super_class, main_class, sub_class "
            f"FROM cid_classification WHERE cid IN ({placeholders})"
        )
        rows.extend(cur.execute(q, chunk).fetchall())
    con.close()

    return pd.DataFrame(rows, columns=["PubChem_CID", "Super_Class", "Main_Class", "Sub_Class"]).astype("string").drop_duplicates("PubChem_CID")


def hmdb_local_hits(hmdb_lite: pd.DataFrame, met_name: str) -> pd.DataFrame:
    met_lc = met_name.lower()

    def row_matches(row) -> bool:
        name_ok = _clean_scalar(row.get("hmdb_metabolite_name"))
        name_match = (name_ok is not None and name_ok.lower() == met_lc)
        syn_match = met_lc in [s.lower() for s in _split_synonyms(row.get("hmdb_metabolite_synonyms_synonym"))]
        return bool(name_match or syn_match)

    mask = hmdb_lite.apply(row_matches, axis=1)
    return hmdb_lite.loc[mask].copy()


def match_metabolite_ids(metabolite_names: Sequence[str], hmdb_lite: pd.DataFrame, sqlite_db: str, sleep: float, pugview_workers: int = 8) -> pd.DataFrame:
    out_rows = []
    for met in metabolite_names:
        hmdb_hits = hmdb_local_hits(hmdb_lite, met)
        hmdb_ids_local = sorted(set([x for x in hmdb_hits.get("hmdb_metabolite_accession", pd.Series([], dtype="string")).dropna().astype(str) if _clean_scalar(x)]))
        pubchem_local = sorted(set([_clean_cid(x) for x in hmdb_hits.get("hmdb_metabolite_pubchem_compound_id", pd.Series([], dtype="string")).dropna().astype(str)]))
        pubchem_local = [x for x in pubchem_local if x]

        pubchem_cids = fetch_cids_by_name_sqlite(sqlite_db, met)

        hmdb_from_pubchem = []
        if pubchem_cids:
            cid_ann = fetch_pugview_ids_for_cids(pubchem_cids, sleep=sleep, max_workers=pugview_workers)
            if not cid_ann.empty:
                for val in cid_ann["HMDB_ID"].astype("string").dropna().tolist():
                    hmdb_from_pubchem.extend(_split_ids(val))

        hmdb_from_pubchem = sorted(set([x for x in hmdb_from_pubchem if x]))

        has_pubchem_cross = bool(pubchem_local and pubchem_cids)
        pubchem_match = (len(set(pubchem_local).intersection(pubchem_cids)) > 0) if has_pubchem_cross else pd.NA

        has_hmdb_cross = bool(hmdb_ids_local and hmdb_from_pubchem)
        hmdb_match = (len(set(hmdb_ids_local).intersection(hmdb_from_pubchem)) > 0) if has_hmdb_cross else pd.NA

        all_na = not any([hmdb_ids_local, pubchem_local, pubchem_cids, hmdb_from_pubchem])
        any_mismatch = (has_pubchem_cross and not pubchem_match) or (has_hmdb_cross and not hmdb_match)

        if all_na:
            conf = "unmapped"
        elif any_mismatch:
            conf = "non_confidence"
        else:
            conf = "confidence"

        out_rows.append({
            "metabolite_query": met,
            "hmdb_id_from_hmdb": ";".join(hmdb_ids_local) if hmdb_ids_local else pd.NA,
            "pubchem_id_from_hmdb": ";".join(pubchem_local) if pubchem_local else pd.NA,
            "pubchem_id_from_pubchem": ";".join(pubchem_cids) if pubchem_cids else pd.NA,
            "hmdb_id_from_pubchem": ";".join(hmdb_from_pubchem) if hmdb_from_pubchem else pd.NA,
            "pubchem_id_match": pubchem_match,
            "hmdb_id_match": hmdb_match,
            "confidence_status": conf,
        })

    out = pd.DataFrame(out_rows)
    for c in ["hmdb_id_from_hmdb", "pubchem_id_from_hmdb", "pubchem_id_from_pubchem", "hmdb_id_from_pubchem"]:
        out[c] = out[c].astype("string")
    return out


def annotate_pubchem_offline(df: pd.DataFrame, sqlite_db: str, max_synonyms_per_cid: int = 200) -> pd.DataFrame:
    out = df.copy()
    if "PubChem_CID" in out.columns:
        out["PubChem_CID"] = out["PubChem_CID"].astype("string").map(_clean_cid).astype("string")

    for c in ["SMILES", "Synonyms", "InChIKey", "Standardized_Name", "Molecular_Formula", "Exact_Mass"]:
        if c not in out.columns:
            out[c] = pd.NA
        out[c] = _naify_string_series(out[c])

    cids = set([x for x in out["PubChem_CID"].astype("string").dropna().tolist()])
    core = fetch_smiles_sqlite(sqlite_db, cids)
    syn = fetch_synonyms_sqlite(sqlite_db, cids, max_synonyms_per_cid=max_synonyms_per_cid)
    inchi = fetch_inchikey_sqlite(sqlite_db, cids)
    title = fetch_pubchem_title_sqlite(sqlite_db, cids)
    iupac = fetch_pubchem_iupac_sqlite(sqlite_db, cids)
    mass = fetch_pubchem_mass_sqlite(sqlite_db, cids)

    out = out.merge(core, on="PubChem_CID", how="left", suffixes=("", "_pc"))
    out = out.merge(syn, on="PubChem_CID", how="left", suffixes=("", "_pc"))
    out = out.merge(inchi, on="PubChem_CID", how="left", suffixes=("", "_pc"))
    out = out.merge(title, on="PubChem_CID", how="left", suffixes=("", "_pc"))
    out = out.merge(iupac, on="PubChem_CID", how="left", suffixes=("", "_pc"))
    out = out.merge(mass, on="PubChem_CID", how="left", suffixes=("", "_pc"))

    if "SMILES_pc" in out.columns:
        out["SMILES"] = _naify_string_series(out["SMILES"]).fillna(_naify_string_series(out["SMILES_pc"]))
        out = out.drop(columns=["SMILES_pc"])
    if "Synonyms_pc" in out.columns:
        out["Synonyms"] = _naify_string_series(out["Synonyms"]).fillna(_naify_string_series(out["Synonyms_pc"]))
        out = out.drop(columns=["Synonyms_pc"])
    if "InChIKey_pc" in out.columns:
        out["InChIKey"] = _naify_string_series(out["InChIKey"]).map(_clean_inchikey).astype("string").fillna(
            _naify_string_series(out["InChIKey_pc"]).map(_clean_inchikey).astype("string")
        )
        out = out.drop(columns=["InChIKey_pc"])
    if "Standardized_Name_pc" in out.columns:
        out["Standardized_Name"] = _naify_string_series(out["Standardized_Name"]).fillna(_naify_string_series(out["Standardized_Name_pc"]))
        out = out.drop(columns=["Standardized_Name_pc"])
    if "IUPAC_Name_pc" in out.columns:
        out["Standardized_Name"] = _naify_string_series(out["Standardized_Name"]).fillna(_naify_string_series(out["IUPAC_Name_pc"]))
        out = out.drop(columns=["IUPAC_Name_pc"])
    if "Molecular_Formula_pc" in out.columns:
        out["Molecular_Formula"] = _naify_string_series(out["Molecular_Formula"]).fillna(_naify_string_series(out["Molecular_Formula_pc"]))
        out = out.drop(columns=["Molecular_Formula_pc"])
    if "Exact_Mass_pc" in out.columns:
        out["Exact_Mass"] = _naify_string_series(out["Exact_Mass"]).fillna(_naify_string_series(out["Exact_Mass_pc"]))
        out = out.drop(columns=["Exact_Mass_pc"])
    return out


def _fetch_pubchem_property_chunk(cids: List[str], timeout: int = 30) -> pd.DataFrame:
    """
    Fetch Title, MolecularFormula, ExactMass from PubChem PUG property API.
    Returns empty DataFrame on failures.
    """
    if not cids:
        return pd.DataFrame(columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass"], dtype="string")

    cid_str = ",".join(cids)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid_str}/property/Title,MolecularFormula,ExactMass/JSON"
    )
    try:
        res = requests.get(url, timeout=timeout)
        if res.status_code != 200:
            return pd.DataFrame(columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass"], dtype="string")
        payload = res.json()
    except Exception:
        return pd.DataFrame(columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass"], dtype="string")

    props = payload.get("PropertyTable", {}).get("Properties", [])
    rows: List[Dict[str, object]] = []
    for rec in props:
        cid = _clean_cid(rec.get("CID"))
        if not cid:
            continue
        rows.append({
            "PubChem_CID": cid,
            "Standardized_Name": rec.get("Title", pd.NA),
            "Molecular_Formula": rec.get("MolecularFormula", pd.NA),
            "Exact_Mass": str(rec.get("ExactMass")) if rec.get("ExactMass") is not None else pd.NA,
        })

    if not rows:
        return pd.DataFrame(columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass"], dtype="string")
    return pd.DataFrame(rows).astype("string").drop_duplicates("PubChem_CID")


def fetch_pubchem_props_api(cids: Sequence[str], max_workers: int = 8, chunk_size: int = 100) -> pd.DataFrame:
    cleaned = sorted({c for c in (_clean_cid(x) for x in cids) if c})
    if not cleaned:
        return pd.DataFrame(columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass"], dtype="string")

    chunks = [cleaned[i:i + chunk_size] for i in range(0, len(cleaned), chunk_size)]
    out_frames: List[pd.DataFrame] = []

    if max_workers <= 1 or len(chunks) == 1:
        for ch in chunks:
            out_frames.append(_fetch_pubchem_property_chunk(ch))
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs = [ex.submit(_fetch_pubchem_property_chunk, ch) for ch in chunks]
            for fut in as_completed(futs):
                try:
                    out_frames.append(fut.result())
                except Exception:
                    continue

    if not out_frames:
        return pd.DataFrame(columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass"], dtype="string")

    out = pd.concat(out_frames, ignore_index=True)
    if out.empty:
        return pd.DataFrame(columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass"], dtype="string")
    return out.astype("string").drop_duplicates("PubChem_CID")


def apply_pubchem_api_enrichment(df: pd.DataFrame, max_workers: int = 8) -> pd.DataFrame:
    out = df.copy()
    for c in ["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass"]:
        if c not in out.columns:
            out[c] = pd.NA

    out["PubChem_CID"] = out["PubChem_CID"].astype("string").map(_clean_cid).astype("string")
    api_df = fetch_pubchem_props_api(out["PubChem_CID"].astype("string").dropna().tolist(), max_workers=max_workers)
    if api_df.empty:
        return out

    out = out.merge(api_df, on="PubChem_CID", how="left", suffixes=("", "_api"))
    for c in ["Standardized_Name", "Molecular_Formula", "Exact_Mass"]:
        api_col = f"{c}_api"
        if api_col in out.columns:
            out[c] = _naify_string_series(out[c]).fillna(_naify_string_series(out[api_col]))
            out = out.drop(columns=[api_col])
    return out


def fill_hmdb_columns(df: pd.DataFrame, hmdb_lite: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    hmdb_map = hmdb_lite.copy()
    hmdb_map["hmdb_metabolite_accession"] = (
        hmdb_map["hmdb_metabolite_accession"].astype("string").str.strip().str.upper()
    )
    hmdb_map = hmdb_map.dropna(subset=["hmdb_metabolite_accession"])
    hmdb_map = hmdb_map[hmdb_map["hmdb_metabolite_accession"] != ""]
    hmdb_map = hmdb_map.drop_duplicates(subset=["hmdb_metabolite_accession"], keep="first")

    idx = out["HMDB_ID"].map(_first_hmdb_id)

    def fill(col: str, hmdb_col: str):
        lookup = pd.Series(
            hmdb_map[hmdb_col].astype("string").values,
            index=hmdb_map["hmdb_metabolite_accession"]
        ).to_dict()
        vals = idx.map(lambda v: lookup.get(v, pd.NA))
        out[col] = _naify_string_series(out[col]).fillna(pd.Series(vals, dtype="string"))

    fill("Standardized_Name", "hmdb_metabolite_name")
    fill("Synonyms", "hmdb_metabolite_synonyms_synonym")
    fill("Molecular_Formula", "hmdb_metabolite_chemical_formula")
    fill("Exact_Mass", "hmdb_metabolite_monisotopic_molecular_weight")
    fill("InChIKey", "hmdb_metabolite_inchikey")
    fill("SMILES", "hmdb_metabolite_smiles")
    fill("KEGG_ID", "hmdb_metabolite_kegg_id")
    fill("ChEBI_ID", "hmdb_metabolite_chebi_id")
    fill("Super_Class", "hmdb_metabolite_taxonomy_super_class")
    fill("Main_Class", "hmdb_metabolite_taxonomy_class")
    fill("Sub_Class", "hmdb_metabolite_taxonomy_sub_class")
    return out


def apply_class_fallback_from_pubchem_db(df: pd.DataFrame, sqlite_db: str) -> pd.DataFrame:
    out = df.copy()
    for c in ["PubChem_CID", "Super_Class", "Main_Class", "Sub_Class"]:
        if c not in out.columns:
            out[c] = pd.NA

    cids = set([x for x in out["PubChem_CID"].astype("string").map(_clean_cid).dropna().tolist()])
    class_df = fetch_pubchem_class_sqlite(sqlite_db, cids)
    if class_df.empty:
        return out

    out = out.merge(class_df, on="PubChem_CID", how="left", suffixes=("", "_pcclass"))
    for c in ["Super_Class", "Main_Class", "Sub_Class"]:
        cc = f"{c}_pcclass"
        if cc in out.columns:
            out[c] = _naify_string_series(out[c]).fillna(_naify_string_series(out[cc]))
            out = out.drop(columns=[cc])
    return out


def _has_valid_smiles(smiles: object) -> bool:
    s = _clean_scalar(smiles)
    if not s:
        return False
    if not HAS_RDKIT:
        # fallback mode: any non-empty smiles-like string is considered usable
        return True
    return Chem.MolFromSmiles(str(s)) is not None


def _smiles_shingles(smiles: str, k: int = 2) -> Set[str]:
    """Lightweight fallback fingerprint when RDKit is unavailable."""
    s = _clean_scalar(smiles)
    if not s:
        return set()
    s = str(s)
    if len(s) < k:
        return {s}
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def _set_tanimoto(a: Set[str], b: Set[str]) -> Optional[float]:
    if not a or not b:
        return None
    inter = len(a.intersection(b))
    union = len(a.union(b))
    if union == 0:
        return None
    return float(inter / union)


def tanimoto_scores(smiles_values: Sequence[str], ref_idx: Optional[int]) -> List[Optional[float]]:
    if ref_idx is None or ref_idx < 0 or ref_idx >= len(smiles_values):
        return [pd.NA] * len(smiles_values)

    if HAS_RDKIT:
        mols = [Chem.MolFromSmiles(s) if _clean_scalar(s) else None for s in smiles_values]
        if mols[ref_idx] is None:
            return [pd.NA] * len(smiles_values)

        fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048) if m is not None else None for m in mols]
        ref_fp = fps[ref_idx]
        if ref_fp is None:
            return [pd.NA] * len(smiles_values)

        out = []
        for fp in fps:
            if fp is None:
                out.append(pd.NA)
            else:
                out.append(float(DataStructs.TanimotoSimilarity(fp, ref_fp)))
        return out

    # RDKit-unavailable fallback: q-gram set Tanimoto over SMILES strings.
    ref_set = _smiles_shingles(smiles_values[ref_idx])
    if not ref_set:
        return [pd.NA] * len(smiles_values)

    out: List[Optional[float]] = []
    for s in smiles_values:
        cur = _smiles_shingles(s)
        score = _set_tanimoto(cur, ref_set)
        out.append(pd.NA if score is None else score)
    return out


def ensure_cols(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    out = df.copy()
    for c in columns:
        if c not in out.columns:
            out[c] = pd.NA
    return out

PUBLIC_STAGE_COLUMNS = [
    "HuMANet_ID", "original_query_name", "Query_Name", "Database_Source", "InChIKey", "HMDB_ID", "PubChem_CID", "KEGG_ID", "ChEBI_ID",
    "Molecular_Formula", "Exact_Mass", "Super_Class", "Main_Class", "Sub_Class", "Standardized_Name", "SMILES",
    "Synonyms", "matched_name", "Study_Folder",
]
PUBLIC_CONFIDENCE_COLUMNS = PUBLIC_STAGE_COLUMNS[:-1] + ["Confidence", "Study_Folder"]
PUBLIC_NONCONF_COLUMNS = PUBLIC_STAGE_COLUMNS[:-1] + ["Tanimoto_vs_Reference", "Confidence", "Study_Folder"]


def _finalize_public_stage_output(df: pd.DataFrame, columns) -> pd.DataFrame:
    out = df.copy()
    if "Study_Folder" not in out.columns:
        if "Input_File" in out.columns:
            out["Study_Folder"] = out["Input_File"].astype("string")
        else:
            out["Study_Folder"] = pd.Series([pd.NA] * len(out), dtype="string")
    else:
        out["Study_Folder"] = out["Study_Folder"].astype("string")
    for col in columns:
        if col not in out.columns:
            out[col] = pd.NA
    return out[columns].copy()


def run_stage2(stage1_unmapped_csv: str, hmdb_lite_csv: str, sqlite_db: str, out_prefix: str,
               r_classyfire_bridge: Optional[str], sleep: float = 0.3, pugview_workers: int = 8,
               skip_classyfire: bool = False, print_classyfire_output: bool = False, skip_stage5_excel: bool = False) -> Dict[str, pd.DataFrame]:
    total_start = time.perf_counter()

    with timed_step("Load stage1 unmapped + HMDB lite"):
        base = pd.read_csv(stage1_unmapped_csv, dtype="string")
        if "original_query_name" not in base.columns:
            base["original_query_name"] = base["Query_Name"].astype("string")
        hmdb_lite = pd.read_csv(hmdb_lite_csv, dtype="string")

    with timed_step("Initial metabolite ID matching"):
        out = match_metabolite_ids(base["Query_Name"].astype(str).tolist(), hmdb_lite, sqlite_db, sleep, pugview_workers=pugview_workers)

    unmapped_out = out[out["confidence_status"] == "unmapped"].copy()
    if not unmapped_out.empty:
        cleaned_names = unmapped_out["metabolite_query"].str.replace(r"\*$", "", regex=True)
        with timed_step("Re-match cleaned unmapped metabolite names"):
            out_stage1b = match_metabolite_ids(cleaned_names.tolist(), hmdb_lite, sqlite_db, sleep, pugview_workers=pugview_workers)

        map_b = {k: v for k, v in zip(out_stage1b["metabolite_query"], out_stage1b.to_dict("records"))}
        cols_to_replace = [
            "hmdb_id_from_hmdb", "pubchem_id_from_hmdb", "pubchem_id_from_pubchem",
            "hmdb_id_from_pubchem", "pubchem_id_match", "hmdb_id_match", "confidence_status"
        ]
        for idx, row in unmapped_out.iterrows():
            key = re.sub(r"\*$", "", row["metabolite_query"])
            if key in map_b:
                for c in cols_to_replace:
                    out.at[idx, c] = map_b[key][c]

    id_cols = ["hmdb_id_from_hmdb", "pubchem_id_from_hmdb", "pubchem_id_from_pubchem", "hmdb_id_from_pubchem"]
    for c in id_cols:
        out[c] = out[c].astype("string").replace({"NA": pd.NA, "": pd.NA, "Na": pd.NA, "null": pd.NA})

    multiple_ids = out[id_cols].fillna("").apply(lambda r: any(";" in x for x in r.values.tolist()), axis=1)
    out.loc[(out["confidence_status"] == "confidence") & multiple_ids, "confidence_status"] = "non_confidence"

    confidence_mets = set(out.loc[out["confidence_status"] == "confidence", "metabolite_query"].tolist())
    nonconfidence_mets = set(out.loc[out["confidence_status"] == "non_confidence", "metabolite_query"].tolist())
    unmapped_mets = set(out.loc[out["confidence_status"].isin(["unmapped"]) | out["confidence_status"].isna(), "metabolite_query"].tolist())

    humanet_confidence = base[base["Query_Name"].isin(confidence_mets)].copy()
    humanet_nonconfidence = base[base["Query_Name"].isin(nonconfidence_mets)].copy()
    humanet_unmapped = base[base["Query_Name"].isin(unmapped_mets)].copy()

    out_lookup = out.set_index("metabolite_query", drop=False)

    def _lookup_first_row(query: str):
        if query not in out_lookup.index:
            return None
        row = out_lookup.loc[query]
        if isinstance(row, pd.DataFrame):
            if row.empty:
                return None
            return row.iloc[0]
        return row

    source_map = {}
    for _, rr in out.iterrows():
        sources = []
        if _clean_scalar(rr.get("hmdb_id_from_hmdb")) or _clean_scalar(rr.get("pubchem_id_from_hmdb")):
            sources.append("HMDB")
        if _clean_scalar(rr.get("pubchem_id_from_pubchem")) or _clean_scalar(rr.get("hmdb_id_from_pubchem")):
            sources.append("PubChem")
        source_map[str(rr["metabolite_query"])] = "|".join(sources) if sources else pd.NA

    def pick_id(query: str, first: str, second: str) -> Optional[str]:
        row = _lookup_first_row(query)
        if row is None:
            return None
        a = _clean_scalar(row[first])
        b = _clean_scalar(row[second])
        return a or b

    humanet_confidence["HMDB_ID"] = humanet_confidence["HMDB_ID"].astype("string").map(_normalize_hmdb_field).astype("string")
    humanet_confidence["PubChem_CID"] = humanet_confidence["PubChem_CID"].astype("string").map(_clean_cid).astype("string")

    humanet_confidence["HMDB_ID"] = humanet_confidence.apply(
        lambda r: r["HMDB_ID"] if _clean_scalar(r["HMDB_ID"]) else _normalize_hmdb_field(pick_id(r["Query_Name"], "hmdb_id_from_hmdb", "hmdb_id_from_pubchem")),
        axis=1,
    )
    humanet_confidence["PubChem_CID"] = humanet_confidence.apply(
        lambda r: r["PubChem_CID"] if _clean_cid(r["PubChem_CID"]) else _clean_cid(pick_id(r["Query_Name"], "pubchem_id_from_pubchem", "pubchem_id_from_hmdb")),
        axis=1,
    )
    humanet_confidence["Database_Source"] = humanet_confidence["Query_Name"].map(source_map).astype("string")
    humanet_confidence["matched_name"] = humanet_confidence["Query_Name"].astype("string")

    for tbl in [humanet_confidence, humanet_nonconfidence, humanet_unmapped]:
        for c in ["HMDB_ID", "PubChem_CID", "KEGG_ID", "ChEBI_ID", "InChIKey", "SMILES", "Synonyms", "Standardized_Name", "Molecular_Formula", "Exact_Mass", "Super_Class", "Main_Class", "Sub_Class"]:
            if c not in tbl.columns:
                tbl[c] = pd.NA
            else:
                tbl[c] = _naify_string_series(tbl[c])

    with timed_step("Confidence table: offline PubChem annotation"):
        humanet_confidence = annotate_pubchem_offline(humanet_confidence, sqlite_db)

    with timed_step("Confidence table: PubChem property API enrichment"):
        humanet_confidence = apply_pubchem_api_enrichment(humanet_confidence, max_workers=pugview_workers)

    with timed_step("Confidence table: PubChem PUG-View ID backfill"):
        pubchem_ann = fetch_pugview_ids_for_cids(
            humanet_confidence["PubChem_CID"].astype("string").dropna().unique().tolist(),
            sleep=sleep,
            max_workers=pugview_workers,
        )
        if not pubchem_ann.empty:
            humanet_confidence = humanet_confidence.merge(pubchem_ann, on="PubChem_CID", how="left", suffixes=("", "_pv"))
            for c in ["HMDB_ID", "KEGG_ID", "ChEBI_ID"]:
                pv = f"{c}_pv"
                if pv in humanet_confidence.columns:
                    humanet_confidence[c] = _naify_string_series(humanet_confidence[c]).fillna(_naify_string_series(humanet_confidence[pv]))
                    humanet_confidence = humanet_confidence.drop(columns=[pv])

    with timed_step("Confidence table: HMDB fallback fills"):
        humanet_confidence = fill_hmdb_columns(humanet_confidence, hmdb_lite)

    if not skip_classyfire:
        with timed_step("Confidence table: ClassyFire annotation"):
            classy = annotate_classyfire(
                humanet_confidence["InChIKey"].astype("string").tolist(),
                r_classyfire_bridge,
                print_output=print_classyfire_output,
            )
            if not classy.empty:
                humanet_confidence["InChIKey"] = humanet_confidence["InChIKey"].astype("string").map(_clean_inchikey).astype("string")
                classy["InChIKey"] = classy["InChIKey"].astype("string").map(_clean_inchikey).astype("string")
                humanet_confidence = humanet_confidence.merge(classy, on="InChIKey", how="left", suffixes=("", "_cf"))
                for c in ["Super_Class", "Main_Class", "Sub_Class"]:
                    cc = f"{c}_cf"
                    if cc in humanet_confidence.columns:
                        humanet_confidence[c] = _naify_string_series(humanet_confidence[c])
                        humanet_confidence[c] = humanet_confidence[c].fillna(_naify_string_series(humanet_confidence[cc]))
                        humanet_confidence = humanet_confidence.drop(columns=[cc])

    with timed_step("Confidence table: class fallback from PubChem local DB"):
        humanet_confidence = apply_class_fallback_from_pubchem_db(humanet_confidence, sqlite_db)

    humanet_confidence["HMDB_ID"] = humanet_confidence["HMDB_ID"].astype("string").map(_first_hmdb_id).astype("string")
    humanet_confidence["InChIKey"] = humanet_confidence["InChIKey"].astype("string").map(_clean_inchikey).astype("string")
    humanet_confidence["ChEBI_ID"] = humanet_confidence["ChEBI_ID"].astype("string").map(_normalize_chebi_field).astype("string")
    humanet_confidence["Confidence"] = "True"
    humanet_unmapped["Database_Source"] = pd.NA
    humanet_unmapped["matched_name"] = pd.NA
    if "matched_name" not in humanet_unmapped.columns:
        humanet_unmapped["matched_name"] = pd.NA

    with timed_step("Build non-confidence expanded rows"):
        expanded = []
        for _, row in humanet_nonconfidence.iterrows():
            q = row["Query_Name"]
            rr = _lookup_first_row(q)
            if rr is None:
                expanded.append(row.to_dict())
                continue

            hmdb_ids = _split_ids(rr["hmdb_id_from_hmdb"])
            pubchem_ids_hmdb = _split_ids(rr["pubchem_id_from_hmdb"])
            pubchem_ids_pc = _split_ids(rr["pubchem_id_from_pubchem"])
            hmdb_ids_pc = _split_ids(rr["hmdb_id_from_pubchem"])

            for i, hid in enumerate(hmdb_ids):
                r = row.to_dict()
                r["Database_Source"] = "HMDB"
                r["HMDB_ID"] = hid
                r["PubChem_CID"] = pubchem_ids_hmdb[i] if i < len(pubchem_ids_hmdb) else pd.NA
                expanded.append(r)

            for i, cid in enumerate(pubchem_ids_pc):
                r = row.to_dict()
                r["Database_Source"] = "PubChem"
                r["PubChem_CID"] = cid
                r["HMDB_ID"] = hmdb_ids_pc[i] if i < len(hmdb_ids_pc) else pd.NA
                expanded.append(r)

            if not hmdb_ids and not pubchem_ids_pc:
                expanded.append(row.to_dict())

        nonconf_expanded = pd.DataFrame(expanded) if expanded else humanet_nonconfidence.copy()
        nonconf_expanded = ensure_cols(nonconf_expanded, list(base.columns) + ["Database_Source", "Tanimoto_vs_Reference", "matched_name"])
        nonconf_expanded["matched_name"] = nonconf_expanded["Query_Name"].astype("string")
        nonconf_expanded["matched_name"] = nonconf_expanded["Query_Name"].astype("string")

    with timed_step("Non-confidence HMDB enrichment"):
        hmdb_rows = nonconf_expanded["Database_Source"] == "HMDB"
        if hmdb_rows.any():
            hmdb_part = fill_hmdb_columns(nonconf_expanded.loc[hmdb_rows].copy(), hmdb_lite)
            nonconf_expanded.loc[hmdb_rows, hmdb_part.columns] = hmdb_part.values

    with timed_step("Non-confidence PubChem enrichment"):
        pubchem_rows = nonconf_expanded["Database_Source"] == "PubChem"
        if pubchem_rows.any():
            pc_part = annotate_pubchem_offline(nonconf_expanded.loc[pubchem_rows].copy(), sqlite_db)
            pc_part = apply_pubchem_api_enrichment(pc_part, max_workers=pugview_workers)
            ann = fetch_pugview_ids_for_cids(
                pc_part["PubChem_CID"].astype("string").dropna().unique().tolist(),
                sleep=sleep,
                max_workers=pugview_workers,
            )
            if not ann.empty:
                pc_part = pc_part.merge(ann, on="PubChem_CID", how="left", suffixes=("", "_pv"))
                for c in ["HMDB_ID", "KEGG_ID", "ChEBI_ID"]:
                    pv = f"{c}_pv"
                    if pv in pc_part.columns:
                        pc_part[c] = _naify_string_series(pc_part[c]).fillna(_naify_string_series(pc_part[pv]))
                        pc_part = pc_part.drop(columns=[pv])

            if not skip_classyfire:
                cf = annotate_classyfire(
                    pc_part["InChIKey"].astype("string").tolist(),
                    r_classyfire_bridge,
                    print_output=print_classyfire_output,
                )
                if not cf.empty:
                    pc_part["InChIKey"] = pc_part["InChIKey"].astype("string").map(_clean_inchikey).astype("string")
                    cf["InChIKey"] = cf["InChIKey"].astype("string").map(_clean_inchikey).astype("string")
                    pc_part = pc_part.merge(cf, on="InChIKey", how="left", suffixes=("", "_cf"))
                    for c in ["Super_Class", "Main_Class", "Sub_Class"]:
                        cc = f"{c}_cf"
                        if cc in pc_part.columns:
                            pc_part[c] = _naify_string_series(pc_part[c])
                            pc_part[c] = pc_part[c].fillna(_naify_string_series(pc_part[cc]))
                            pc_part = pc_part.drop(columns=[cc])

            pc_part = apply_class_fallback_from_pubchem_db(pc_part, sqlite_db)
            nonconf_expanded.loc[pubchem_rows, pc_part.columns] = pc_part.values

    with timed_step("Non-confidence HMDB fallback fills"):
        # Apply HMDB fallback to all rows (including PubChem rows) for any fields still missing.
        # This mirrors the confidence flow where HMDB can backfill name/formula/mass after PubChem lookup.
        nonconf_expanded = fill_hmdb_columns(nonconf_expanded, hmdb_lite)

    nonconf_expanded["HMDB_ID"] = nonconf_expanded["HMDB_ID"].astype("string").map(_first_hmdb_id).astype("string")
    nonconf_expanded["InChIKey"] = nonconf_expanded["InChIKey"].astype("string").map(_clean_inchikey).astype("string")
    nonconf_expanded["ChEBI_ID"] = nonconf_expanded["ChEBI_ID"].astype("string").map(_normalize_chebi_field).astype("string")
    nonconf_expanded["Confidence"] = "False"

    with timed_step("Non-confidence Tanimoto scoring"):
        if "HuMANet_ID" in nonconf_expanded.columns:
            scores = []
            for _, grp in nonconf_expanded.groupby("HuMANet_ID", dropna=False):
                ref_idx = None
                hmdb_valid = grp.index[(grp["Database_Source"] == "HMDB") & grp["SMILES"].map(_has_valid_smiles)]
                if len(hmdb_valid) > 0:
                    ref_idx = grp.index.get_loc(hmdb_valid[0])
                else:
                    any_valid = grp.index[grp["SMILES"].map(_has_valid_smiles)]
                    if len(any_valid) > 0:
                        ref_idx = grp.index.get_loc(any_valid[0])
                vals = tanimoto_scores(grp["SMILES"].tolist(), ref_idx)
                scores.extend(list(zip(grp.index.tolist(), vals)))
            for ix, sc in scores:
                nonconf_expanded.at[ix, "Tanimoto_vs_Reference"] = sc

    humanet_confidence = _finalize_public_stage_output(humanet_confidence, PUBLIC_CONFIDENCE_COLUMNS)
    nonconf_expanded = _finalize_public_stage_output(nonconf_expanded, PUBLIC_NONCONF_COLUMNS)
    humanet_unmapped = _finalize_public_stage_output(humanet_unmapped, PUBLIC_STAGE_COLUMNS)

    outputs = {
        "confidence": humanet_confidence,
        "nonconfidence_expanded": nonconf_expanded,
        "unmapped": humanet_unmapped,
        "id_tracking_table": out,
    }

    with timed_step("Write output CSV files"):
        for name, df in outputs.items():
            df.to_csv(f"{out_prefix}_{name}.csv", index=False)

    excel_sheets = {
        "confidence": humanet_confidence,
        "nonconfidence": nonconf_expanded,
        "unmapped": humanet_unmapped,
    }
    if not skip_stage5_excel:
        stage5_sheets = build_stage5_annotation_sheets(
            pd.concat([humanet_confidence, nonconf_expanded], ignore_index=True, sort=False)
        )
        excel_sheets.update(stage5_sheets)

    with timed_step("Write Stage2 Excel workbook"):
        write_excel_workbook(excel_sheets, f"{out_prefix}.xlsx")

    print(f"[TIMER] TOTAL run_stage2: {time.perf_counter() - total_start:.2f}s")
    return outputs


def main():
    parser = argparse.ArgumentParser(description="HuMANet Stage2 Python pipeline (hybrid offline + optional R classyfire bridge).")
    parser.add_argument("--stage1_unmapped", default="Refmet_unmapped.csv", help="Stage1 unmapped CSV")
    parser.add_argument("--out_prefix", default="stage2", help="Output file prefix")
    parser.add_argument("--sleep", type=float, default=0.3, help="Sleep seconds between PubChem PUG-View requests")
    parser.add_argument("--pugview_workers", type=int, default=8, help="Thread workers for PubChem PUG-View lookups")
    parser.add_argument("--skip_classyfire", action="store_true", help="Disable ClassyFire bridge and use only local fallback for class fields")
    parser.add_argument("--print_classyfire_output", action="store_true", help="Print ClassyFire result rows to stdout for debugging")
    parser.add_argument("--skip_stage5_excel", action="store_true", help="Skip adding Stage5 sheets to the Excel workbook")
    args = parser.parse_args()

    bridge = CLASSYFIRE_BRIDGE_R if HAS_RPY2 else None
    run_stage2(
        stage1_unmapped_csv=args.stage1_unmapped,
        hmdb_lite_csv=HMDB_LITE_CSV,
        sqlite_db=PUBCHEM_OFFLINE_SQLITE,
        out_prefix=args.out_prefix,
        r_classyfire_bridge=bridge,
        sleep=args.sleep,
        pugview_workers=args.pugview_workers,
        skip_classyfire=args.skip_classyfire,
        print_classyfire_output=args.print_classyfire_output,
        skip_stage5_excel=args.skip_stage5_excel,
    )


if __name__ == "__main__":
    main()
