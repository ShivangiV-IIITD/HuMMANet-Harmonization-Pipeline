import argparse
import re
import sqlite3
import time
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd

import r_environment_bootstrap  # noqa: F401
from knowledge_annotation_step import build_stage5_annotation_sheets
from pipeline_utils import (
    HAS_RPY2,
    _clean_cid,
    _clean_hmdb_id,
    _clean_inchikey,
    _clean_scalar,
    _first_hmdb_id,
    _normalize_hmdb_field,
    _split_ids,
    _split_synonyms,
    _sqlite_table_exists,
    extract_ids_from_pugview,
    fetch_cids_by_name_sqlite,
    fetch_pugview_ids_for_cids,
    map_refmet_queries,
    prepare_stage1_input_dataframe,
    timed_step,
    write_excel_workbook,
)
from resource_config import HMDB_LITE_CSV, PUBCHEM_OFFLINE_SQLITE, REFMET_BRIDGE_R

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except Exception:
    HAS_RDKIT = False


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

    return pd.DataFrame(rows, columns=["PubChem_CID", "SMILES"]).astype("string")


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

    return pd.DataFrame(rows, columns=["PubChem_CID", "InChIKey"]).astype("string")


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


def fetch_pubchem_class_sqlite(db_path: str, cids: Set[str]) -> pd.DataFrame:
    """
    Optional offline class fallback if db contains table `cid_classification` with columns:
      cid, super_class, main_class, sub_class
    """
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

    return pd.DataFrame(rows, columns=["PubChem_CID", "Super_Class", "Main_Class", "Sub_Class"]).astype("string")


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

        if all_na or any_mismatch:
            conf = "unmapped"
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
    cids = set([x for x in out["PubChem_CID"].astype("string").dropna().tolist()])
    core = fetch_smiles_sqlite(sqlite_db, cids)
    syn = fetch_synonyms_sqlite(sqlite_db, cids, max_synonyms_per_cid=max_synonyms_per_cid)
    inchi = fetch_inchikey_sqlite(sqlite_db, cids)
    out = out.merge(core, on="PubChem_CID", how="left", suffixes=("", "_pc"))
    out = out.merge(syn, on="PubChem_CID", how="left", suffixes=("", "_pc"))
    out = out.merge(inchi, on="PubChem_CID", how="left", suffixes=("", "_pc"))

    if "SMILES_pc" in out.columns:
        out["SMILES"] = out["SMILES"].astype("string").fillna(out["SMILES_pc"].astype("string"))
        out = out.drop(columns=["SMILES_pc"])
    if "Synonyms_pc" in out.columns:
        out["Synonyms"] = out["Synonyms"].astype("string").fillna(out["Synonyms_pc"].astype("string"))
        out = out.drop(columns=["Synonyms_pc"])
    if "InChIKey_pc" in out.columns:
        out["InChIKey"] = out["InChIKey"].astype("string").fillna(out["InChIKey_pc"].astype("string"))
        out = out.drop(columns=["InChIKey_pc"])
    return out


def fill_hmdb_columns(df: pd.DataFrame, hmdb_lite: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    hmdb_map = hmdb_lite.dropna(subset=["hmdb_metabolite_accession"]).drop_duplicates(
        subset=["hmdb_metabolite_accession"], keep="first"
    )
    idx = out["HMDB_ID"].map(_first_hmdb_id)
    map_df = hmdb_map.set_index("hmdb_metabolite_accession", drop=False)

    def fill(col: str, hmdb_col: str):
        vals = idx.map(lambda v: map_df.at[v, hmdb_col] if v in map_df.index else pd.NA)
        out[col] = out[col].astype("string").fillna(pd.Series(vals, dtype="string"))

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
            out[c] = out[c].astype("string").fillna(out[cc].astype("string"))
            out = out.drop(columns=[cc])
    return out


def _has_valid_smiles(smiles: object) -> bool:
    if not _clean_scalar(smiles):
        return False
    if not HAS_RDKIT:
        return True
    return Chem.MolFromSmiles(str(smiles)) is not None


def _smiles_shingles(smiles: object, k: int = 3) -> Set[str]:
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


def _naify_stage1_series(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .replace({"": pd.NA, "-": pd.NA, "NA": pd.NA, "Na": pd.NA, "null": pd.NA, "NA_character_": pd.NA, "<NA>": pd.NA})
    )


def smart_split_metabolite(x: object) -> List[str]:
    s = _clean_scalar(x)
    if not s:
        return []

    parts, buf = [], []
    depth_sq = 0
    depth_par = 0
    for ch in str(s).strip():
        if ch == "[":
            depth_sq += 1
        elif ch == "]":
            depth_sq = max(0, depth_sq - 1)
        elif ch == "(":
            depth_par += 1
        elif ch == ")":
            depth_par = max(0, depth_par - 1)

        if ch == "/" and depth_sq == 0 and depth_par == 0:
            part = "".join(buf).strip()
            if part:
                parts.append(part)
            buf = []
        else:
            buf.append(ch)

    tail = "".join(buf).strip()
    if tail:
        parts.append(tail)

    out = []
    seen = set()
    for part in parts:
        cleaned = re.sub(r"\s+", " ", part.strip())
        if cleaned and cleaned not in seen:
            seen.add(cleaned)
            out.append(cleaned)
    return out


def _first_existing(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    cols = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        hit = cols.get(candidate.lower())
        if hit is not None:
            return hit
    return None


def _join_input_rows_with_refmet(base: pd.DataFrame, refmet_df: pd.DataFrame) -> pd.DataFrame:
    left = base.copy()
    left["Database_Source"] = "RefMet"

    if refmet_df.empty:
        return left

    work = refmet_df.copy()
    input_name_col = _first_existing(work, ["Input.name", "Input_name", "input.name"])
    if input_name_col is None:
        raise ValueError("RefMet bridge output is missing the input-name column.")

    keep_first = work.drop_duplicates(subset=[input_name_col], keep="first").copy()
    merged = left.merge(keep_first, left_on="Query_Name", right_on=input_name_col, how="left")

    if input_name_col in merged.columns:
        merged = merged.drop(columns=[input_name_col])

    def coalesce_into_target(target: str, candidates: Sequence[str]) -> None:
        present = [col for col in candidates if col in merged.columns]
        if not present:
            return
        collapsed = merged[present[0]].astype("string")
        for col in present[1:]:
            collapsed = collapsed.fillna(merged[col].astype("string"))
        for col in present:
            if col in merged.columns:
                merged.drop(columns=[col], inplace=True)
        merged[target] = collapsed

    coalesce_into_target("InChIKey", ["InChIKey", "INCHI_KEY", "InChIKey_x", "InChIKey_y"])
    coalesce_into_target("HMDB_ID", ["HMDB_ID", "HMDB_ID_x", "HMDB_ID_y"])
    coalesce_into_target("PubChem_CID", ["PubChem_CID", "PubChem CID", "PubChem_CID_x", "PubChem_CID_y"])
    coalesce_into_target("KEGG_ID", ["KEGG_ID", "Kegg.Id", "Kegg ID", "KEGG_ID_x", "KEGG_ID_y"])
    coalesce_into_target("ChEBI_ID", ["ChEBI_ID", "CHEBI_ID", "CheBI_ID", "ChEBI_ID_x", "ChEBI_ID_y"])
    coalesce_into_target("Molecular_Formula", ["Molecular_Formula", "Molecular.Formula", "Formula", "Molecular_Formula_x", "Molecular_Formula_y"])
    coalesce_into_target("Exact_Mass", ["Exact_Mass", "Exact.mass", "Exact.Mass", "Exact_Mass_x", "Exact_Mass_y"])
    coalesce_into_target("Super_Class", ["Super_Class", "Super.class", "Super_Class_x", "Super_Class_y"])
    coalesce_into_target("Main_Class", ["Main_Class", "Main.class", "Main_Class_x", "Main_Class_y"])
    coalesce_into_target("Sub_Class", ["Sub_Class", "Sub.class", "Sub_Class_x", "Sub_Class_y"])
    coalesce_into_target("Standardized_Name", ["Standardized_Name", "Standardized.name", "Standardized.Name", "Metabolite_Name", "Standardized_Name_x", "Standardized_Name_y"])
    coalesce_into_target("SMILES", ["SMILES", "SMILES_x", "SMILES_y"])
    coalesce_into_target("Synonyms", ["Synonyms", "Synonyms_x", "Synonyms_y"])

    if "LM_ID" in merged.columns:
        merged = merged.drop(columns=["LM_ID"])
    if "RefMet_ID" in merged.columns:
        merged = merged.drop(columns=["RefMet_ID"])

    if "Standardized_Name" not in merged.columns:
        std_col = _first_existing(merged, ["Standardized.Name", "Standardized_Name", "Metabolite_Name"])
        if std_col:
            merged = merged.rename(columns={std_col: "Standardized_Name"})

    dedupe_targets = [
        "InChIKey",
        "HMDB_ID",
        "PubChem_CID",
        "KEGG_ID",
        "ChEBI_ID",
        "Molecular_Formula",
        "Exact_Mass",
        "Super_Class",
        "Main_Class",
        "Sub_Class",
        "Standardized_Name",
        "SMILES",
        "Synonyms",
    ]
    for target in dedupe_targets:
        dupes = merged.loc[:, merged.columns == target]
        if dupes.shape[1] > 1:
            collapsed = dupes.iloc[:, 0].astype("string")
            for j in range(1, dupes.shape[1]):
                collapsed = collapsed.fillna(dupes.iloc[:, j].astype("string"))
            merged = merged.drop(columns=target)
            merged[target] = collapsed

    return merged


PUBLIC_STAGE_COLUMNS = [
    "HuMANet_ID", "original_query_name", "Query_Name", "Database_Source", "InChIKey", "HMDB_ID", "PubChem_CID", "KEGG_ID", "ChEBI_ID",
    "Molecular_Formula", "Exact_Mass", "Super_Class", "Main_Class", "Sub_Class", "Standardized_Name", "SMILES",
    "Synonyms", "matched_name", "Study_Folder",
]


def _finalize_public_stage_output(df: pd.DataFrame, columns=PUBLIC_STAGE_COLUMNS) -> pd.DataFrame:
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


def run_stage1(stage1_unmapped_csv: str, hmdb_lite_csv: str, sqlite_db: str, out_prefix: str,
               r_classyfire_bridge: Optional[str], sleep: float = 0.3, pugview_workers: int = 8,
               skip_classyfire: bool = False, skip_stage5_excel: bool = False) -> Dict[str, pd.DataFrame]:
    total_start = time.perf_counter()

    with timed_step("Load stage1 input + HMDB lite"):
        base = pd.read_csv(stage1_unmapped_csv, dtype="string")
        base = prepare_stage1_input_dataframe(base, stage1_unmapped_csv)
        if "original_query_name" not in base.columns:
            base["original_query_name"] = base["Query_Name"].astype("string")
        hmdb_lite = pd.read_csv(hmdb_lite_csv, dtype="string")

    with timed_step("Split slash-delimited metabolite names"):
        expanded_rows = []
        for _, row in base.iterrows():
            original_query = row["original_query_name"] if "original_query_name" in row.index else row["Query_Name"]
            parts = smart_split_metabolite(row["Query_Name"])
            if not parts:
                parts = [row["Query_Name"]]
            for part in parts:
                record = row.to_dict()
                record["original_query_name"] = original_query
                record["Query_Name"] = part
                expanded_rows.append(record)
        query_table = pd.DataFrame(expanded_rows).astype("string")
        query_table = query_table.dropna(subset=["Query_Name"])
        query_table = query_table[query_table["Query_Name"].astype("string").str.strip().ne("")]
        print(f"[INFO] Stage1 query rows after split: {len(query_table)}")

    with timed_step("Map metabolites via RefMet bridge"):
        refmet_mapped = map_refmet_queries(query_table["Query_Name"].tolist(), REFMET_BRIDGE_R if HAS_RPY2 else None)
        print(f"[INFO] RefMet rows returned: {len(refmet_mapped)}")

    with timed_step("Join RefMet output back to HuMANet table"):
        humanet_confidence = _join_input_rows_with_refmet(query_table, refmet_mapped)

    with timed_step("Normalize identifier columns"):
        for col in ["HMDB_ID", "PubChem_CID", "KEGG_ID", "ChEBI_ID", "InChIKey", "SMILES", "Synonyms", "Standardized_Name", "Molecular_Formula", "Exact_Mass", "Super_Class", "Main_Class", "Sub_Class"]:
            if col not in humanet_confidence.columns:
                humanet_confidence[col] = pd.NA
            else:
                humanet_confidence[col] = _naify_stage1_series(humanet_confidence[col])
        humanet_confidence["HMDB_ID"] = humanet_confidence["HMDB_ID"].astype("string").map(_normalize_hmdb_field).astype("string")
        humanet_confidence["PubChem_CID"] = humanet_confidence["PubChem_CID"].astype("string").map(_clean_cid).astype("string")
        humanet_confidence["KEGG_ID"] = _naify_stage1_series(humanet_confidence["KEGG_ID"])
        humanet_confidence["ChEBI_ID"] = _naify_stage1_series(humanet_confidence["ChEBI_ID"])

    with timed_step("Confidence table: offline PubChem annotation"):
        humanet_confidence = annotate_pubchem_offline(humanet_confidence, sqlite_db)

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
                    humanet_confidence[c] = humanet_confidence[c].astype("string").fillna(humanet_confidence[pv].astype("string"))
                    humanet_confidence = humanet_confidence.drop(columns=[pv])

    with timed_step("Confidence table: HMDB fallback fills"):
        humanet_confidence = fill_hmdb_columns(humanet_confidence, hmdb_lite)

    with timed_step("Confidence table: class fallback from PubChem local DB"):
        humanet_confidence = apply_class_fallback_from_pubchem_db(humanet_confidence, sqlite_db)

    with timed_step("Split mapped and unmapped outputs"):
        id_presence = (
            humanet_confidence["PubChem_CID"].astype("string").map(_clean_cid).notna()
            | humanet_confidence["ChEBI_ID"].astype("string").notna()
            | humanet_confidence["HMDB_ID"].astype("string").map(_normalize_hmdb_field).notna()
            | humanet_confidence["KEGG_ID"].astype("string").notna()
        )
        mapped = humanet_confidence.loc[id_presence].copy()
        unmapped = humanet_confidence.loc[~id_presence].copy()

        mapped["HMDB_ID"] = mapped["HMDB_ID"].astype("string").map(_first_hmdb_id).astype("string")
        mapped["matched_name"] = mapped["Query_Name"].astype("string")

        unmapped.loc[:, mapped.columns.difference(["HuMANet_ID", "original_query_name", "Query_Name", "Input_File"], sort=False)] = pd.NA
        unmapped["Database_Source"] = pd.NA
        unmapped["matched_name"] = pd.NA

    id_tracking_table = ensure_cols(
        humanet_confidence.copy(),
        ["HuMANet_ID", "original_query_name", "Query_Name", "matched_name", "Input_File", "PubChem_CID", "HMDB_ID", "KEGG_ID", "ChEBI_ID"],
    )

    mapped = _finalize_public_stage_output(mapped)
    unmapped = _finalize_public_stage_output(unmapped)

    outputs = {
        "mapped": mapped,
        "unmapped": unmapped,
        "id_tracking_table": id_tracking_table,
    }

    with timed_step("Write output CSV files"):
        for name, df in outputs.items():
            df.to_csv(f"{out_prefix}_{name}.csv", index=False)

    excel_sheets = {
        "mapped": mapped,
        "unmapped": unmapped,
    }
    if not skip_stage5_excel:
        excel_sheets.update(build_stage5_annotation_sheets(mapped))

    with timed_step("Write Stage1 Excel workbook"):
        write_excel_workbook(excel_sheets, f"{out_prefix}.xlsx")

    print(f"[TIMER] TOTAL run_stage1: {time.perf_counter() - total_start:.2f}s")
    return outputs


run_stage2 = run_stage1


def main():
    parser = argparse.ArgumentParser(description="HuMANet Stage1 Python pipeline (RefMet-backed Stage1 workflow).")
    parser.add_argument("--stage1_unmapped", default="Refmet_unmapped.csv", help="Stage1 unmapped CSV")
    parser.add_argument("--out_prefix", default="stage1", help="Output file prefix")
    parser.add_argument("--sleep", type=float, default=0.3, help="Sleep seconds between PubChem PUG-View requests")
    parser.add_argument("--pugview_workers", type=int, default=8, help="Thread workers for PubChem PUG-View lookups")
    parser.add_argument("--skip_stage5_excel", action="store_true", help="Skip adding Stage5 sheets to the Excel workbook")
    args = parser.parse_args()

    run_stage1(
        stage1_unmapped_csv=args.stage1_unmapped,
        hmdb_lite_csv=HMDB_LITE_CSV,
        sqlite_db=PUBCHEM_OFFLINE_SQLITE,
        out_prefix=args.out_prefix,
        r_classyfire_bridge=None,
        sleep=args.sleep,
        pugview_workers=args.pugview_workers,
        skip_classyfire=True,
        skip_stage5_excel=args.skip_stage5_excel,
    )


if __name__ == "__main__":
    main()
