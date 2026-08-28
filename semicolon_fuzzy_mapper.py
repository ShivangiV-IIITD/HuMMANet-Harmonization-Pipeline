import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set

import pandas as pd

from extended_annotation_step import (
    _fetch_sqlite_map,
    _fetch_synonyms_sqlite,
    _first_existing_column,
    _lower_name,
    _naify,
    _normalize_chebi,
    _normalize_name,
    _refmet_bridge_query,
    fetch_cids_by_name_sqlite,
    fetch_pubchem_props_api,
)
from pipeline_utils import (
    HAS_RPY2,
    _clean_cid,
    _clean_inchikey,
    _clean_scalar,
    _split_synonyms,
    fetch_pugview_ids_for_cids,
    map_refmet_queries,
)
from resource_config import MW_UNMAPPED_DATABASE_CSV, PUBCHEM_OFFLINE_SQLITE, REFMET_BRIDGE_R

OUTPUT_COLUMNS = [
    "Query_Name",
    "Matched_Query_Part",
    "matched_name",
    "PubChem_CID",
    "Database_Source",
    "Input_File",
    "Input_Database_Source",
    "Study_Folder",
    "HuMANet_ID",
    "Original_HuMANet_ID",
    "HMDB_ID",
    "KEGG_ID",
    "ChEBI_ID",
    "Standardized_Name",
    "IUPAC_Name",
    "Molecular_Formula",
    "Exact_Mass",
    "InChIKey",
    "SMILES",
    "Synonyms",
    "Super_Class",
    "Main_Class",
    "Sub_Class",
]


def _split_query_parts(x: object) -> List[str]:
    s = _clean_scalar(x)
    if not s:
        return []
    parts = [p.strip() for p in s.split(';') if p.strip()]
    out: List[str] = []
    seen: Set[str] = set()
    for part in parts:
        if part not in seen:
            seen.add(part)
            out.append(part)
    return out


def _load_refmet_maps() -> Dict[str, str]:
    mw = pd.read_csv(MW_UNMAPPED_DATABASE_CSV, dtype="string")
    mw.columns = [str(c).strip() for c in mw.columns]
    name_col = _first_existing_column(
        mw.columns.tolist(),
        ["Metabolite.Name", "Metabolite_Name", "Metabolite Name", "Query_Name", "name", "title"],
    )
    refmet_col = _first_existing_column(
        mw.columns.tolist(),
        [
            "RefMet.Name.Standardized.Name.",
            "RefMet_Name_Standardized_Name",
            "RefMet.Standardized.Name",
            "RefMet.Name.Standardized.Name",
            "RefMet_Standardized_Name",
        ],
    )
    if not name_col or not refmet_col:
        return {}

    work = mw[[name_col, refmet_col]].copy()
    work["__query_key"] = work[name_col].map(_normalize_name)
    work[refmet_col] = work[refmet_col].astype("string").str.replace("*", "", regex=False)
    work[refmet_col] = work[refmet_col].astype("string").str.replace("&", "", regex=False)
    work[refmet_col] = work[refmet_col].map(_lower_name)
    work = work.dropna(subset=["__query_key", refmet_col])
    work = work.drop_duplicates(subset=["__query_key"], keep="first")
    return dict(zip(work["__query_key"], work[refmet_col].astype("string")))


def _refmet_bridge_lookup(parts: Sequence[str]) -> Dict[str, str]:
    if not parts or not HAS_RPY2:
        return {}
    queries = [_refmet_bridge_query(p) for p in parts]
    queries = [q for q in queries if _clean_scalar(q)]
    if not queries:
        return {}
    bridge_df = map_refmet_queries(queries, REFMET_BRIDGE_R)
    if bridge_df.empty:
        return {}

    cols = [str(c) for c in bridge_df.columns]
    input_col = _first_existing_column(cols, ["Input.name", "Input_name", "input_name", "query_name"])
    std_col = _first_existing_column(cols, ["Standardized.name", "Standardized_name", "standardized_name"])
    if not input_col or not std_col:
        return {}

    work = bridge_df[[input_col, std_col]].copy()
    work["__query_key"] = work[input_col].map(_refmet_bridge_query)
    work[std_col] = work[std_col].astype("string").str.replace("*", "", regex=False)
    work[std_col] = work[std_col].astype("string").str.replace("&", "", regex=False)
    work[std_col] = work[std_col].map(_lower_name)
    work = work.dropna(subset=["__query_key", std_col])
    work = work.drop_duplicates(subset=["__query_key"], keep="first")
    return dict(zip(work["__query_key"], work[std_col].astype("string")))


def _resolve_row(query_name: object, refmet_map: Dict[str, str]) -> Dict[str, object]:
    parts = _split_query_parts(query_name)
    bridge_map: Dict[str, str] = {}
    out = {
        "Matched_Query_Part": pd.NA,
        "matched_name": pd.NA,
        "PubChem_CID": pd.NA,
    }
    if not parts:
        return out

    for part in parts:
        qkey = _normalize_name(part)
        refmet_name = refmet_map.get(qkey) if qkey else None
        if not refmet_name:
            if not bridge_map:
                bridge_map = _refmet_bridge_lookup(parts)
            refmet_name = bridge_map.get(_refmet_bridge_query(part))

        search_names = []
        if _clean_scalar(refmet_name):
            search_names.append(refmet_name)
        if _clean_scalar(part):
            search_names.append(part)

        seen_names: Set[str] = set()
        ordered_names: List[str] = []
        for name in search_names:
            key = _lower_name(name)
            if key and key not in seen_names:
                seen_names.add(key)
                ordered_names.append(name)

        for lookup_name in ordered_names:
            cids = fetch_cids_by_name_sqlite(PUBCHEM_OFFLINE_SQLITE, lookup_name, limit=20)
            cid = next((c for c in cids if _clean_cid(c)), None)
            if cid:
                out["Matched_Query_Part"] = part
                out["matched_name"] = lookup_name
                out["PubChem_CID"] = cid
                return out

    return out


def _ensure_output_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in OUTPUT_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    return out


def run_mapper(input_csv: str, output_csv: str, query_col: str = "Query_Name", pubchem_workers: int = 16) -> pd.DataFrame:
    df = pd.read_csv(input_csv, dtype="string")
    if query_col not in df.columns:
        raise ValueError(f'Missing query column: {query_col}')

    refmet_map = _load_refmet_maps()

    results = []
    for _, row in df.iterrows():
        base = {col: row[col] if col in row.index else pd.NA for col in df.columns}
        base.setdefault("Query_Name", row[query_col])
        resolved = _resolve_row(row[query_col], refmet_map)
        base.update(resolved)
        results.append(base)

    out = pd.DataFrame(results).astype("string")
    out = _ensure_output_columns(out)
    out["PubChem_CID"] = out["PubChem_CID"].map(_clean_cid).astype("string")

    cids = {c for c in out["PubChem_CID"].dropna().tolist() if _clean_cid(c)}
    if cids:
        title_map = _fetch_sqlite_map(PUBCHEM_OFFLINE_SQLITE, "cid_title", "cid", "title", cids)
        iupac_map = _fetch_sqlite_map(PUBCHEM_OFFLINE_SQLITE, "cid_iupac", "cid", "iupac_name", cids)
        smiles_map = _fetch_sqlite_map(PUBCHEM_OFFLINE_SQLITE, "cid_smiles", "cid", "smiles", cids)
        inchi_map = _fetch_sqlite_map(PUBCHEM_OFFLINE_SQLITE, "cid_inchi", "cid", "inchikey", cids)
        synonyms_map = _fetch_synonyms_sqlite(PUBCHEM_OFFLINE_SQLITE, cids)

        out["Standardized_Name"] = _naify(out["Standardized_Name"]).combine_first(out["PubChem_CID"].map(title_map))
        out["IUPAC_Name"] = _naify(out["IUPAC_Name"]).combine_first(out["PubChem_CID"].map(iupac_map))
        out["SMILES"] = _naify(out["SMILES"]).combine_first(out["PubChem_CID"].map(smiles_map))
        out["InChIKey"] = _naify(out["InChIKey"]).combine_first(out["PubChem_CID"].map(inchi_map))
        out["Synonyms"] = _naify(out["Synonyms"]).combine_first(out["PubChem_CID"].map(synonyms_map))

        api = fetch_pubchem_props_api(sorted(cids), workers=pubchem_workers)
        if not api.empty:
            out = out.merge(api, on="PubChem_CID", how="left", suffixes=("", "_api"))
            for col in ["Standardized_Name", "Molecular_Formula", "Exact_Mass", "InChIKey", "SMILES"]:
                api_col = f"{col}_api"
                if api_col in out.columns:
                    out[col] = _naify(out[col]).combine_first(_naify(out[api_col]))
                    out = out.drop(columns=[api_col])

        ann = fetch_pugview_ids_for_cids(sorted(cids), sleep=0.2, max_workers=pubchem_workers)
        if not ann.empty:
            out = out.merge(ann, on="PubChem_CID", how="left", suffixes=("", "_pv"))
            for col in ["HMDB_ID", "KEGG_ID", "ChEBI_ID"]:
                pv_col = f"{col}_pv"
                if pv_col in out.columns:
                    out[col] = _naify(out[col]).combine_first(_naify(out[pv_col]))
                    out = out.drop(columns=[pv_col])

    out["InChIKey"] = out["InChIKey"].map(_clean_inchikey).astype("string")
    out["ChEBI_ID"] = out["ChEBI_ID"].map(_normalize_chebi).astype("string")
    out["Database_Source"] = _naify(out["Database_Source"])
    out.loc[out["Database_Source"].isna() & out["PubChem_CID"].notna(), "Database_Source"] = "PubChem"
    out["Super_Class"] = pd.NA
    out["Main_Class"] = pd.NA
    out["Sub_Class"] = pd.NA

    out = out[OUTPUT_COLUMNS]
    out.to_csv(output_csv, index=False)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Semicolon-aware lightweight Stage3-style mapper")
    parser.add_argument("--input_csv", required=True, help="Input CSV with Query_Name column")
    parser.add_argument("--output_csv", required=True, help="Output CSV path")
    parser.add_argument("--query_col", default="Query_Name", help="Query-name column")
    parser.add_argument("--pubchem_workers", type=int, default=16, help="Workers for PubChem enrichment")
    args = parser.parse_args()

    out = run_mapper(
        input_csv=args.input_csv,
        output_csv=args.output_csv,
        query_col=args.query_col,
        pubchem_workers=args.pubchem_workers,
    )
    print(f"Rows written: {len(out)}")
    print(f"Saved to: {args.output_csv}")


if __name__ == "__main__":
    main()
