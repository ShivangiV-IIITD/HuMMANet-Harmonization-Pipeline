import argparse
import time
from typing import Dict, Optional, Sequence, Tuple

import pandas as pd

from knowledge_annotation_step import build_stage5_annotation_sheets
from pipeline_utils import _clean_scalar, _split_synonyms, prepare_stage1_input_dataframe, timed_step, write_excel_workbook
from resource_config import HUMANNET_LIBRARY_CSV


PUBLIC_STAGE_COLUMNS = [
    "HuMANet_ID", "original_query_name", "Query_Name", "Database_Source", "InChIKey", "HMDB_ID", "PubChem_CID", "KEGG_ID", "ChEBI_ID",
    "Molecular_Formula", "Exact_Mass", "Super_Class", "Main_Class", "Sub_Class", "Standardized_Name", "SMILES",
    "Synonyms", "matched_name", "Study_Folder",
]

LIBRARY_TRANSFER_COLUMNS = [
    "Database_Source", "InChIKey", "HMDB_ID", "PubChem_CID", "KEGG_ID", "ChEBI_ID",
    "Molecular_Formula", "Exact_Mass", "Super_Class", "Main_Class", "Sub_Class",
    "Standardized_Name", "SMILES", "Synonyms",
]

SOURCE_PRIORITY = {
    "all_studies_mapped": 0,
    "HuMANet_mimedb_all": 1,
    "HuMANet_monoculture_all": 2,
    "all_studies_nonconfidence": 3,
}


def _norm_query_name(x: object) -> object:
    s = _clean_scalar(x)
    if not s:
        return pd.NA
    return " ".join(str(s).strip().lower().split())


def _ensure_public_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
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
    return out[list(columns)].copy()


def _split_library_entries(x: object) -> Sequence[str]:
    s = _clean_scalar(x)
    if not s:
        return []
    vals = _split_synonyms(s)
    return vals if vals else [str(s).strip()]


def _prepare_library_rows(library_df: pd.DataFrame) -> pd.DataFrame:
    work = library_df.copy()
    if "Query_Name" not in work.columns:
        raise ValueError('HuMANet library must contain a "Query_Name" column.')

    work = work.loc[:, [c for c in work.columns if not str(c).startswith("Unnamed:")]].copy()
    work["__query_key"] = work["Query_Name"].map(_norm_query_name)
    work = work.dropna(subset=["__query_key"]).copy()

    score_cols = [
        "PubChem_CID", "HMDB_ID", "KEGG_ID", "ChEBI_ID", "InChIKey", "SMILES",
        "Standardized_Name", "Molecular_Formula", "Exact_Mass", "Synonyms",
    ]
    for col in score_cols:
        if col not in work.columns:
            work[col] = pd.NA
    work["__score"] = work[score_cols].notna().sum(axis=1)

    if "source_priority_used" not in work.columns:
        work["source_priority_used"] = pd.NA

    work["__source_priority_rank"] = work["source_priority_used"].map(
        lambda x: SOURCE_PRIORITY.get(_clean_scalar(x), 99)
    )
    work = work.sort_values(
        by=["__query_key", "__source_priority_rank", "__score", "source_priority_used"],
        ascending=[True, True, False, True],
        na_position="last",
        kind="stable",
    )
    return work.copy()


def _build_lookup_map(library_df: pd.DataFrame, column: str) -> Dict[str, Dict[str, object]]:
    if column not in library_df.columns:
        return {}

    out: Dict[str, Dict[str, object]] = {}
    for _, row in library_df.iterrows():
        for token in _split_library_entries(row[column]):
            key = _norm_query_name(token)
            key = _clean_scalar(key)
            if key and key not in out:
                out[key] = row.to_dict()
    return out


def _lookup_library_row(
    query_key: object,
    query_map: Dict[str, Dict[str, object]],
    standardized_map: Dict[str, Dict[str, object]],
    synonym_map: Dict[str, Dict[str, object]],
) -> Tuple[Optional[str], Optional[Dict[str, object]]]:
    key = _clean_scalar(query_key)
    if not key:
        return None, None
    if key in query_map:
        return "Query_Name", query_map[key]
    if key in standardized_map:
        return "Standardized_Name", standardized_map[key]
    if key in synonym_map:
        return "Synonyms", synonym_map[key]
    return None, None


def run_library_lookup(input_csv: str, library_csv: str, out_prefix: str) -> Dict[str, pd.DataFrame]:
    total_start = time.perf_counter()

    with timed_step("Load input + HuMANet library"):
        base = pd.read_csv(input_csv, dtype="string")
        base = prepare_stage1_input_dataframe(base, input_csv)
        if "original_query_name" not in base.columns:
            base["original_query_name"] = base["Query_Name"].astype("string")
        library_df = pd.read_csv(library_csv, dtype="string")
        library_df = _prepare_library_rows(library_df)

    with timed_step("Library lookup by query, standardized name, then synonyms"):
        base["__query_key"] = base["Query_Name"].map(_norm_query_name)
        query_map = _build_lookup_map(library_df, "Query_Name")
        standardized_map = _build_lookup_map(library_df, "Standardized_Name")
        synonym_map = _build_lookup_map(library_df, "Synonyms")

        hit_sources = []
        hit_rows = []
        for query_key in base["__query_key"].tolist():
            hit_source, hit_row = _lookup_library_row(query_key, query_map, standardized_map, synonym_map)
            hit_sources.append(hit_source)
            hit_rows.append(hit_row or {})

        merged = base.copy()
        merged["__lookup_source"] = pd.Series(hit_sources, index=merged.index, dtype="string")
        hit_df = pd.DataFrame(hit_rows, index=merged.index)

        for col in LIBRARY_TRANSFER_COLUMNS:
            if col in hit_df.columns:
                if col not in merged.columns:
                    merged[col] = pd.NA
                merged[col] = merged[col].astype("string").fillna(hit_df[col].astype("string"))

        merged["matched_name"] = pd.Series([pd.NA] * len(merged), dtype="string")
        library_hit = merged["__lookup_source"].notna()
        mapped = merged.loc[library_hit].copy()
        unmapped = merged.loc[~library_hit].copy()

        if not mapped.empty:
            mapped["matched_name"] = mapped["Query_Name"].astype("string")
        if not unmapped.empty:
            unmapped["matched_name"] = pd.NA

        mapped = mapped.drop(columns=["__query_key", "__lookup_source"], errors="ignore")
        unmapped = unmapped.drop(columns=["__query_key", "__lookup_source"], errors="ignore")

    with timed_step("Finalize outputs"):
        mapped_public = _ensure_public_columns(mapped, PUBLIC_STAGE_COLUMNS)
        unmapped_public = _ensure_public_columns(unmapped, PUBLIC_STAGE_COLUMNS)

    with timed_step("Write output CSV files"):
        mapped_public.to_csv(f"{out_prefix}_mapped.csv", index=False)
        unmapped_public.to_csv(f"{out_prefix}_unmapped.csv", index=False)

    with timed_step("Write Library Lookup Excel workbook"):
        workbook = {"mapped": mapped_public, "unmapped": unmapped_public}
        if not mapped_public.empty:
            workbook.update(build_stage5_annotation_sheets(mapped_public))
        else:
            workbook.update(
                {
                    "species": pd.DataFrame(),
                    "disease": pd.DataFrame(),
                    "pathways": pd.DataFrame(),
                }
            )
        write_excel_workbook(workbook, f"{out_prefix}.xlsx")

    print(f"[TIMER] TOTAL run_library_lookup: {time.perf_counter() - total_start:.2f}s")
    return {"mapped": mapped_public, "unmapped": unmapped_public}


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the initial HuMANet library query-name lookup stage")
    parser.add_argument("--input_csv", required=True, help="Raw pipeline input CSV")
    parser.add_argument("--library_csv", default=HUMANNET_LIBRARY_CSV, help="HuMANet library CSV")
    parser.add_argument("--out_prefix", required=True, help="Output file prefix")
    args = parser.parse_args()

    run_library_lookup(
        input_csv=args.input_csv,
        library_csv=args.library_csv,
        out_prefix=args.out_prefix,
    )


if __name__ == "__main__":
    main()
