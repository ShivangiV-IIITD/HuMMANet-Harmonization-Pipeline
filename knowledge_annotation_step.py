import argparse
import os
from pathlib import Path
from functools import lru_cache
from glob import glob
from typing import Iterable, List, Optional, Sequence

import pandas as pd

from drug_similarity_step import build_drug_annotation_sheets
from resource_config import HMDB_LITE_CSV, SMPDB_DIR, SMPDB_MAPPING_CACHE, SPECIES_INFO_CSV


NULL_STRINGS = {"", "nan", "na", "n/a", "null", "<na>", "-", "none"}
DISEASE_SHEET_COLUMNS = [
    "HuMANet_ID",
    "Query_Name",
    "Study_Folder",
    "Annotation_ID_Type",
    "Annotation_ID_Value",
    "InChIKey",
    "PubChem_CID",
    "HMDB_ID",
    "HMDB_Biofluid",
    "HMDB_Disease",
]
SPECIES_SHEET_COLUMNS = [
    "HuMANet_ID",
    "Query_Name",
    "Study_Folder",
    "Annotation_ID_Type",
    "Annotation_ID_Value",
    "InChIKey",
    "PubChem_CID",
    "HMDB_ID",
    "Species",
    "Species_Biofluids",
]
PATHWAY_SHEET_COLUMNS = [
    "HuMANet_ID",
    "Query_Name",
    "Study_Folder",
    "Annotation_ID_Type",
    "Annotation_ID_Value",
    "InChIKey",
    "PubChem_CID",
    "HMDB_ID",
    "SMPDB_ID",
    "Pathway_Name",
    "Pathway_Subject",
    "SMPDB_Metabolite_ID",
]


def clean_scalar(value: object) -> Optional[str]:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text.lower() in NULL_STRINGS:
        return None
    return text


def clean_series(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .replace(
            {
                "": pd.NA,
                "nan": pd.NA,
                "NaN": pd.NA,
                "NA": pd.NA,
                "N/A": pd.NA,
                "null": pd.NA,
                "<NA>": pd.NA,
                "-": pd.NA,
                "None": pd.NA,
            }
        )
    )


def unique_join(values: Iterable[object], sep: str = "|") -> str:
    seen = set()
    ordered: List[str] = []
    for value in values:
        text = clean_scalar(value)
        if text and text not in seen:
            seen.add(text)
            ordered.append(text)
    return sep.join(ordered)


def split_multi_value(value: object) -> List[str]:
    text = clean_scalar(value)
    if not text:
        return []
    parts = [part.strip() for part in text.split(";")]
    return [part for part in parts if clean_scalar(part)]


def collect_values(df: pd.DataFrame, column: str) -> List[str]:
    if column not in df.columns or df.empty:
        return []

    seen = set()
    ordered: List[str] = []
    for value in df[column]:
        for item in split_multi_value(value):
            if item not in seen:
                seen.add(item)
                ordered.append(item)
    return ordered


def ensure_columns(df: pd.DataFrame, columns: Sequence[str], name: str) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")


def build_lookup(meta: pd.DataFrame, column: str) -> pd.DataFrame:
    usable = meta[meta[column].notna()].copy()
    return usable.set_index(column, drop=False)


def fetch_matches(lookup: pd.DataFrame, key: Optional[str]) -> Optional[pd.DataFrame]:
    if not key or key not in lookup.index:
        return None

    matched = lookup.loc[key]
    if isinstance(matched, pd.Series):
        matched = matched.to_frame().T
    return matched.reset_index(drop=True)


def aggregate_hmdb_annotations(hmdb_df: pd.DataFrame, key_column: str, value_column: str) -> pd.DataFrame:
    usable = hmdb_df[[key_column, value_column]].copy()
    usable[key_column] = clean_series(usable[key_column])
    usable = usable.dropna(subset=[key_column])
    if usable.empty:
        return pd.DataFrame(columns=[key_column, value_column])

    records = []
    for key, group in usable.groupby(key_column, sort=False):
        values = collect_values(group, value_column)
        if values:
            records.append({key_column: key, value_column: "|".join(values)})

    return pd.DataFrame(records)


@lru_cache(maxsize=4)
def _load_species_lookup(species_file: str) -> pd.DataFrame:
    species_df = pd.read_csv(
        species_file,
        usecols=["InChIKey", "Species", "biofluids"],
        low_memory=False,
    )
    return aggregate_species_annotations(species_df).set_index("InChIKey", drop=False)


@lru_cache(maxsize=4)
def _load_hmdb_raw_for_sheets(hmdb_path: str) -> pd.DataFrame:
    return pd.read_csv(
        hmdb_path,
        usecols=[
            "biofluid",
            "hmdb_metabolite_accession",
            "hmdb_metabolite_inchikey",
            "hmdb_metabolite_pubchem_compound_id",
            "hmdb_metabolite_diseases_disease_name",
        ],
        low_memory=False,
    )


@lru_cache(maxsize=4)
def _load_species_raw_for_sheets(species_file: str) -> pd.DataFrame:
    return pd.read_csv(
        species_file,
        usecols=["InChIKey", "Species", "biofluids"],
        low_memory=False,
    )


def _build_smpdb_mapping_cache(smpdb_dir: str) -> pd.DataFrame:
    patterns = [
        os.path.join(smpdb_dir, "*.csv"),
        os.path.join(smpdb_dir, "*", "*.csv"),
        os.path.join(smpdb_dir, "**", "*.csv"),
    ]
    files: List[str] = []
    for pattern in patterns:
        files.extend(glob(pattern, recursive=True))
    files = sorted(set(files))

    print(f"Total SMPDB files found: {len(files)}")
    if not files:
        raise FileNotFoundError(f"No SMPDB metabolite files found under: {smpdb_dir}")

    keep_cols = [
        "SMPDB ID",
        "Pathway Name",
        "Pathway Subject",
        "Metabolite ID",
        "HMDB ID",
        "KEGG ID",
        "ChEBI ID",
    ]
    frames: List[pd.DataFrame] = []
    for i, file_path in enumerate(files, 1):
        if i % 1000 == 0:
            print(f"Processed {i}/{len(files)} files")
        df = pd.read_csv(file_path, usecols=keep_cols, low_memory=False)
        ensure_columns(df, keep_cols, os.path.basename(file_path))
        frames.append(df)

    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=keep_cols)
    for col in ["HMDB ID", "KEGG ID", "ChEBI ID"]:
        out[col] = clean_series(out[col])
    return out


@lru_cache(maxsize=2)
def _load_smpdb_mapping_table(smpdb_dir: str, cache_path: str = SMPDB_MAPPING_CACHE) -> pd.DataFrame:
    cache_file = Path(cache_path)
    if cache_file.exists():
        print(f"Loading cached SMPDB mapping: {cache_file}")
        return pd.read_pickle(cache_file)

    print("Building cached SMPDB mapping table...")
    mapping = _build_smpdb_mapping_cache(smpdb_dir)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_pickle(cache_file)
    print(f"Saved SMPDB mapping cache: {cache_file}")
    return mapping


def aggregate_species_annotations(species_df: pd.DataFrame) -> pd.DataFrame:
    usable = species_df[["InChIKey", "Species", "biofluids"]].copy()
    usable["InChIKey"] = clean_series(usable["InChIKey"])
    usable = usable.dropna(subset=["InChIKey"])
    if usable.empty:
        return pd.DataFrame(columns=["InChIKey", "Species", "biofluids"])

    records = []
    for key, group in usable.groupby("InChIKey", sort=False):
        records.append(
            {
                "InChIKey": key,
                "Species": "|".join(collect_values(group, "Species")),
                "biofluids": "|".join(collect_values(group, "biofluids")),
            }
        )

    return pd.DataFrame(records)


@lru_cache(maxsize=4)
def load_hmdb_annotation_tables(hmdb_path: str) -> dict:
    hmdb = pd.read_csv(
        hmdb_path,
        usecols=[
            "biofluid",
            "hmdb_metabolite_accession",
            "hmdb_metabolite_inchikey",
            "hmdb_metabolite_pubchem_compound_id",
            "hmdb_metabolite_diseases_disease_name",
        ],
        low_memory=False,
    )

    by_accession = aggregate_hmdb_annotations(
        hmdb,
        "hmdb_metabolite_accession",
        "biofluid",
    ).merge(
        aggregate_hmdb_annotations(
            hmdb,
            "hmdb_metabolite_accession",
            "hmdb_metabolite_diseases_disease_name",
        ),
        on="hmdb_metabolite_accession",
        how="outer",
    )

    by_inchikey = aggregate_hmdb_annotations(
        hmdb,
        "hmdb_metabolite_inchikey",
        "biofluid",
    ).merge(
        aggregate_hmdb_annotations(
            hmdb,
            "hmdb_metabolite_inchikey",
            "hmdb_metabolite_diseases_disease_name",
        ),
        on="hmdb_metabolite_inchikey",
        how="outer",
    )

    by_pubchem = aggregate_hmdb_annotations(
        hmdb,
        "hmdb_metabolite_pubchem_compound_id",
        "biofluid",
    ).merge(
        aggregate_hmdb_annotations(
            hmdb,
            "hmdb_metabolite_pubchem_compound_id",
            "hmdb_metabolite_diseases_disease_name",
        ),
        on="hmdb_metabolite_pubchem_compound_id",
        how="outer",
    )

    return {
        "by_accession": by_accession.set_index("hmdb_metabolite_accession", drop=False),
        "by_inchikey": by_inchikey.set_index("hmdb_metabolite_inchikey", drop=False),
        "by_pubchem": by_pubchem.set_index("hmdb_metabolite_pubchem_compound_id", drop=False),
    }


def collect_hmdb_info(
    hmdb_tables: dict,
    hmdb_ids: Sequence[Optional[str]],
    inchikeys: Sequence[Optional[str]],
    pubchem_ids: Sequence[Optional[str]],
) -> tuple[str, str]:
    biofluids: List[str] = []
    diseases: List[str] = []

    def add_from_table(table: pd.DataFrame, keys: Sequence[Optional[str]]) -> None:
        nonlocal biofluids, diseases
        for key in keys:
            cleaned = clean_scalar(key)
            if not cleaned or cleaned not in table.index:
                continue

            matched = table.loc[cleaned]
            if isinstance(matched, pd.Series):
                matched = matched.to_frame().T

            biofluids.extend(collect_values(matched, "biofluid"))
            diseases.extend(collect_values(matched, "hmdb_metabolite_diseases_disease_name"))

    add_from_table(hmdb_tables["by_accession"], hmdb_ids)
    add_from_table(hmdb_tables["by_inchikey"], inchikeys)
    add_from_table(hmdb_tables["by_pubchem"], pubchem_ids)

    return unique_join(biofluids), unique_join(diseases)


def collect_species_info(species_lookup: pd.DataFrame, inchikeys: Sequence[Optional[str]]) -> tuple[str, str]:
    species_values: List[str] = []
    biofluid_values: List[str] = []

    for key in inchikeys:
        cleaned = clean_scalar(key)
        if not cleaned or cleaned not in species_lookup.index:
            continue

        matched = species_lookup.loc[cleaned]
        if isinstance(matched, pd.Series):
            matched = matched.to_frame().T

        species_values.extend(collect_values(matched, "Species"))
        biofluid_values.extend(collect_values(matched, "biofluids"))

    return unique_join(species_values), unique_join(biofluid_values)


def _split_pipe_value(value: object) -> List[str]:
    text = clean_scalar(value)
    if not text:
        return []
    return [part.strip() for part in str(text).split("|") if clean_scalar(part)]


def _unique_join_with_sep(values: Iterable[object], sep: str = ";") -> Optional[str]:
    seen = set()
    ordered: List[str] = []
    for value in values:
        text = clean_scalar(value)
        if text and text not in seen:
            seen.add(text)
            ordered.append(text)
    return sep.join(ordered) if ordered else None


def _split_sep_values(values: Iterable[object], sep: str = ";") -> List[str]:
    out: List[str] = []
    for value in values:
        text = clean_scalar(value)
        if not text:
            continue
        out.extend([part.strip() for part in str(text).split(sep) if clean_scalar(part)])
    return out


def _build_species_sheet(stage5_df: pd.DataFrame, species_df: pd.DataFrame) -> pd.DataFrame:
    species_df = species_df.copy()
    species_df["InChIKey"] = clean_series(species_df["InChIKey"])
    species_df["Species"] = clean_series(species_df["Species"])
    species_df["biofluids"] = clean_series(species_df["biofluids"])

    rows = []
    for _, row in stage5_df.iterrows():
        inchikeys = _split_pipe_value(row.get("InChIKey"))
        if not inchikeys:
            continue

        matched = species_df[species_df["InChIKey"].isin(inchikeys)].copy()
        matched = matched.dropna(subset=["Species"])
        if matched.empty:
            continue

        for species_name, grp in matched.groupby("Species", sort=False):
            rows.append(
                {
                    "HuMANet_ID": row["HuMANet_ID"],
                    "Query_Name": row["Query_Name"],
                    "Study_Folder": row["Study_Folder"],
                    "Annotation_ID_Type": row["Annotation_ID_Type"],
                    "Annotation_ID_Value": row["Annotation_ID_Value"],
                    "InChIKey": row["InChIKey"],
                    "PubChem_CID": row["PubChem_CID"],
                    "HMDB_ID": row["HMDB_ID"],
                    "Species": species_name,
                    "Species_Biofluids": _unique_join_with_sep(grp["biofluids"].tolist(), sep=";"),
                }
            )

    if not rows:
        return pd.DataFrame(columns=SPECIES_SHEET_COLUMNS)

    species = pd.DataFrame(rows)
    group_cols = [
        "HuMANet_ID",
        "Query_Name",
        "Study_Folder",
        "InChIKey",
        "PubChem_CID",
        "HMDB_ID",
        "Species",
    ]
    species = (
        species.groupby(group_cols, dropna=False, as_index=False)
        .agg(
            {
                "Annotation_ID_Type": lambda s: _unique_join_with_sep(s.tolist(), sep=";"),
                "Annotation_ID_Value": lambda s: _unique_join_with_sep(s.tolist(), sep=";"),
                "Species_Biofluids": lambda s: _unique_join_with_sep(_split_sep_values(s.tolist(), sep=";"), sep=";"),
            }
        )
    )
    return species[SPECIES_SHEET_COLUMNS].drop_duplicates().reset_index(drop=True)


def _build_disease_sheet(stage5_df: pd.DataFrame, hmdb_df: pd.DataFrame) -> pd.DataFrame:
    hmdb_df = hmdb_df.copy()
    for col in [
        "biofluid",
        "hmdb_metabolite_accession",
        "hmdb_metabolite_inchikey",
        "hmdb_metabolite_pubchem_compound_id",
        "hmdb_metabolite_diseases_disease_name",
    ]:
        hmdb_df[col] = clean_series(hmdb_df[col])

    rows = []
    for _, row in stage5_df.iterrows():
        hmdb_ids = _split_pipe_value(row.get("HMDB_ID"))
        inchikeys = _split_pipe_value(row.get("InChIKey"))
        pubchem_ids = _split_pipe_value(row.get("PubChem_CID"))

        matched_frames = []
        if hmdb_ids:
            matched_frames.append(hmdb_df[hmdb_df["hmdb_metabolite_accession"].isin(hmdb_ids)])
        if inchikeys:
            matched_frames.append(hmdb_df[hmdb_df["hmdb_metabolite_inchikey"].isin(inchikeys)])
        if pubchem_ids:
            matched_frames.append(hmdb_df[hmdb_df["hmdb_metabolite_pubchem_compound_id"].isin(pubchem_ids)])
        if not matched_frames:
            continue

        matched = pd.concat(matched_frames, ignore_index=True).drop_duplicates()
        matched = matched.dropna(subset=["hmdb_metabolite_diseases_disease_name"])
        if matched.empty:
            continue

        for disease_name, grp in matched.groupby("hmdb_metabolite_diseases_disease_name", sort=False):
            rows.append(
                {
                    "HuMANet_ID": row["HuMANet_ID"],
                    "Query_Name": row["Query_Name"],
                    "Study_Folder": row["Study_Folder"],
                    "Annotation_ID_Type": row["Annotation_ID_Type"],
                    "Annotation_ID_Value": row["Annotation_ID_Value"],
                    "InChIKey": row["InChIKey"],
                    "PubChem_CID": row["PubChem_CID"],
                    "HMDB_ID": row["HMDB_ID"],
                    "HMDB_Biofluid": _unique_join_with_sep(grp["biofluid"].tolist(), sep=";"),
                    "HMDB_Disease": disease_name,
                }
            )

    if not rows:
        return pd.DataFrame(columns=DISEASE_SHEET_COLUMNS)

    disease = pd.DataFrame(rows)
    group_cols = [
        "HuMANet_ID",
        "Query_Name",
        "Study_Folder",
        "InChIKey",
        "PubChem_CID",
        "HMDB_ID",
        "HMDB_Disease",
    ]
    disease = (
        disease.groupby(group_cols, dropna=False, as_index=False)
        .agg(
            {
                "Annotation_ID_Type": lambda s: _unique_join_with_sep(s.tolist(), sep=";"),
                "Annotation_ID_Value": lambda s: _unique_join_with_sep(s.tolist(), sep=";"),
                "HMDB_Biofluid": lambda s: _unique_join_with_sep(_split_sep_values(s.tolist(), sep=";"), sep=";"),
            }
        )
    )
    return disease[DISEASE_SHEET_COLUMNS].drop_duplicates().reset_index(drop=True)


def default_output_path(input_path: str) -> str:
    base_dir = os.path.dirname(input_path)
    stem = os.path.splitext(os.path.basename(input_path))[0]
    return os.path.join(base_dir, f"{stem}_stage5_smpdb_mapping.csv")


def prepare_input_meta(meta: pd.DataFrame) -> pd.DataFrame:
    meta = meta.copy()
    if "Study_Folder" not in meta.columns and "Input_File" in meta.columns:
        meta["Study_Folder"] = meta["Input_File"]
    ensure_columns(
        meta,
        ["HuMANet_ID", "Query_Name", "Study_Folder", "InChIKey", "PubChem_CID", "HMDB_ID", "KEGG_ID", "ChEBI_ID"],
        "Input metadata",
    )

    for column in ["HuMANet_ID", "Query_Name", "Study_Folder", "InChIKey", "PubChem_CID", "HMDB_ID", "KEGG_ID", "ChEBI_ID"]:
        meta[column] = clean_series(meta[column])

    if not meta["HuMANet_ID"].dropna().str.startswith("HMAN").all():
        raise AssertionError("HuMANet_ID column contains non-HMAN values")
    return meta


def build_stage5_mapping_from_df(
    meta: pd.DataFrame,
    smpdb_dir: str,
    hmdb_lite_file: str,
    species_file: str,
) -> pd.DataFrame:
    meta = prepare_input_meta(meta)

    hmdb_lookup = build_lookup(meta, "HMDB_ID")
    kegg_lookup = build_lookup(meta, "KEGG_ID")
    chebi_lookup = build_lookup(meta, "ChEBI_ID")

    print("Loaded HuMANet mappings:")
    print("HMDB:", hmdb_lookup.shape[0])
    print("KEGG:", kegg_lookup.shape[0])
    print("ChEBI:", chebi_lookup.shape[0])

    print("\nLoading HMDB Lite annotations...")
    hmdb_tables = load_hmdb_annotation_tables(hmdb_lite_file)

    print("Loading species annotations...")
    species_lookup = _load_species_lookup(species_file)

    smpdb_df = _load_smpdb_mapping_table(smpdb_dir)

    hmdb_keys = set(hmdb_lookup.index.astype(str)) if not hmdb_lookup.empty else set()
    kegg_keys = set(kegg_lookup.index.astype(str)) if not kegg_lookup.empty else set()
    chebi_keys = set(chebi_lookup.index.astype(str)) if not chebi_lookup.empty else set()

    matched_frames = []
    if hmdb_keys:
        matched_frames.append(smpdb_df[smpdb_df["HMDB ID"].isin(hmdb_keys)])
    if kegg_keys:
        matched_frames.append(smpdb_df[smpdb_df["KEGG ID"].isin(kegg_keys)])
    if chebi_keys:
        matched_frames.append(smpdb_df[smpdb_df["ChEBI ID"].isin(chebi_keys)])

    if matched_frames:
        smpdb_df = pd.concat(matched_frames, ignore_index=True).drop_duplicates()
    else:
        smpdb_df = smpdb_df.iloc[0:0].copy()

    print(f"Filtered SMPDB rows to evaluate: {len(smpdb_df)}")
    results = []

    for _, row in smpdb_df.iterrows():
            matched_types: List[str] = []
            matched_values: List[str] = []
            matched_q_frames: List[pd.DataFrame] = []

            hmdb = clean_scalar(row["HMDB ID"])
            hmdb_match = fetch_matches(hmdb_lookup, hmdb)
            if hmdb_match is not None:
                matched_types.append("HMDB")
                matched_values.append(hmdb)
                matched_q_frames.append(hmdb_match)

            kegg = clean_scalar(row["KEGG ID"])
            kegg_match = fetch_matches(kegg_lookup, kegg)
            if kegg_match is not None:
                matched_types.append("KEGG")
                matched_values.append(kegg)
                matched_q_frames.append(kegg_match)

            chebi = clean_scalar(row["ChEBI ID"])
            chebi_match = fetch_matches(chebi_lookup, chebi)
            if chebi_match is not None:
                matched_types.append("ChEBI")
                matched_values.append(chebi)
                matched_q_frames.append(chebi_match)

            if not matched_q_frames:
                continue

            all_q = pd.concat(matched_q_frames, ignore_index=True)

            humanet_ids = sorted(all_q["HuMANet_ID"].dropna().astype(str).unique())
            query_names = sorted(all_q["Query_Name"].dropna().astype(str).unique())
            input_files = sorted(all_q["Study_Folder"].dropna().astype(str).unique())
            query_inchikeys = sorted(all_q["InChIKey"].dropna().astype(str).unique())
            query_pubchems = sorted(all_q["PubChem_CID"].dropna().astype(str).unique())
            query_hmdbs = sorted(all_q["HMDB_ID"].dropna().astype(str).unique())

            assert all(h.startswith("HMAN") for h in humanet_ids)
            assert not any("\n" in h for h in humanet_ids)

            hmdb_ids_for_enrichment = query_hmdbs + ([hmdb] if hmdb else [])
            inchikeys_for_enrichment = query_inchikeys
            pubchems_for_enrichment = query_pubchems

            hmdb_biofluid, hmdb_diseases = collect_hmdb_info(
                hmdb_tables=hmdb_tables,
                hmdb_ids=hmdb_ids_for_enrichment,
                inchikeys=inchikeys_for_enrichment,
                pubchem_ids=pubchems_for_enrichment,
            )

            species_names, species_biofluids = collect_species_info(
                species_lookup=species_lookup,
                inchikeys=inchikeys_for_enrichment,
            )

            results.append(
                [
                    "|".join(humanet_ids),
                    "|".join(query_names),
                    "|".join(input_files),
                    "|".join(sorted(set(matched_types))),
                    "|".join(sorted(set(matched_values))),
                    unique_join(query_inchikeys),
                    unique_join(query_pubchems),
                    unique_join(query_hmdbs),
                    hmdb_biofluid,
                    hmdb_diseases,
                    species_names,
                    species_biofluids,
                    row["SMPDB ID"],
                    row["Pathway Name"],
                    row["Pathway Subject"],
                    row["Metabolite ID"],
                ]
            )

    return pd.DataFrame(
        results,
        columns=[
            "HuMANet_ID",
            "Query_Name",
            "Study_Folder",
            "Annotation_ID_Type",
            "Annotation_ID_Value",
            "InChIKey",
            "PubChem_CID",
            "HMDB_ID",
            "HMDB_Biofluid",
            "HMDB_Disease",
            "Species",
            "Species_Biofluids",
            "SMPDB_ID",
            "Pathway_Name",
            "Pathway_Subject",
            "SMPDB_Metabolite_ID",
        ],
    )


def run_stage5(
    input_file: str,
    smpdb_dir: str,
    hmdb_lite_file: str,
    species_file: str,
    output_file: Optional[str] = None,
) -> str:
    output_file = output_file or default_output_path(input_file)

    meta = pd.read_csv(input_file, low_memory=False)
    if "Study_Folder" not in meta.columns and "Input_File" in meta.columns:
        meta["Study_Folder"] = meta["Input_File"]
    meta = meta[[
        "HuMANet_ID",
        "Query_Name",
        "Study_Folder",
        "InChIKey",
        "PubChem_CID",
        "HMDB_ID",
        "KEGG_ID",
        "ChEBI_ID",
    ]].copy()

    out_df = build_stage5_mapping_from_df(
        meta=meta,
        smpdb_dir=smpdb_dir,
        hmdb_lite_file=hmdb_lite_file,
        species_file=species_file,
    )

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    out_df.to_csv(output_file, index=False)

    print("\nFULL MAPPING COMPLETE")
    print("Rows written:", len(out_df))
    print("Saved to:", output_file)
    return output_file


def _explode_pipe_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        if column not in out.columns:
            out[column] = pd.NA
        out[column] = out[column].map(
            lambda value: [item.strip() for item in str(value).split("|") if item.strip()]
            if clean_scalar(value)
            else [pd.NA]
        )
        out = out.explode(column, ignore_index=True)
    return out


def build_stage5_annotation_sheets(mapped_df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    include_drug_annotations = os.getenv("HUMANNET_ENABLE_DRUG_ANNOTATIONS", "0").strip().lower() in {"1", "true", "yes", "y"}
    if mapped_df.empty:
        sheets = {
            "species": pd.DataFrame(columns=SPECIES_SHEET_COLUMNS),
            "disease": pd.DataFrame(columns=DISEASE_SHEET_COLUMNS),
            "pathways": pd.DataFrame(columns=PATHWAY_SHEET_COLUMNS),
        }
        if include_drug_annotations:
            sheets.update(build_drug_annotation_sheets(mapped_df))
        return sheets

    stage5_df = build_stage5_mapping_from_df(
        meta=mapped_df,
        smpdb_dir=SMPDB_DIR,
        hmdb_lite_file=HMDB_LITE_CSV,
        species_file=SPECIES_INFO_CSV,
    )

    hmdb_raw = _load_hmdb_raw_for_sheets(HMDB_LITE_CSV)
    species_raw = _load_species_raw_for_sheets(SPECIES_INFO_CSV)

    disease = _build_disease_sheet(stage5_df, hmdb_raw)
    species = _build_species_sheet(stage5_df, species_raw)
    pathways = stage5_df[PATHWAY_SHEET_COLUMNS].copy()

    disease = disease.drop_duplicates().reset_index(drop=True)
    species = species.drop_duplicates().reset_index(drop=True)
    pathways = pathways.drop_duplicates().reset_index(drop=True)

    sheets = {
        "species": species,
        "disease": disease,
        "pathways": pathways,
    }
    if include_drug_annotations:
        sheets.update(build_drug_annotation_sheets(mapped_df))
    return sheets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map HuMANet confidence-stage metabolites to SMPDB and enrich with HMDB/species annotations."
    )
    parser.add_argument(
        "input_file_positional",
        nargs="?",
        default=None,
        help="Optional positional input CSV path.",
    )
    parser.add_argument(
        "--input-file",
        default=None,
        help="Stage 2 confidence CSV containing HuMANet metadata.",
    )
    parser.add_argument(
        "--output-file",
        default=None,
        help="Optional output CSV path. Defaults to <input>_stage5_smpdb_mapping.csv.",
    )
    args = parser.parse_args()
    if args.input_file_positional:
        args.input_file = args.input_file_positional
    return args


if __name__ == "__main__":
    args = parse_args()
    run_stage5(
        input_file=args.input_file,
        smpdb_dir=SMPDB_DIR,
        hmdb_lite_file=HMDB_LITE_CSV,
        species_file=SPECIES_INFO_CSV,
        output_file=args.output_file,
    )
