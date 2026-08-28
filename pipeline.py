import argparse
import os
import sys
from typing import List, Optional

import pandas as pd

from extended_annotation_step import run_stage3
from library_lookup_step import run_library_lookup
from pubchem_hmdb_reconciliation_step import run_stage2
from refmet_harmonization_step import run_stage1
from resource_config import (
    CLASSYFIRE_BRIDGE_R,
    HMDB_LITE_CSV,
    HUMANNET_LIBRARY_CSV,
    MW_DATABASE_CSV,
    MW_UNMAPPED_DATABASE_CSV,
    PUBCHEM_OFFLINE_SQLITE,
)
from pipeline_utils import HAS_RPY2
from semicolon_fuzzy_mapper import run_mapper as run_fuzzy_mapper


VALID_STAGES = {"1", "2", "3", "4"}


def parse_stage_selection(raw: str) -> List[str]:
    stages = []
    for part in raw.split(","):
        value = part.strip()
        if value and value not in stages:
            stages.append(value)

    invalid = [stage for stage in stages if stage not in VALID_STAGES]
    if not stages or invalid:
        raise ValueError("Choose stages using only 1, 2, 3, 4 or comma-separated combinations like 1,2 or 1,2,3,4.")
    return stages


def prompt_nonempty(message: str) -> str:
    while True:
        value = input(message).strip()
        if value:
            return value
        print("A value is required.")


def prompt_optional(message: str, default: str) -> str:
    value = input(f"{message} [{default}]: ").strip()
    return value or default


def prompt_yes_no(message: str, default: bool = False) -> bool:
    suffix = "Y/n" if default else "y/N"
    value = input(f"{message} [{suffix}]: ").strip().lower()
    if not value:
        return default
    return value in {"y", "yes"}


def is_interactive_stdin() -> bool:
    try:
        return sys.stdin.isatty()
    except Exception:
        return False


def default_prefix_for_input(input_path: str, stage_name: str) -> str:
    base_dir = os.path.dirname(os.path.abspath(input_path))
    stem = os.path.splitext(os.path.basename(input_path))[0]
    return os.path.join(base_dir, f"{stem}_{stage_name}")


def _canonical_original_query(df: pd.DataFrame) -> pd.Series:
    if "original_query_name" in df.columns:
        out = df["original_query_name"].astype("string")
    elif "Query_Name" in df.columns:
        out = df["Query_Name"].astype("string")
    else:
        out = pd.Series([pd.NA] * len(df), dtype="string")
    return out.astype("string").str.strip()


def build_final_unmapped(
    library_prefix: Optional[str],
    stage2_prefix: Optional[str],
    stage3_prefix: Optional[str],
    stage4_prefix: str,
) -> str:
    mapped_paths = []
    if library_prefix:
        mapped_paths.append(f"{library_prefix}_mapped.csv")
    if stage2_prefix:
        mapped_paths.append(f"{stage2_prefix}_mapped.csv")
    if stage3_prefix:
        mapped_paths.extend(
            [
                f"{stage3_prefix}_confidence.csv",
                f"{stage3_prefix}_nonconfidence_expanded.csv",
            ]
        )
    mapped_paths.append(f"{stage4_prefix}_mapped.csv")
    final_unmapped_path = f"{stage4_prefix}_final_unmapped.csv"
    stage4_unmapped_path = f"{stage4_prefix}_unmapped.csv"

    mapped_originals = set()
    for path in mapped_paths:
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, dtype="string")
        mapped_originals.update(
            x for x in _canonical_original_query(df).dropna().tolist() if x
        )

    if not os.path.exists(stage4_unmapped_path):
        raise FileNotFoundError(f"Stage 4 unmapped file not found: {stage4_unmapped_path}")

    stage4_unmapped = pd.read_csv(stage4_unmapped_path, dtype="string")
    keep_mask = ~_canonical_original_query(stage4_unmapped).isin(mapped_originals)
    final_unmapped = stage4_unmapped.loc[keep_mask].copy()
    final_unmapped.to_csv(final_unmapped_path, index=False)
    return final_unmapped_path


def _has_primary_mapping_identifier(df: pd.DataFrame) -> pd.Series:
    present = [c for c in ["PubChem_CID", "HMDB_ID", "KEGG_ID"] if c in df.columns]
    if not present:
        return pd.Series([False] * len(df), index=df.index)
    out = df[present[0]].astype("string").notna()
    for col in present[1:]:
        out = out | df[col].astype("string").notna()
    return out


def run_fuzzy_matching(input_path: str, out_prefix: str, workers: int = 8) -> str:
    all_path = f"{out_prefix}_all.csv"
    mapped_path = f"{out_prefix}_mapped.csv"
    unmapped_path = f"{out_prefix}_unmapped.csv"

    out = run_fuzzy_mapper(
        input_csv=input_path,
        output_csv=all_path,
        query_col="Query_Name",
        pubchem_workers=workers,
    )
    mapped = out.loc[_has_primary_mapping_identifier(out)].copy()
    unmapped = out.loc[~_has_primary_mapping_identifier(out)].copy()
    mapped.to_csv(mapped_path, index=False)
    unmapped.to_csv(unmapped_path, index=False)
    return unmapped_path


def run_stage1_from_pipeline(input_path: str, out_prefix: str) -> None:
    run_library_lookup(
        input_csv=input_path,
        library_csv=HUMANNET_LIBRARY_CSV,
        out_prefix=out_prefix,
    )


def run_stage2_from_pipeline(input_path: str, out_prefix: str, workers: int = 8) -> None:
    bridge = CLASSYFIRE_BRIDGE_R if HAS_RPY2 else None
    run_stage1(
        stage1_unmapped_csv=input_path,
        hmdb_lite_csv=HMDB_LITE_CSV,
        sqlite_db=PUBCHEM_OFFLINE_SQLITE,
        out_prefix=out_prefix,
        r_classyfire_bridge=bridge,
        sleep=0.3,
        pugview_workers=workers,
        skip_classyfire=False,
    )


def run_stage3_from_pipeline(input_path: str, out_prefix: str, workers: int = 8) -> None:
    bridge = CLASSYFIRE_BRIDGE_R if HAS_RPY2 else None
    run_stage2(
        stage1_unmapped_csv=input_path,
        hmdb_lite_csv=HMDB_LITE_CSV,
        sqlite_db=PUBCHEM_OFFLINE_SQLITE,
        out_prefix=out_prefix,
        r_classyfire_bridge=bridge,
        sleep=0.3,
        pugview_workers=workers,
        skip_classyfire=False,
        print_classyfire_output=False,
    )


def run_stage4_from_pipeline(input_path: str, out_prefix: str, workers: int = 8) -> None:
    bridge = CLASSYFIRE_BRIDGE_R if HAS_RPY2 else None
    run_stage3(
        stage2_unmapped_csv=input_path,
        mw_database_csv=MW_DATABASE_CSV,
        out_prefix=out_prefix,
        hmdb_lite_csv=HMDB_LITE_CSV,
        sqlite_db=PUBCHEM_OFFLINE_SQLITE,
        pubchem_workers=workers,
        r_classyfire_bridge=bridge,
        skip_classyfire=False,
        print_classyfire_output=False,
        sleep=0.2,
        sqlite_name_fallback=True,
        mw_unmapped_db_csv=MW_UNMAPPED_DATABASE_CSV,
    )


def run_all_four_chained(
    workers: int = 8,
    use_fuzzy_matching: bool = False,
    input_csv: Optional[str] = None,
    base_prefix: Optional[str] = None,
) -> None:
    stage1_input = input_csv or prompt_nonempty("Enter the raw input file for Stage 1 library lookup (for example a CSV with Database source and metabolite_name): ")
    base_prefix = base_prefix or prompt_optional("Enter a base output prefix", default_prefix_for_input(stage1_input, "pipeline"))

    stage1_prefix = f"{base_prefix}_stage1_library"
    print(f"\nRunning Stage 1 library lookup with input: {stage1_input}")
    run_stage1_from_pipeline(stage1_input, stage1_prefix)

    stage2_input = f"{stage1_prefix}_unmapped.csv"
    stage2_prefix = f"{base_prefix}_stage2"
    print(f"\nRunning Stage 2 with input: {stage2_input}")
    run_stage2_from_pipeline(stage2_input, stage2_prefix, workers=workers)

    stage3_input = f"{stage2_prefix}_unmapped.csv"
    stage3_prefix = f"{base_prefix}_stage3"
    print(f"\nRunning Stage 3 with input: {stage3_input}")
    run_stage3_from_pipeline(stage3_input, stage3_prefix, workers=workers)

    stage4_input = f"{stage3_prefix}_unmapped.csv"
    stage4_prefix = f"{base_prefix}_stage4"
    print(f"\nRunning Stage 4 with input: {stage4_input}")
    run_stage4_from_pipeline(stage4_input, stage4_prefix, workers=workers)

    final_unmapped_path = build_final_unmapped(stage1_prefix, stage2_prefix, stage3_prefix, stage4_prefix)
    fuzzy_unmapped_path = None
    if use_fuzzy_matching:
        fuzzy_prefix = f"{base_prefix}_fuzzy"
        print(f"\nRunning semicolon-based alternative-name matching with input: {final_unmapped_path}")
        fuzzy_unmapped_path = run_fuzzy_matching(final_unmapped_path, fuzzy_prefix, workers=workers)

    print("\nPipeline complete.")
    print(f"Stage 1 library outputs: {stage1_prefix}_*.csv and {stage1_prefix}.xlsx")
    print(f"Stage 2 outputs: {stage2_prefix}_*.csv and {stage2_prefix}.xlsx")
    print(f"Stage 3 outputs: {stage3_prefix}_*.csv and {stage3_prefix}.xlsx")
    print(f"Stage 4 outputs: {stage4_prefix}_*.csv and {stage4_prefix}.xlsx")
    print(f"Final unmapped output: {final_unmapped_path}")
    if fuzzy_unmapped_path:
        print(f"Fuzzy matching outputs: {fuzzy_prefix}_*.csv")
        print(f"Final unmapped after semicolon-based alternative-name matching: {fuzzy_unmapped_path}")


def run_selected_stages_independently(stages: List[str], workers: int = 8, use_fuzzy_matching: bool = False) -> None:
    stage1_prefix_run: Optional[str] = None
    stage2_prefix_run: Optional[str] = None
    stage3_prefix_run: Optional[str] = None
    stage4_prefix_run: Optional[str] = None

    for stage in stages:
        if stage == "1":
            input_path = prompt_nonempty("\nEnter the raw input file for Stage 1 library lookup: ")
            out_prefix = prompt_optional("Enter the output prefix for Stage 1 library lookup", default_prefix_for_input(input_path, "stage1_library"))
            print(f"Running Stage 1 library lookup with input: {input_path}")
            run_stage1_from_pipeline(input_path, out_prefix)
            stage1_prefix_run = out_prefix
        elif stage == "2":
            input_path = prompt_nonempty("\nEnter the Stage 1 library-unmapped CSV for Stage 2: ")
            out_prefix = prompt_optional("Enter the output prefix for Stage 2", default_prefix_for_input(input_path, "stage2"))
            print(f"Running Stage 2 with input: {input_path}")
            run_stage2_from_pipeline(input_path, out_prefix, workers=workers)
            stage2_prefix_run = out_prefix
        elif stage == "3":
            input_path = prompt_nonempty("\nEnter the Stage 2 unmapped CSV for Stage 3: ")
            out_prefix = prompt_optional("Enter the output prefix for Stage 3", default_prefix_for_input(input_path, "stage3"))
            print(f"Running Stage 3 with input: {input_path}")
            run_stage3_from_pipeline(input_path, out_prefix, workers=workers)
            stage3_prefix_run = out_prefix
        elif stage == "4":
            input_path = prompt_nonempty("\nEnter the Stage 3 unmapped CSV for Stage 4: ")
            out_prefix = prompt_optional("Enter the output prefix for Stage 4", default_prefix_for_input(input_path, "stage4"))
            print(f"Running Stage 4 with input: {input_path}")
            run_stage4_from_pipeline(input_path, out_prefix, workers=workers)
            stage4_prefix_run = out_prefix

    if stage4_prefix_run:
        final_unmapped_path = build_final_unmapped(stage1_prefix_run, stage2_prefix_run, stage3_prefix_run, stage4_prefix_run)
        print(f"\nFinal unmapped output: {final_unmapped_path}")
        if use_fuzzy_matching:
            fuzzy_prefix = f"{stage4_prefix_run}_fuzzy"
            fuzzy_unmapped_path = run_fuzzy_matching(final_unmapped_path, fuzzy_prefix, workers=workers)
            print(f"Fuzzy matching outputs: {fuzzy_prefix}_*.csv")
            print(f"Final unmapped after semicolon-based alternative-name matching: {fuzzy_unmapped_path}")

    print("\nSelected stage runs complete.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Interactive HuMANet stage pipeline runner. Stage 1 performs the HuMANet library lookup before the original Stage 1-3 flow."
    )
    parser.add_argument(
        "--stages",
        default=None,
        help="Optional stage selection like 1, 2, 3, 4, 1,2, or 1,2,3,4. If omitted, you will be prompted.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=min(32, os.cpu_count() or 8),
        help="Worker count to pass through to the pipeline stages.",
    )
    parser.add_argument(
        "--use_fuzzy_matching",
        action="store_true",
        help="Run semicolon-based alternative-name matching on the final unmapped file after Stage 4.",
    )
    parser.add_argument(
        "--input_csv",
        default=None,
        help="Non-interactive raw input CSV for the full chained 1,2,3,4 run.",
    )
    parser.add_argument(
        "--base_prefix",
        default=None,
        help="Non-interactive base output prefix for the full chained 1,2,3,4 run.",
    )
    args = parser.parse_args()

    if args.stages:
        raw_selection = args.stages
    else:
        raw_selection = prompt_nonempty("Which stages do you want to run? Choose from 1, 2, 3, 4 or combinations like 1,2 or 1,2,3,4: ")

    stages = parse_stage_selection(raw_selection)
    use_fuzzy_matching = False
    if "4" in stages:
        if args.use_fuzzy_matching:
            use_fuzzy_matching = True
        elif is_interactive_stdin():
            use_fuzzy_matching = prompt_yes_no(
                "Do you want to run semicolon-based alternative-name matching after Stage 4?",
                default=False,
            )
        else:
            use_fuzzy_matching = False
    if stages == ["1", "2", "3", "4"]:
        run_all_four_chained(
            workers=args.workers,
            use_fuzzy_matching=use_fuzzy_matching,
            input_csv=args.input_csv,
            base_prefix=args.base_prefix,
        )
    else:
        run_selected_stages_independently(stages, workers=args.workers, use_fuzzy_matching=use_fuzzy_matching)


if __name__ == "__main__":
    main()
