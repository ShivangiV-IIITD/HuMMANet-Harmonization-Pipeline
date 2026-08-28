# HuMANet Final Pipeline

This repository contains the HuMANet metabolite harmonization pipeline code, the bundled HuMANet library used for the initial lookup stage, and helper scripts to bootstrap the local PubChem SQLite database.

Large reference resources are not stored in this repository. The non-PubChem databases are expected to be downloaded from the Zenodo record `10.5281/zenodo.22146911`, from where users can access the DBs, and placed in the HuMANet resource directory, while the PubChem SQLite database is built locally from the provided download and build scripts.

## Package Layout

The repository is organized as:

* `pipeline.py`: main command-line runner for the full HuMANet workflow
* `library_lookup_step.py`: initial HuMANet library lookup stage
* `refmet_harmonization_step.py`: RefMet-backed name standardization stage
* `pubchem_hmdb_reconciliation_step.py`: PubChem/HMDB cross-database reconciliation stage
* `extended_annotation_step.py`: Metabolomics Workbench plus offline PubChem recovery stage
* `knowledge_annotation_step.py`: downstream knowledge-layer annotation sheets
* `drug_similarity_step.py`: DrugBank and DrugCentral structural similarity annotation
* `download_pubchem_source_files.py`, `build_pubchem_sqlite.py`, `bootstrap_pubchem_sqlite.py`: local PubChem SQLite bootstrap workflow
* `Databases/Humannet_Library_V1_ungrouped.csv`: bundled HuMANet lookup library

Installed command-line entry points:

* `humanet-pipeline`
* `humanet-bootstrap-pubchem-sqlite`
* `humanet-download-pubchem-source-files`
* `humanet-build-pubchem-sqlite`
* `humanet-semicolon-stage3-mapper`

## Resource Setup

HuMANet expects a resource directory containing:

* `Humannet_Library_V1_ungrouped.csv`
* `HMDB_Lite.csv`
* `MW_Database.csv`
* `MW_unmapped_Database.csv`
* `Metabolites_associated_with_species_all_info.csv`
* `smpdb_metabolite_mapping.pkl`
* `drug_annotation/`
* `pubchem_offline.sqlite`

Only `Humannet_Library_V1_ungrouped.csv` is bundled with the repository/package. All other non-PubChem resources should be downloaded from Zenodo `10.5281/zenodo.22146911`, from where users can access the DBs, and placed in a resource directory of your choice.

By default the pipeline looks for resources in:

```text
<repo-or-install-root>/Databases
```

You can point HuMANet to a different resource directory by setting:

```bash
export HUMANNET_RESOURCE_DIR=/path/to/HuMANet_resources
```

The directory should then contain the Zenodo-downloaded resources plus the locally built `pubchem_offline.sqlite`.

## PubChem SQLite Setup

The PubChem SQLite database is built locally and is not distributed in this repository.

1. Download the PubChem source files:

```bash
humanet-download-pubchem-source-files --output-dir /path/to/pubchem_source_gz
```

2. Build the SQLite database:

```bash
humanet-build-pubchem-sqlite \
  --input-dir /path/to/pubchem_source_gz \
  --output-db /path/to/HuMANet_resources/pubchem_offline.sqlite
```

Or run the two-step bootstrap wrapper:

```bash
humanet-bootstrap-pubchem-sqlite \
  --download-dir /path/to/pubchem_source_gz \
  --output-db /path/to/HuMANet_resources/pubchem_offline.sqlite
```

## What The Pipeline Does

Given a CSV of study-specific metabolite names, HuMANet runs:

1. library lookup against `Humannet_Library_V1_ungrouped.csv`
2. RefMet-backed harmonization
3. PubChem/HMDB reconciliation
4. Metabolomics Workbench plus offline PubChem recovery
5. optional downstream knowledge-layer annotation sheets

The initial library lookup resolves duplicate names by the following source priority:

1. `all_studies_mapped`
2. `HuMANet_mimedb_all`
3. `HuMANet_monoculture_all`
4. `all_studies_nonconfidence`

## Main Runtime Command

When the package is installed, the main command is:

```bash
humanet-pipeline
```

The runner supports chained or individual execution of the HuMANet stages and can be used interactively or with explicit file paths and output prefixes.

## Full Runtime Example

```bash
humanet-pipeline \
  --input_csv /path/to/stage1_input.csv \
  --base_prefix /path/to/results/run1 \
  --workers 32 \
  --stages 1,2,3,4 \
  --use_fuzzy_matching
```

If you are running the Python script directly:

```bash
python pipeline.py
```

## Runtime Input Requirements

The runtime input must be a CSV with at least:

* `Query_Name`

Recommended additional columns:

* `HuMANet_ID`
* `Study_Folder`
* `Database_Source`

If `HuMANet_ID` is missing, the preparation utilities will create HuMANet-compatible IDs. If `Study_Folder` is missing, the input filename is used as the source label.

## Outputs

Depending on the selected stages, HuMANet writes:

* mapped CSVs
* confidence CSVs
* non-confidence CSVs
* unmapped CSVs
* Excel workbooks for stage-wise review
* optional species, disease, pathway, DrugBank, and DrugCentral annotation sheets

HuMANet also preserves:

* `original_query_name`: the original study-facing input name
* `matched_name`: the name that led to the successful match when applicable

## Development Install

Recommended:

```bash
conda env create -f environment.package.yml
conda activate humanet_refmet
pip install -e .
```

## Notes

* The repository intentionally excludes large databases from version control.
* The bundled HuMANet library is ungrouped so study-specific aliases are retained for the first lookup stage.
* DrugBank structural similarity uses `drug_annotation/structures.sdf`.
