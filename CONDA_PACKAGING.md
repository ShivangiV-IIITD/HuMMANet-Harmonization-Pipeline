# Conda Packaging

This folder can be built as a conda package without shipping the local PubChem SQLite database.

## Included in the package

- pipeline scripts
- R bridge files
- guide/documentation
- packaged CSV / PKL resources in `Databases/`

## Excluded from the package

- `Databases/pubchem_offline.sqlite`
- downloaded PubChem source `.gz` files in `Databases/pubchem_source_gz/`
- Python cache files

Users build the SQLite database locally after installation:

```bash
humanet-bootstrap-pubchem-sqlite
```

## Build environment

Use a compatible conda environment with the dependencies required by the package.

An example environment file is included:

`environment.package.yml`

## Create a packaging env (optional)

```bash
conda env create -f environment.package.yml
```

## Build the conda package

```bash
conda build conda-recipe
```

## Installed commands

- `humanet-pipeline`
- `humanet-bootstrap-pubchem-sqlite`
- `humanet-download-pubchem-source-files`
- `humanet-build-pubchem-sqlite`
- `humanet-semicolon-stage3-mapper`

## Notes

- The recipe currently skips Windows.
- The package depends on Python/R runtime dependencies, but users may still need the R packages required by `RefMet` and optional `classyfireR` functionality in their target environment.
