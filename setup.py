from setuptools import setup


setup(
    name="humanet-final-pipeline",
    version="0.1.0",
    description="HuMANet metabolite harmonization pipeline with local PubChem SQLite bootstrap workflow",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    python_requires=">=3.10,<3.11",
    py_modules=[
        "bootstrap_pubchem_sqlite",
        "build_pubchem_sqlite",
        "download_pubchem_source_files",
        "drug_similarity_step",
        "extended_annotation_step",
        "knowledge_annotation_step",
        "library_lookup_step",
        "pipeline",
        "pipeline_utils",
        "pubchem_hmdb_reconciliation_step",
        "r_environment_bootstrap",
        "refmet_harmonization_step",
        "resource_config",
        "semicolon_fuzzy_mapper",
    ],
    entry_points={
        "console_scripts": [
            "humanet-pipeline=pipeline:main",
            "humanet-bootstrap-pubchem-sqlite=bootstrap_pubchem_sqlite:main",
            "humanet-download-pubchem-source-files=download_pubchem_source_files:main",
            "humanet-build-pubchem-sqlite=build_pubchem_sqlite:main",
            "humanet-semicolon-stage3-mapper=semicolon_fuzzy_mapper:main",
        ]
    },
)
