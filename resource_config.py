import os
from pathlib import Path


RESOURCE_DIR = Path(__file__).resolve().parent
DATABASE_DIR = Path(os.getenv("HUMANNET_RESOURCE_DIR", str(RESOURCE_DIR / "Databases"))).resolve()

HMDB_LITE_CSV = str(DATABASE_DIR / "HMDB_Lite.csv")
HUMANNET_LIBRARY_CSV = str(DATABASE_DIR / "Humannet_Library_V1_ungrouped.csv")
MW_DATABASE_CSV = str(DATABASE_DIR / "MW_Database.csv")
MW_UNMAPPED_DATABASE_CSV = str(DATABASE_DIR / "MW_unmapped_Database.csv")
PUBCHEM_OFFLINE_SQLITE = str(DATABASE_DIR / "pubchem_offline.sqlite")
SPECIES_INFO_CSV = str(DATABASE_DIR / "Metabolites_associated_with_species_all_info.csv")
SMPDB_DIR = str(DATABASE_DIR / "smpdb_metabolites")
SMPDB_MAPPING_CACHE = str(DATABASE_DIR / "smpdb_metabolite_mapping.pkl")
DRUG_ANNOTATION_DIR = DATABASE_DIR / "drug_annotation"
DRUGBANK_SDF = str(DRUG_ANNOTATION_DIR / "structures.sdf")
DRUGBANK_XML = str(DRUG_ANNOTATION_DIR / "full_database.xml")
DRUGCENTRAL_TSV = str(DRUG_ANNOTATION_DIR / "drugcentral.tsv")

CLASSYFIRE_BRIDGE_R = str(RESOURCE_DIR / "stage2_classyfire_bridge.R")
REFMET_BRIDGE_R = str(RESOURCE_DIR / "refmet_bridge.R")
