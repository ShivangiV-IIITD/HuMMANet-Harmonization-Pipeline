import argparse
import subprocess
import sys
from pathlib import Path


RESOURCE_DIR = Path(__file__).resolve().parent
DEFAULT_DOWNLOAD_SCRIPT = RESOURCE_DIR / "download_pubchem_source_files.py"
DEFAULT_BUILD_SCRIPT = RESOURCE_DIR / "build_pubchem_sqlite.py"
DEFAULT_DOWNLOAD_DIR = RESOURCE_DIR / "Databases" / "pubchem_source_gz"
DEFAULT_OUT_DB = RESOURCE_DIR / "Databases" / "pubchem_offline.sqlite"


def _run(cmd):
    print(f"[BOOTSTRAP] Running: {' '.join(str(x) for x in cmd)}", flush=True)
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="One-step bootstrap for downloading PubChem source files and building the offline SQLite database."
    )
    parser.add_argument("--python", default=sys.executable, help="Python executable to use for child scripts.")
    parser.add_argument("--download_dir", default=str(DEFAULT_DOWNLOAD_DIR), help="Directory for downloaded PubChem .gz files.")
    parser.add_argument("--out_db", default=str(DEFAULT_OUT_DB), help="Output SQLite database path.")
    parser.add_argument("--batch_size", type=int, default=250_000, help="Batch size for SQLite build inserts.")
    parser.add_argument("--tmp_build", action="store_true", help="Build SQLite into a temporary file before replacing the final DB.")
    parser.add_argument("--force_download", action="store_true", help="Redownload PubChem source files even if they already exist.")
    args = parser.parse_args()

    download_cmd = [
        args.python,
        str(DEFAULT_DOWNLOAD_SCRIPT),
        "--out_dir",
        args.download_dir,
    ]
    if args.force_download:
        download_cmd.append("--force")

    build_cmd = [
        args.python,
        str(DEFAULT_BUILD_SCRIPT),
        "--download_dir",
        args.download_dir,
        "--out_db",
        args.out_db,
        "--batch_size",
        str(args.batch_size),
    ]
    if args.tmp_build:
        build_cmd.append("--tmp_build")

    print("[BOOTSTRAP] Step 1/2: download PubChem source files", flush=True)
    _run(download_cmd)

    print("[BOOTSTRAP] Step 2/2: build offline PubChem SQLite", flush=True)
    _run(build_cmd)

    print("[BOOTSTRAP] Done.", flush=True)
    print(f"[BOOTSTRAP] SQLite DB: {args.out_db}", flush=True)


if __name__ == "__main__":
    main()
