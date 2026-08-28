import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path


RESOURCE_DIR = Path(__file__).resolve().parent
DATABASE_DIR = RESOURCE_DIR / "Databases"
DEFAULT_OUT_DIR = DATABASE_DIR / "pubchem_source_gz"
PUBCHEM_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras"

PUBCHEM_FILES = {
    "CID-SMILES.gz": f"{PUBCHEM_BASE_URL}/CID-SMILES.gz",
    "CID-Synonym-unfiltered.gz": f"{PUBCHEM_BASE_URL}/CID-Synonym-unfiltered.gz",
    "CID-InChI-Key.gz": f"{PUBCHEM_BASE_URL}/CID-InChI-Key.gz",
    "CID-Title.gz": f"{PUBCHEM_BASE_URL}/CID-Title.gz",
    "CID-IUPAC.gz": f"{PUBCHEM_BASE_URL}/CID-IUPAC.gz",
}


def _format_size(num_bytes: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(num_bytes)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.1f}{unit}"
        size /= 1024
    return f"{num_bytes}B"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, destination: Path, force: bool = False) -> None:
    tmp_path = destination.with_suffix(destination.suffix + ".part")
    if destination.exists() and not force:
        print(f"[SKIP] {destination.name} already exists: {destination}")
        return

    if tmp_path.exists():
        tmp_path.unlink()

    print(f"[DOWNLOAD] {destination.name}")
    print(f"  URL: {url}")
    print(f"  OUT: {destination}")

    with urllib.request.urlopen(url) as response, tmp_path.open("wb") as out:
        total = response.headers.get("Content-Length")
        total_bytes = int(total) if total else None
        downloaded = 0
        next_report = 0

        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            downloaded += len(chunk)

            if total_bytes:
                pct = (downloaded / total_bytes) * 100
                if pct >= next_report:
                    print(
                        f"  Progress: {pct:5.1f}% "
                        f"({_format_size(downloaded)} / {_format_size(total_bytes)})",
                        flush=True,
                    )
                    next_report += 5
            elif downloaded // (1024 * 1024 * 500) > (downloaded - len(chunk)) // (1024 * 1024 * 500):
                print(f"  Downloaded {_format_size(downloaded)}", flush=True)

    tmp_path.replace(destination)
    print(f"  Done: {_format_size(destination.stat().st_size)}")
    print(f"  SHA256: {_sha256(destination)}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Download PubChem source .gz files needed to build the offline SQLite database.")
    parser.add_argument("--out_dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--force", action="store_true", help="Redownload files even if they already exist.")
    parser.add_argument(
        "--files",
        nargs="*",
        default=list(PUBCHEM_FILES.keys()),
        help="Subset of filenames to download. Default: all required files.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    requested = []
    for name in args.files:
        if name not in PUBCHEM_FILES:
            print(f"[ERROR] Unsupported file name: {name}", file=sys.stderr)
            print("Supported files:", ", ".join(PUBCHEM_FILES.keys()), file=sys.stderr)
            sys.exit(1)
        requested.append(name)

    print(f"[INFO] Output directory: {out_dir}")
    print(f"[INFO] Files to download: {', '.join(requested)}")

    for name in requested:
        _download(PUBCHEM_FILES[name], out_dir / name, force=args.force)

    print("[DONE] All requested PubChem source files are available.")


if __name__ == "__main__":
    main()
