import argparse
import gzip
import re
import sqlite3
import time
from pathlib import Path

from typing import Iterable, List, Tuple, Optional


RESOURCE_DIR = Path(__file__).resolve().parent
DATABASE_DIR = RESOURCE_DIR / "Databases"
DEFAULT_DOWNLOAD_DIR = DATABASE_DIR / "pubchem_source_gz"
DEFAULT_OUT_DB = DATABASE_DIR / "pubchem_offline.sqlite"


def aggressive_name_key(x: object) -> Optional[str]:
    if x is None:
        return None
    s = str(x).strip()
    if not s:
        return None
    s = re.sub(r"[‐‑‒–—−]", "-", s)
    s = s.replace('"', "").replace("*", "").replace("&", "")
    s = re.sub(r"\s+", " ", s).strip().lower()
    s = re.sub(r"[^a-z0-9]", "", s)
    return s or None


def iter_two_cols_gz(gz_path: str, sep: Optional[str] = None) -> Iterable[Tuple[str, str]]:
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split(sep, 1) if sep is not None else line.split(None, 1)
            if len(parts) == 2:
                yield parts[0], parts[1]


def iter_cid_inchi_key_gz(gz_path: str) -> Iterable[Tuple[str, str]]:
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                cid = parts[0]
                inchikey = parts[-1]
                if cid and inchikey:
                    yield cid, inchikey


def load_cid_inchi_key_gz(
    cur: sqlite3.Cursor,
    gz_path: str,
    insert_sql: str,
    batch_size: int = 250_000,
    progress_interval: int = 1_000_000,
):
    start = time.time()
    batch: List[Tuple[str, str]] = []
    total = 0

    print(f"\nLoading {gz_path}")

    for row in iter_cid_inchi_key_gz(gz_path):
        batch.append(row)

        if len(batch) >= batch_size:
            cur.executemany(insert_sql, batch)
            total += len(batch)
            batch.clear()

            if total % progress_interval < batch_size:
                elapsed = time.time() - start
                rate = total / elapsed if elapsed else 0
                print(
                    f"  {total:,} rows inserted "
                    f"({rate:,.0f} rows/sec, {elapsed/60:.1f} min elapsed)"
                )

    if batch:
        cur.executemany(insert_sql, batch)
        total += len(batch)

    elapsed = time.time() - start
    print(
        f"Finished {gz_path} -> {total:,} rows "
        f"in {elapsed/60:.2f} minutes "
        f"({total/elapsed:,.0f} rows/sec)"
    )


def load_two_col_gz(
    cur: sqlite3.Cursor,
    gz_path: str,
    insert_sql: str,
    sep: Optional[str] = None,
    batch_size: int = 250_000,
    progress_interval: int = 1_000_000,
):
    start = time.time()
    batch: List[Tuple[str, str]] = []
    total = 0

    print(f"\nLoading {gz_path}")

    for row in iter_two_cols_gz(gz_path, sep=sep):
        batch.append(row)

        if len(batch) >= batch_size:
            cur.executemany(insert_sql, batch)
            total += len(batch)
            batch.clear()

            if total % progress_interval < batch_size:
                elapsed = time.time() - start
                rate = total / elapsed if elapsed else 0
                print(
                    f"  {total:,} rows inserted "
                    f"({rate:,.0f} rows/sec, {elapsed/60:.1f} min elapsed)"
                )

    if batch:
        cur.executemany(insert_sql, batch)
        total += len(batch)

    elapsed = time.time() - start
    print(
        f"Finished {gz_path} → {total:,} rows "
        f"in {elapsed/60:.2f} minutes "
        f"({total/elapsed:,.0f} rows/sec)"
    )


def build_cid_syn_aggr(
    con: sqlite3.Connection,
    batch_size: int = 250_000,
    progress_interval: int = 1_000_000,
):
    start = time.time()
    batch: List[Tuple[str, str, str]] = []
    total = 0
    read_cur = con.cursor()
    write_cur = con.cursor()

    print("\nBuilding cid_syn_aggr from cid_syn")

    for cid, syn in read_cur.execute("SELECT cid, syn FROM cid_syn"):
        key = aggressive_name_key(syn)
        if not key:
            continue
        batch.append((key, cid, syn))

        if len(batch) >= batch_size:
            write_cur.executemany("INSERT INTO cid_syn_aggr VALUES (?,?,?)", batch)
            total += len(batch)
            batch.clear()

            if total % progress_interval < batch_size:
                elapsed = time.time() - start
                rate = total / elapsed if elapsed else 0
                print(
                    f"  {total:,} rows inserted "
                    f"({rate:,.0f} rows/sec, {elapsed/60:.1f} min elapsed)"
                )

    if batch:
        write_cur.executemany("INSERT INTO cid_syn_aggr VALUES (?,?,?)", batch)
        total += len(batch)

    elapsed = time.time() - start
    print(
        f"Finished cid_syn_aggr -> {total:,} rows "
        f"in {elapsed/60:.2f} minutes "
        f"({total/elapsed:,.0f} rows/sec)"
    )


def main():

    start_total = time.time()

    ap = argparse.ArgumentParser(description="Build offline PubChem SQLite DB")
    ap.add_argument("--download_dir", default=str(DEFAULT_DOWNLOAD_DIR))
    ap.add_argument("--cid_smiles_gz", default=None)
    ap.add_argument("--cid_syn_gz", default=None)
    ap.add_argument("--cid_inchi_gz", default=None)
    ap.add_argument("--cid_title_gz", default=None)
    ap.add_argument("--cid_iupac_gz", default=None)
    ap.add_argument("--out_db", default=str(DEFAULT_OUT_DB))
    ap.add_argument("--batch_size", type=int, default=250_000)
    ap.add_argument("--tmp_build", action="store_true")

    args = ap.parse_args()

    download_dir = Path(args.download_dir)
    args.cid_smiles_gz = args.cid_smiles_gz or str(download_dir / "CID-SMILES.gz")
    args.cid_syn_gz = args.cid_syn_gz or str(download_dir / "CID-Synonym-unfiltered.gz")
    args.cid_inchi_gz = args.cid_inchi_gz or str(download_dir / "CID-InChI-Key.gz")
    args.cid_title_gz = args.cid_title_gz or str(download_dir / "CID-Title.gz")
    args.cid_iupac_gz = args.cid_iupac_gz or str(download_dir / "CID-IUPAC.gz")

    inputs = [
        args.cid_smiles_gz,
        args.cid_syn_gz,
        args.cid_inchi_gz,
        args.cid_title_gz,
        args.cid_iupac_gz,
    ]

    for f in inputs:
        if not Path(f).exists():
            raise FileNotFoundError(f"Missing {f}")

    out_db = Path(args.out_db)
    out_db.parent.mkdir(parents=True, exist_ok=True)
    build_db = out_db.with_suffix(out_db.suffix + ".tmp") if args.tmp_build else out_db

    for p in [build_db, Path(str(build_db) + "-wal"), Path(str(build_db) + "-shm")]:
        if p.exists():
            p.unlink()

    print("\nOpening database:", build_db)

    con = sqlite3.connect(str(build_db), timeout=120)
    con.execute("PRAGMA busy_timeout=120000")

    print("Applying performance PRAGMAs")

    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA temp_store=MEMORY")
    con.execute("PRAGMA cache_size=-500000")
    con.execute("PRAGMA locking_mode=EXCLUSIVE")

    cur = con.cursor()

    print("\nCreating tables")

    with con:

        cur.executescript("""
            DROP TABLE IF EXISTS cid_smiles;
            DROP TABLE IF EXISTS cid_syn;
            DROP TABLE IF EXISTS cid_syn_aggr;
            DROP TABLE IF EXISTS cid_inchi;
            DROP TABLE IF EXISTS cid_title;
            DROP TABLE IF EXISTS cid_iupac;

            CREATE TABLE cid_smiles (cid TEXT PRIMARY KEY, smiles TEXT);
            CREATE TABLE cid_syn    (cid TEXT, syn TEXT);
            CREATE TABLE cid_syn_aggr (syn_aggr TEXT, cid TEXT, syn TEXT);
            CREATE TABLE cid_inchi  (cid TEXT PRIMARY KEY, inchikey TEXT);
            CREATE TABLE cid_title  (cid TEXT PRIMARY KEY, title TEXT);
            CREATE TABLE cid_iupac  (cid TEXT PRIMARY KEY, iupac_name TEXT);
        """)

        load_two_col_gz(
            cur,
            args.cid_smiles_gz,
            "INSERT OR REPLACE INTO cid_smiles VALUES (?,?)",
            batch_size=args.batch_size,
        )

        load_two_col_gz(
            cur,
            args.cid_syn_gz,
            "INSERT INTO cid_syn VALUES (?,?)",
            sep="\t",
            batch_size=args.batch_size,
        )

        build_cid_syn_aggr(
            con,
            batch_size=args.batch_size,
        )

        load_cid_inchi_key_gz(
            cur,
            args.cid_inchi_gz,
            "INSERT OR REPLACE INTO cid_inchi VALUES (?,?)",
            batch_size=args.batch_size,
        )

        load_two_col_gz(
            cur,
            args.cid_title_gz,
            "INSERT OR REPLACE INTO cid_title VALUES (?,?)",
            batch_size=args.batch_size,
        )

        load_two_col_gz(
            cur,
            args.cid_iupac_gz,
            "INSERT OR REPLACE INTO cid_iupac VALUES (?,?)",
            batch_size=args.batch_size,
        )

        print("\nBuilding indexes")

        idx_start = time.time()

        cur.executescript("""
            CREATE INDEX idx_syn_cid     ON cid_syn(cid);
            CREATE INDEX idx_syn_lower   ON cid_syn(lower(syn));
            CREATE INDEX idx_cid_syn_aggr_key ON cid_syn_aggr(syn_aggr);
            CREATE INDEX idx_cid_syn_aggr_cid ON cid_syn_aggr(cid);
            CREATE INDEX idx_inchi_cid   ON cid_inchi(cid);
            CREATE INDEX idx_inchi_upper ON cid_inchi(upper(inchikey));
            CREATE INDEX idx_title_cid   ON cid_title(cid);
            CREATE INDEX idx_title_lower ON cid_title(lower(title));
            CREATE INDEX idx_iupac_cid   ON cid_iupac(cid);
            CREATE INDEX idx_iupac_lower ON cid_iupac(lower(iupac_name));
        """)


        print(f"Indexes built in {(time.time()-idx_start)/60:.2f} minutes")

    con.close()

    if args.tmp_build:
        build_db.replace(out_db)

    print("\nDatabase built:", out_db)
    print(f"Total runtime: {(time.time()-start_total)/60:.2f} minutes")


if __name__ == "__main__":
    main()
