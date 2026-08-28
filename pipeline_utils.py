import json
import os
import re
import sqlite3
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from threading import Lock
from typing import Dict, List, Optional, Sequence

import pandas as pd
import requests

import r_environment_bootstrap  # noqa: F401

try:
    from rpy2.robjects import conversion, default_converter, pandas2ri
    from rpy2.robjects.packages import STAP
    HAS_RPY2 = True
except Exception:
    HAS_RPY2 = False


_PUGVIEW_ID_CACHE: Dict[str, Dict[str, Optional[str]]] = {}
_PUGVIEW_CACHE_LOCK = Lock()


def _clean_scalar(x: object) -> Optional[str]:
    if x is None or pd.isna(x):
        return None
    s = str(x).strip()
    if s in {"", "NA", "Na", "null", "-", "nan", "NA_character_", "<NA>"}:
        return None
    return s


def _split_synonyms(x: object) -> List[str]:
    s = _clean_scalar(x)
    if not s:
        return []
    return [p.strip() for p in str(s).split(";") if p.strip()]


def _clean_cid(cid: object) -> Optional[str]:
    s = _clean_scalar(cid)
    if not s:
        return None
    return s if re.fullmatch(r"\d+", s) else None


def _clean_hmdb_id(hmdb_id: object) -> Optional[str]:
    s = _clean_scalar(hmdb_id)
    if not s:
        return None
    m = re.search(r"HMDB\d{7}", s.upper())
    return m.group(0) if m else None


def _split_ids(val: object) -> List[str]:
    s = _clean_scalar(val)
    if not s:
        return []
    return [v.strip() for v in s.split(";") if _clean_scalar(v)]


def _normalize_hmdb_field(val: object) -> Optional[str]:
    ids = [_clean_hmdb_id(x) for x in _split_ids(val)]
    ids = [x for x in ids if x]
    if not ids:
        return None
    seen = set()
    out = []
    for x in ids:
        if x in seen:
            continue
        seen.add(x)
        out.append(x)
    return ";".join(out)


def _first_hmdb_id(val: object) -> Optional[str]:
    norm = _normalize_hmdb_field(val)
    if not norm:
        return None
    return norm.split(";", 1)[0]


def _clean_inchikey(val: object) -> Optional[str]:
    s = _clean_scalar(val)
    if not s:
        return None
    m = re.search(r"[A-Z]{14}-[A-Z]{10}-[A-Z0-9]", s.upper())
    return m.group(0) if m else None


def _sqlite_table_exists(db_path: str, table_name: str) -> bool:
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    row = cur.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
        (table_name,),
    ).fetchone()
    con.close()
    return row is not None


def _first_existing_column(cols: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    lookup = {c.lower(): c for c in cols}
    for c in candidates:
        if c.lower() in lookup:
            return lookup[c.lower()]
    return None


def prepare_stage1_input_dataframe(df: pd.DataFrame, input_path: str) -> pd.DataFrame:
    """
    Normalize Stage 1 input so the pipeline accepts either:
    1. canonical HuMANet input with Query_Name/HuMANet_ID
    2. raw metabolite CSVs with Database source/metabolite_name
    """
    out = df.copy()
    query_col = _first_existing_column(out.columns, ["Query_Name", "metabolite_name"])
    if query_col is None:
        raise ValueError(
            'Stage 1 input must contain either "Query_Name" or "metabolite_name".'
        )

    annotation_db_col = _first_existing_column(out.columns, ["Database_Source"])
    input_db_col = _first_existing_column(out.columns, ["Input_Database_Source", "Database source"])
    input_file_value = os.path.splitext(os.path.basename(input_path))[0]

    input_file_series = (
        out["Input_File"].astype("string")
        if "Input_File" in out.columns
        else (
            out["Study_Folder"].astype("string")
            if "Study_Folder" in out.columns
            else pd.Series([input_file_value] * len(out), dtype="string")
        )
    )

    normalized = pd.DataFrame(
        {
            "Query_Name": out[query_col].astype("string"),
            "Input_File": input_file_series,
            "Database_Source": (
                out[annotation_db_col].astype("string")
                if annotation_db_col is not None
                else pd.Series([pd.NA] * len(out), dtype="string")
            ),
            "Input_Database_Source": (
                out[input_db_col].astype("string")
                if input_db_col is not None
                else pd.Series([pd.NA] * len(out), dtype="string")
            ),
        }
    )

    if "HuMANet_ID" in out.columns:
        normalized["HuMANet_ID"] = out["HuMANet_ID"].astype("string")
    else:
        normalized["HuMANet_ID"] = pd.Series(
            [f"HMAN{idx:06d}" for idx in range(1, len(out) + 1)],
            dtype="string",
        )

    optional_columns = [
        "HMDB_ID",
        "PubChem_CID",
        "KEGG_ID",
        "ChEBI_ID",
        "InChIKey",
        "SMILES",
        "Synonyms",
        "Standardized_Name",
        "Molecular_Formula",
        "Exact_Mass",
        "Super_Class",
        "Main_Class",
        "Sub_Class",
        "Study_Folder",
        "original_query_name",
        "matched_name",
    ]
    for column in optional_columns:
        if column in out.columns:
            normalized[column] = out[column].astype("string")
        else:
            normalized[column] = pd.Series([pd.NA] * len(out), dtype="string")

    normalized = normalized.dropna(subset=["Query_Name"]).copy()
    normalized["Query_Name"] = normalized["Query_Name"].astype("string").str.strip()
    normalized = normalized[normalized["Query_Name"].ne("")]
    normalized = normalized.reset_index(drop=True)

    if "HuMANet_ID" not in out.columns:
        normalized["HuMANet_ID"] = pd.Series(
            [f"HMAN{idx:06d}" for idx in range(1, len(normalized) + 1)],
            dtype="string",
        )

    return normalized


def fetch_cids_by_name_sqlite(db_path: str, name: str, limit: int = 50) -> List[str]:
    qname = (_clean_scalar(name) or "").lower()
    if not qname:
        return []
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows = cur.execute(
        "SELECT DISTINCT cid FROM cid_syn WHERE lower(syn)=? LIMIT ?",
        (qname, limit),
    ).fetchall()
    con.close()
    out = [_clean_cid(r[0]) for r in rows]
    return [x for x in out if x is not None]


def extract_ids_from_pugview(cid: str, sleep: float = 0.2) -> Dict[str, Optional[str]]:
    clean_cid = _clean_cid(cid)
    if not clean_cid:
        return {"HMDB_ID": None, "KEGG_ID": None, "ChEBI_ID": None}

    with _PUGVIEW_CACHE_LOCK:
        cached = _PUGVIEW_ID_CACHE.get(clean_cid)
    if cached is not None:
        return dict(cached)

    time.sleep(sleep)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{clean_cid}/JSON"
    try:
        res = requests.get(url, timeout=20)
        if res.status_code != 200:
            out = {"HMDB_ID": None, "KEGG_ID": None, "ChEBI_ID": None}
            with _PUGVIEW_CACHE_LOCK:
                _PUGVIEW_ID_CACHE[clean_cid] = out
            return dict(out)
        payload = res.json()
    except Exception:
        out = {"HMDB_ID": None, "KEGG_ID": None, "ChEBI_ID": None}
        with _PUGVIEW_CACHE_LOCK:
            _PUGVIEW_ID_CACHE[clean_cid] = out
        return dict(out)

    txt = json.dumps(payload.get("Record", {}).get("Section", {}))
    hmdb = sorted(set(re.findall(r"HMDB\d{7}", txt)))
    kegg = sorted(set(re.findall(r"\bC\d{5}\b", txt)))
    chebi = sorted(set(re.findall(r"CHEBI:\d+", txt)))

    out = {
        "HMDB_ID": ";".join(hmdb) if hmdb else None,
        "KEGG_ID": ";".join(kegg) if kegg else None,
        "ChEBI_ID": ";".join(chebi) if chebi else None,
    }
    with _PUGVIEW_CACHE_LOCK:
        _PUGVIEW_ID_CACHE[clean_cid] = out
    return dict(out)


def fetch_pugview_ids_for_cids(cids: Sequence[str], sleep: float = 0.2, max_workers: int = 8) -> pd.DataFrame:
    clean_unique = sorted({c for c in (_clean_cid(cid) for cid in cids) if c})
    if not clean_unique:
        return pd.DataFrame(columns=["PubChem_CID", "HMDB_ID", "KEGG_ID", "ChEBI_ID"], dtype="string")

    rows: List[Dict[str, Optional[str]]] = []
    if max_workers <= 1 or len(clean_unique) == 1:
        for cid in clean_unique:
            rec = extract_ids_from_pugview(cid, sleep=sleep)
            rec["PubChem_CID"] = cid
            rows.append(rec)
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            fut_map = {ex.submit(extract_ids_from_pugview, cid, sleep): cid for cid in clean_unique}
            for fut in as_completed(fut_map):
                cid = fut_map[fut]
                try:
                    rec = fut.result()
                except Exception:
                    rec = {"HMDB_ID": None, "KEGG_ID": None, "ChEBI_ID": None}
                rec["PubChem_CID"] = cid
                rows.append(rec)

    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["PubChem_CID", "HMDB_ID", "KEGG_ID", "ChEBI_ID"], dtype="string")
    return out.astype("string")


def load_classyfire_bridge(r_bridge_path: str):
    if not HAS_RPY2:
        return None
    with open(r_bridge_path, "r", encoding="utf-8") as f:
        r_code = f.read()
    mod = STAP(r_code, "stage2_classyfire_bridge")
    return mod.annotate_inchikeys_bridge


def load_refmet_bridge(r_bridge_path: str):
    if not HAS_RPY2:
        return None
    with open(r_bridge_path, "r", encoding="utf-8") as f:
        r_code = f.read()
    mod = STAP(r_code, "refmet_bridge")
    return mod.refmet_map_bridge


def map_refmet_queries(query_names: Sequence[str], r_bridge_path: Optional[str]) -> pd.DataFrame:
    empty = pd.DataFrame()
    if not query_names or not r_bridge_path:
        return empty

    fn = load_refmet_bridge(r_bridge_path)
    if fn is None:
        raise RuntimeError("RefMet bridge could not be loaded.")

    preview = list(query_names)[:10]
    print(f"[REFMET] Query count: {len(query_names)}")
    print(f"[REFMET] Query preview: {preview}")

    with conversion.localconverter(default_converter + pandas2ri.converter):
        r_vec = conversion.py2rpy(pd.Series(list(query_names), dtype="string"))

    r_df = fn(r_vec)
    with conversion.localconverter(default_converter + pandas2ri.converter):
        out = conversion.rpy2py(r_df)
    if not isinstance(out, pd.DataFrame):
        out = pd.DataFrame(out)
    out = out.astype("string")

    print(f"[REFMET] Returned columns: {list(out.columns)}")
    print(f"[REFMET] Returned rows: {len(out)}")
    if not out.empty:
        print(out.head(10).to_string(index=False))

    return out


def annotate_classyfire(
    inchikeys: Sequence[str],
    r_bridge_path: Optional[str],
    print_output: bool = False,
) -> pd.DataFrame:
    cleaned = pd.Series(inchikeys, dtype="string").map(_clean_inchikey)
    uniq = [x for x in cleaned.dropna().unique().tolist() if x]
    empty = pd.DataFrame(columns=["InChIKey", "Super_Class", "Main_Class", "Sub_Class"], dtype="string")

    if not uniq or not r_bridge_path:
        return empty

    try:
        fn = load_classyfire_bridge(r_bridge_path)
        if fn is None:
            return empty

        with conversion.localconverter(default_converter + pandas2ri.converter):
            r_vec = conversion.py2rpy(pd.Series(uniq, dtype="string"))

        r_df = fn(r_vec)
        with conversion.localconverter(default_converter + pandas2ri.converter):
            out = conversion.rpy2py(r_df)
        if not isinstance(out, pd.DataFrame):
            out = pd.DataFrame(out)
        out = out.astype("string")
        if print_output:
            print(f"[CLASSYFIRE] Returned rows: {len(out)} for keys: {len(uniq)}")
            print(out.head(20).to_string(index=False))
        return out
    except Exception as e:
        print(f"[WARN] ClassyFire annotation unavailable, continuing without it: {e}")
        return empty


@contextmanager
def timed_step(label: str):
    start = time.perf_counter()
    print(f"[TIMER] START: {label}")
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        print(f"[TIMER] END: {label} | {elapsed:.2f}s")


def _safe_sheet_name(name: str) -> str:
    cleaned = re.sub(r"[:\\\\/?*\\[\\]]", "_", str(name)).strip()
    return (cleaned or "Sheet")[:31]


def write_excel_workbook(sheets: Dict[str, pd.DataFrame], output_path: str) -> str:
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with pd.ExcelWriter(output_path) as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=_safe_sheet_name(sheet_name), index=False)
    return output_path
