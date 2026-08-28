import argparse
import os
import sqlite3
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
import requests

from resource_config import (
    CLASSYFIRE_BRIDGE_R,
    HMDB_LITE_CSV,
    MW_DATABASE_CSV,
    MW_UNMAPPED_DATABASE_CSV,
    PUBCHEM_OFFLINE_SQLITE,
    REFMET_BRIDGE_R,
)
from knowledge_annotation_step import build_stage5_annotation_sheets
from pipeline_utils import (
    HAS_RPY2,
    _clean_cid,
    _clean_inchikey,
    _clean_scalar,
    _first_existing_column,
    _split_synonyms,
    _sqlite_table_exists,
    annotate_classyfire,
    extract_ids_from_pugview,
    fetch_cids_by_name_sqlite,
    fetch_pugview_ids_for_cids,
    load_classyfire_bridge,
    map_refmet_queries,
    timed_step,
    write_excel_workbook,
)


def _clean_kegg_ids(x: object) -> List[str]:
    s = _clean_scalar(x)
    if not s:
        return []
    s = s.replace('"', "")
    parts = re.split(r"[;,/|]", s)
    out, seen = [], set()
    for p in parts:
        p = p.strip()
        if re.fullmatch(r"C\d{5}", p) and p not in seen:
            seen.add(p)
            out.append(p)
    return out


def _normalize_name(x: object) -> Optional[str]:
    s = _clean_scalar(x)
    if not s:
        return None
    s = re.sub(r"[‐‑‒–—−]", "-", s)
    s = s.replace('"', "").replace("*", "").replace("&", "")
    s = re.sub(r"\s+", " ", s.strip().lower())
    s = re.sub(r"\s*\(\s*", "(", s)
    s = re.sub(r"\s*\)\s*", ")", s)
    s = re.sub(r"\s*-\s*", "-", s)
    return s


def _lower_name(x: object) -> Optional[str]:
    s = _clean_scalar(x)
    if not s:
        return None
    s = s.replace("*", "").replace("&", "")
    return re.sub(r"\s+", " ", s.lower()).strip() or None


def _refmet_bridge_query(x: object) -> Optional[str]:
    s = _clean_scalar(x)
    if not s:
        return None
    s = re.sub(r"[‐‑‒–—−]", "-", s)
    s = s.replace('"', "").replace("*", "").replace("&", "")
    s = re.sub(r"\s+", " ", s).strip()
    return s or None


def _aggressive_name_key(x: object) -> Optional[str]:
    s = _clean_scalar(x)
    if not s:
        return None
    s = re.sub(r"[‐‑‒–—−]", "-", s)
    s = s.replace('"', "").replace("*", "").replace("&", "")
    s = re.sub(r"\s+", " ", s).strip().lower()
    s = re.sub(r"[^a-z0-9]", "", s)
    return s or None


def _candidate_lookup_names(x: object) -> List[str]:
    s = _clean_scalar(x)
    if not s:
        return []

    s = re.sub(r"[‐‑‒–—−]", "-", s)
    s = s.replace('"', "").replace("*", "").replace("&", "")
    s = re.sub(r"\s+", " ", s).strip()

    variants: List[str] = [s]
    split_parts = [p.strip() for p in re.split(r"\s+(?:or|and)\s+", s, flags=re.I) if p.strip()]
    variants.extend(split_parts)

    expanded: List[str] = []
    for cand in variants:
        expanded.append(cand)

        swapped = re.sub(r";O([23])(?=($|/))", lambda m: f";{m.group(1)}O", cand, flags=re.I)
        if swapped != cand:
            expanded.append(swapped)

        if re.match(r"^\s*SM\b", cand, flags=re.I):
            paren_match = re.search(r"\(([^)]*)\)\s*$", cand)
            bracket_match = re.search(r"\[([^]]*)\]\s*$", cand)
            brace_match = re.search(r"\{([^}]*)\}\s*$", cand)
            tail_value = None
            if paren_match:
                tail_value = paren_match.group(1)
            elif bracket_match:
                tail_value = bracket_match.group(1)
            elif brace_match:
                tail_value = brace_match.group(1)

            # Do not collapse lipid shorthand like "SM (d14:0/18:2)" to bare "SM":
            # that generic key matches many unrelated PubChem synonyms and can
            # incorrectly assign all SM species to one arbitrary CID.
            is_specific_sm_lipid = bool(
                tail_value
                and (
                    "/" in tail_value
                    or ":" in tail_value
                    or re.search(r"\b[dtmeop]\d", tail_value, flags=re.I)
                )
            )

            if not is_specific_sm_lipid:
                trimmed = re.sub(r"\s*\([^)]*\)\s*$", "", cand).strip()
                trimmed = re.sub(r"\s*\[[^]]*\]\s*$", "", trimmed).strip()
                trimmed = re.sub(r"\s*\{[^}]*\}\s*$", "", trimmed).strip()
                if trimmed and trimmed != cand:
                    expanded.append(trimmed)
                    swapped_trimmed = re.sub(r";O([23])(?=($|/))", lambda m: f";{m.group(1)}O", trimmed, flags=re.I)
                    if swapped_trimmed != trimmed:
                        expanded.append(swapped_trimmed)

    out: List[str] = []
    seen = set()
    for cand in expanded:
        norm = re.sub(r"\s+", " ", cand).strip()
        if norm and norm not in seen:
            seen.add(norm)
            out.append(norm)
    return out


def _build_name_to_cid_map_aggressive(df: pd.DataFrame, name_col: str, cid_col: str = "PubChem_CID") -> Dict[str, str]:
    out: Dict[str, str] = {}
    if name_col not in df.columns or cid_col not in df.columns:
        return out
    for _, rr in df[[name_col, cid_col]].iterrows():
        key = _aggressive_name_key(rr[name_col])
        cid = _clean_cid(rr[cid_col])
        if key and cid and key not in out:
            out[key] = cid
    return out


def _build_synonym_to_cid_map_aggressive(df: pd.DataFrame, synonym_col: str = "synonyms", cid_col: str = "PubChem_CID") -> Dict[str, str]:
    out: Dict[str, str] = {}
    if synonym_col not in df.columns or cid_col not in df.columns:
        return out
    for _, rr in df[[cid_col, synonym_col]].iterrows():
        cid = _clean_cid(rr[cid_col])
        if not cid:
            continue
        for syn in _split_synonyms(rr[synonym_col]):
            key = _aggressive_name_key(syn)
            if key and key not in out:
                out[key] = cid
    return out


def _resolve_candidate_cid(
    name_value: object,
    lower_maps: Sequence[Dict[str, str]],
    aggressive_maps: Sequence[Dict[str, str]],
) -> Optional[str]:
    for cand in _candidate_lookup_names(name_value):
        key_lc = _lower_name(cand)
        if key_lc:
            for mp in lower_maps:
                if key_lc in mp:
                    return mp[key_lc]

        key_aggr = _aggressive_name_key(cand)
        if key_aggr:
            for mp in aggressive_maps:
                if key_aggr in mp:
                    return mp[key_aggr]
    return None


def _parallel_map_ordered(
    values: Sequence[object],
    func,
    max_workers: int,
) -> List[object]:
    if not values:
        return []
    if max_workers <= 1 or len(values) == 1:
        return [func(v) for v in values]

    out: List[object] = [None] * len(values)
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(func, value): idx for idx, value in enumerate(values)}
        for fut in as_completed(futures):
            idx = futures[fut]
            out[idx] = fut.result()
    return out


def _relaxed_name_keys(x: object) -> List[str]:
    n = _normalize_name(x)
    if not n:
        return []
    keys = [n]

    no_pos = re.sub(r"^\s*[12]-", "", n)
    if no_pos != n:
        keys.append(no_pos)

    compact = re.sub(r"[-\s]+", "", n)
    if compact and compact != n:
        keys.append(compact)

    compact_no_pos = re.sub(r"[-\s]+", "", no_pos)
    if compact_no_pos and compact_no_pos not in keys:
        keys.append(compact_no_pos)

    out, seen = [], set()
    for k in keys:
        if k and k not in seen:
            seen.add(k)
            out.append(k)
    return out


def _smart_split_metabolite(x: object) -> List[str]:
    s = _clean_scalar(x)
    if not s:
        return []

    parts, buf = [], []
    depth_sq = 0
    depth_par = 0
    for ch in s:
        if ch == "[":
            depth_sq += 1
        elif ch == "]":
            depth_sq = max(0, depth_sq - 1)
        elif ch == "(":
            depth_par += 1
        elif ch == ")":
            depth_par = max(0, depth_par - 1)

        if ch == "/" and depth_sq == 0 and depth_par == 0:
            part = "".join(buf).strip()
            if part:
                parts.append(part)
            buf = []
        else:
            buf.append(ch)

    tail = "".join(buf).strip()
    if tail:
        parts.append(tail)

    out, seen = [], set()
    for p in parts:
        n = _normalize_name(p)
        if n and n not in seen:
            seen.add(n)
            out.append(n)
    return out


def fetch_cids_by_name_sqlite(db_path: str, name: str, limit: int = 50) -> List[str]:
    if not name:
        return []

    candidates = _candidate_lookup_names(name)
    if not candidates:
        return []

    aggr_keys = []
    seen_keys: Set[str] = set()
    for candidate in candidates:
        key = _aggressive_name_key(candidate)
        if key and key not in seen_keys:
            seen_keys.add(key)
            aggr_keys.append(key)

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows = []

    if aggr_keys and _sqlite_table_exists(db_path, "cid_syn_aggr"):
        for i in range(0, len(aggr_keys), 900):
            chunk = aggr_keys[i:i + 900]
            placeholders = ",".join(["?"] * len(chunk))
            rows.extend(
                cur.execute(
                    f"SELECT DISTINCT cid FROM cid_syn_aggr WHERE syn_aggr IN ({placeholders}) LIMIT ?",
                    (*chunk, limit),
                ).fetchall()
            )
    elif _sqlite_table_exists(db_path, "cid_syn"):
        norm_names = []
        seen_norm: Set[str] = set()
        for candidate in candidates:
            qname = _normalize_name(candidate)
            if qname and qname not in seen_norm:
                seen_norm.add(qname)
                norm_names.append(qname)
        for i in range(0, len(norm_names), 900):
            chunk = norm_names[i:i + 900]
            placeholders = ",".join(["?"] * len(chunk))
            rows.extend(
                cur.execute(
                    f"SELECT DISTINCT cid FROM cid_syn WHERE lower(syn) IN ({placeholders}) LIMIT ?",
                    (*chunk, limit),
                ).fetchall()
            )

    con.close()
    out = []
    seen_cids: Set[str] = set()
    for row in rows:
        cid = _clean_cid(row[0])
        if cid and cid not in seen_cids:
            seen_cids.add(cid)
            out.append(cid)
            if len(out) >= limit:
                break
    return out


def _naify(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().replace(
        {"": pd.NA, "-": pd.NA, "NA": pd.NA, "Na": pd.NA, "null": pd.NA, "NA_character_": pd.NA, "<NA>": pd.NA}
    )


def _ensure_cols(df: pd.DataFrame, cols: Sequence[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        if c not in out.columns:
            out[c] = pd.NA
        out[c] = _naify(out[c])
    return out


def _has_primary_mapping_identifier(df: pd.DataFrame) -> pd.Series:
    id_candidates = ["PubChem_CID", "KEGG_ID", "HMDB_ID"]
    present = [c for c in id_candidates if c in df.columns]
    if not present:
        return pd.Series([False] * len(df), index=df.index)

    bits = []
    for c in present:
        bits.append(_naify(df[c]).notna())
    out = bits[0].copy()
    for b in bits[1:]:
        out = out | b
    return out


def _first_valid_cid(series: pd.Series) -> Optional[str]:
    for v in series.tolist():
        c = _clean_cid(v)
        if c:
            return c
    return None


def _first_valid_scalar(series: pd.Series) -> Optional[str]:
    for v in series.tolist():
        c = _clean_scalar(v)
        if c:
            return c
    return None


def _build_mw_name_maps(
    mw_df: pd.DataFrame,
    name_col: str,
    pubchem_col: Optional[str],
    kegg_col: Optional[str],
) -> Tuple[Dict[str, str], Dict[str, str]]:
    tmp = mw_df.copy()
    tmp["__norm_name"] = tmp[name_col].map(_normalize_name)
    tmp = tmp.dropna(subset=["__norm_name"])

    pubchem_map: Dict[str, str] = {}
    kegg_map: Dict[str, str] = {}

    grouped = tmp.groupby("__norm_name", dropna=False, sort=False)
    for nm, grp in grouped:
        if pubchem_col and pubchem_col in grp.columns:
            pc = _first_valid_cid(grp[pubchem_col])
            if pc:
                pubchem_map[str(nm)] = pc
        if kegg_col and kegg_col in grp.columns:
            kg = _first_valid_scalar(grp[kegg_col])
            if kg:
                kegg_map[str(nm)] = kg
    return pubchem_map, kegg_map


def _build_name_to_cid_map_lower(df: pd.DataFrame, name_col: str, cid_col: str) -> Dict[str, str]:
    work = df[[name_col, cid_col]].copy()
    work["_name_lc"] = work[name_col].map(_lower_name)
    work["_cid"] = work[cid_col].map(_clean_cid)
    work = work.dropna(subset=["_name_lc", "_cid"])
    if work.empty:
        return {}
    work = work.drop_duplicates(subset=["_name_lc"], keep="first")
    return dict(zip(work["_name_lc"], work["_cid"]))


def _build_synonym_to_cid_map_lower(df: pd.DataFrame, synonym_col: str = "synonyms", cid_col: str = "PubChem_CID") -> Dict[str, str]:
    out: Dict[str, str] = {}
    if synonym_col not in df.columns or cid_col not in df.columns:
        return out
    for _, rr in df[[cid_col, synonym_col]].iterrows():
        cid = _clean_cid(rr[cid_col])
        if not cid:
            continue
        for syn in _split_synonyms(rr[synonym_col]):
            key = _lower_name(syn)
            if key and key not in out:
                out[key] = cid
    return out


def _coalesce_columns(df: pd.DataFrame, candidates: Sequence[str]) -> pd.Series:
    out = pd.Series([pd.NA] * len(df), dtype="string")
    for c in candidates:
        if c in df.columns:
            s = df[c].astype("string")
            out = out.where(~out.isna(), s)
    return out


def _extract_pubchem_from_kegg(kegg_id: str, timeout: int = 20) -> Optional[str]:
    kid = _clean_scalar(kegg_id)
    if not kid:
        return None
    try:
        res = requests.get(f"https://rest.kegg.jp/get/{kid}", timeout=timeout)
        if res.status_code != 200:
            return None
        txt = res.text
    except Exception:
        return None

    for line in txt.splitlines():
        if "PubChem" not in line:
            continue
        m = re.search(r"PubChem[^0-9]*(\d+)", line)
        if m:
            return m.group(1)
    return None


def _first_mapping_from_split_key(df: pd.DataFrame, key_col: str, val_col: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for _, row in df[[key_col, val_col]].dropna().iterrows():
        val = _clean_cid(row[val_col])
        if not val:
            continue
        for kid in _clean_kegg_ids(row[key_col]):
            if kid not in mapping:
                mapping[kid] = val
    return mapping


def _fetch_sqlite_map(db_path: str, table: str, cid_col: str, val_col: str, cids: Set[str]) -> Dict[str, str]:
    if not cids or not _sqlite_table_exists(db_path, table):
        return {}
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    out: Dict[str, str] = {}
    cid_list = list(cids)
    for i in range(0, len(cid_list), 900):
        chunk = cid_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = f"SELECT {cid_col}, {val_col} FROM {table} WHERE {cid_col} IN ({placeholders})"
        for cid, val in cur.execute(q, chunk).fetchall():
            cc = _clean_cid(cid)
            vv = _clean_scalar(val)
            if cc and vv:
                out[cc] = vv
    con.close()
    return out


def _fetch_synonyms_sqlite(db_path: str, cids: Set[str], max_synonyms_per_cid: int = 200) -> Dict[str, str]:
    if not cids or not _sqlite_table_exists(db_path, "cid_syn"):
        return {}
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    agg = {cid: [] for cid in cids}
    seen = {cid: set() for cid in cids}
    cid_list = list(cids)
    for i in range(0, len(cid_list), 900):
        chunk = cid_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        for cid, syn in cur.execute(f"SELECT cid, syn FROM cid_syn WHERE cid IN ({placeholders})", chunk).fetchall():
            cid = _clean_cid(cid)
            syn = _clean_scalar(syn)
            if not cid or not syn:
                continue
            if len(agg[cid]) >= max_synonyms_per_cid or syn in seen[cid]:
                continue
            seen[cid].add(syn)
            agg[cid].append(syn)
    con.close()
    return {k: "; ".join(v) for k, v in agg.items() if v}


def _fetch_synonym_name_to_cid_sqlite(db_path: str, names_lc: Set[str]) -> Dict[str, str]:
    if not names_lc:
        return {}
    if _sqlite_table_exists(db_path, "cid_syn_aggr"):
        out: Dict[str, str] = {}
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        key_to_names: Dict[str, List[str]] = {}
        for name in names_lc:
            key = _aggressive_name_key(name)
            if key:
                key_to_names.setdefault(key, []).append(name)
        key_list = sorted(key_to_names)
        for i in range(0, len(key_list), 900):
            chunk = key_list[i:i + 900]
            placeholders = ",".join(["?"] * len(chunk))
            q = f"SELECT cid, syn_aggr FROM cid_syn_aggr WHERE syn_aggr IN ({placeholders})"
            for cid, syn_aggr in cur.execute(q, chunk).fetchall():
                c = _clean_cid(cid)
                for original_name in key_to_names.get(_clean_scalar(syn_aggr) or "", []):
                    if c and original_name and original_name not in out:
                        out[original_name] = c
        con.close()
        return out

    if not _sqlite_table_exists(db_path, "cid_syn"):
        return {}

    out: Dict[str, str] = {}
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    name_list = sorted([n for n in names_lc if n])
    for i in range(0, len(name_list), 900):
        chunk = name_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = f"SELECT cid, syn FROM cid_syn WHERE lower(syn) IN ({placeholders})"
        for cid, syn in cur.execute(q, chunk).fetchall():
            c = _clean_cid(cid)
            k = _lower_name(syn)
            if c and k and k not in out:
                out[k] = c
    con.close()
    return out


def _fetch_synonym_candidates_sqlite(db_path: str, value: object, limit: int = 50) -> List[str]:
    return fetch_cids_by_name_sqlite(db_path, _clean_scalar(value), limit=limit)


def _fetch_synonym_name_to_cid_sqlite_aggressive(db_path: str, names_aggr: Set[str]) -> Dict[str, str]:
    if not names_aggr:
        return {}
    if not _sqlite_table_exists(db_path, "cid_syn_aggr"):
        return {}

    out: Dict[str, str] = {}
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    key_list = sorted([n for n in names_aggr if n])
    for i in range(0, len(key_list), 900):
        chunk = key_list[i:i + 900]
        placeholders = ",".join(["?"] * len(chunk))
        q = f"SELECT cid, syn_aggr FROM cid_syn_aggr WHERE syn_aggr IN ({placeholders})"
        for cid, syn_aggr in cur.execute(q, chunk).fetchall():
            c = _clean_cid(cid)
            k = _clean_scalar(syn_aggr)
            if c and k and k not in out:
                out[k] = c
    con.close()
    return out


def _load_pubchem_names_from_sqlite(
    sqlite_db: str,
    needed_cids: Optional[Set[str]] = None,
    needed_names_norm: Optional[Set[str]] = None,
) -> pd.DataFrame:
    if not _sqlite_table_exists(sqlite_db, "cid_title") and not _sqlite_table_exists(sqlite_db, "cid_iupac"):
        raise ValueError("SQLite DB is missing both cid_title and cid_iupac tables")

    con = sqlite3.connect(sqlite_db)
    title_exists = _sqlite_table_exists(sqlite_db, "cid_title")
    iupac_exists = _sqlite_table_exists(sqlite_db, "cid_iupac")

    needed_cids = {c for c in (needed_cids or set()) if _clean_cid(c)}
    needed_names_norm = {n for n in (needed_names_norm or set()) if n}

    rows_by_cid: Dict[str, Dict[str, object]] = {}

    def _upsert(cid: object, title: object = pd.NA, iupac: object = pd.NA):
        cc = _clean_cid(cid)
        if not cc:
            return
        if cc not in rows_by_cid:
            rows_by_cid[cc] = {"PubChem_CID": cc, "title": pd.NA, "iupac": pd.NA, "synonyms": pd.NA}
        if _clean_scalar(title):
            rows_by_cid[cc]["title"] = title
        if _clean_scalar(iupac):
            rows_by_cid[cc]["iupac"] = iupac

    def _query_in_chunks(sql: str, vals: List[str], chunk_size: int = 900):
        if not vals:
            return
        cur = con.cursor()
        for i in range(0, len(vals), chunk_size):
            chunk = vals[i:i + chunk_size]
            placeholders = ",".join(["?"] * len(chunk))
            q = sql.format(placeholders=placeholders)
            for rec in cur.execute(q, chunk).fetchall():
                yield rec

    if needed_cids:
        cid_list = sorted(needed_cids)
        if title_exists:
            for cid, title in _query_in_chunks(
                "SELECT cid, title FROM cid_title WHERE cid IN ({placeholders})",
                cid_list,
            ):
                _upsert(cid, title=title)
        if iupac_exists:
            for cid, iupac in _query_in_chunks(
                "SELECT cid, iupac_name FROM cid_iupac WHERE cid IN ({placeholders})",
                cid_list,
            ):
                _upsert(cid, iupac=iupac)

    if needed_names_norm:
        name_list = sorted(needed_names_norm)
        if title_exists:
            for cid, title in _query_in_chunks(
                "SELECT cid, title FROM cid_title WHERE lower(title) IN ({placeholders})",
                name_list,
            ):
                _upsert(cid, title=title)
        if iupac_exists:
            for cid, iupac in _query_in_chunks(
                "SELECT cid, iupac_name FROM cid_iupac WHERE lower(iupac_name) IN ({placeholders})",
                name_list,
            ):
                _upsert(cid, iupac=iupac)

    out = pd.DataFrame(rows_by_cid.values())
    con.close()

    if out.empty:
        out = pd.DataFrame(columns=["PubChem_CID", "title", "iupac", "synonyms"])
    out = out.astype("string")
    for c in ["PubChem_CID", "title", "iupac", "synonyms"]:
        if c not in out.columns:
            out[c] = pd.NA
    return out[["PubChem_CID", "title", "iupac", "synonyms"]]


def load_pubchem_names_index(
    sqlite_db: str,
    needed_cids: Optional[Set[str]] = None,
    needed_names_norm: Optional[Set[str]] = None,
) -> pd.DataFrame:
    return _load_pubchem_names_from_sqlite(
        sqlite_db,
        needed_cids=needed_cids,
        needed_names_norm=needed_names_norm,
    )


def _fetch_pubchem_props_chunk(cids: List[str], timeout: int = 30) -> pd.DataFrame:
    if not cids:
        return pd.DataFrame(
            columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass", "InChIKey", "SMILES"],
            dtype="string",
        )
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{','.join(cids)}/property/Title,MolecularFormula,ExactMass,InChIKey,CanonicalSMILES/JSON"
    )
    try:
        res = requests.get(url, timeout=timeout)
        if res.status_code != 200:
            return pd.DataFrame(
                columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass", "InChIKey", "SMILES"],
                dtype="string",
            )
        props = res.json().get("PropertyTable", {}).get("Properties", [])
    except Exception:
        return pd.DataFrame(
            columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass", "InChIKey", "SMILES"],
            dtype="string",
        )

    rows = []
    for r in props:
        cid = _clean_cid(r.get("CID"))
        if not cid:
            continue
        rows.append(
            {
                "PubChem_CID": cid,
                "Standardized_Name": _clean_scalar(r.get("Title")),
                "Molecular_Formula": _clean_scalar(r.get("MolecularFormula")),
                "Exact_Mass": str(r.get("ExactMass")) if r.get("ExactMass") is not None else pd.NA,
                "InChIKey": _clean_scalar(r.get("InChIKey")),
                "SMILES": _clean_scalar(r.get("CanonicalSMILES")),
            }
        )
    if not rows:
        return pd.DataFrame(
            columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass", "InChIKey", "SMILES"],
            dtype="string",
        )
    return pd.DataFrame(rows).astype("string").drop_duplicates("PubChem_CID")


def fetch_pubchem_props_api(cids: Sequence[str], workers: int = 8, chunk_size: int = 100) -> pd.DataFrame:
    cleaned = sorted({c for c in (_clean_cid(x) for x in cids) if c})
    if not cleaned:
        return pd.DataFrame(
            columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass", "InChIKey", "SMILES"],
            dtype="string",
        )
    chunks = [cleaned[i:i + chunk_size] for i in range(0, len(cleaned), chunk_size)]
    frames: List[pd.DataFrame] = []
    if workers <= 1 or len(chunks) == 1:
        frames = [_fetch_pubchem_props_chunk(c) for c in chunks]
    else:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futures = [ex.submit(_fetch_pubchem_props_chunk, c) for c in chunks]
            for fut in as_completed(futures):
                try:
                    frames.append(fut.result())
                except Exception:
                    continue
    if not frames:
        return pd.DataFrame(
            columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass", "InChIKey", "SMILES"],
            dtype="string",
        )
    out = pd.concat(frames, ignore_index=True)
    if out.empty:
        return pd.DataFrame(
            columns=["PubChem_CID", "Standardized_Name", "Molecular_Formula", "Exact_Mass", "InChIKey", "SMILES"],
            dtype="string",
        )
    return out.astype("string").drop_duplicates("PubChem_CID")


def _normalize_chebi(val: object) -> Optional[str]:
    s = _clean_scalar(val)
    if not s:
        return None
    m = re.search(r"(?:CHEBI:)?(\d+)", s, flags=re.IGNORECASE)
    return m.group(1) if m else None


PUBLIC_STAGE_COLUMNS = [
    "HuMANet_ID", "original_query_name", "Query_Name", "Database_Source", "InChIKey", "HMDB_ID", "PubChem_CID", "KEGG_ID", "ChEBI_ID",
    "Molecular_Formula", "Exact_Mass", "Super_Class", "Main_Class", "Sub_Class", "Standardized_Name", "SMILES",
    "Synonyms", "matched_name", "Study_Folder",
]


def _finalize_public_stage_output(df: pd.DataFrame, columns=PUBLIC_STAGE_COLUMNS) -> pd.DataFrame:
    out = df.copy()
    if "Study_Folder" not in out.columns:
        if "Input_File" in out.columns:
            out["Study_Folder"] = out["Input_File"].astype("string")
        else:
            out["Study_Folder"] = pd.Series([pd.NA] * len(out), dtype="string")
    else:
        out["Study_Folder"] = out["Study_Folder"].astype("string")
    for col in columns:
        if col not in out.columns:
            out[col] = pd.NA
    return out[columns].copy()


def run_stage3(
    stage2_unmapped_csv: str,
    mw_database_csv: str,
    out_prefix: str,
    hmdb_lite_csv: Optional[str] = None,
    sqlite_db: Optional[str] = None,
    pubchem_workers: int = 8,
    r_classyfire_bridge: Optional[str] = None,
    skip_classyfire: bool = False,
    print_classyfire_output: bool = False,
    sleep: float = 0.2,
    sqlite_name_fallback: bool = True,
    mw_unmapped_db_csv: Optional[str] = None,
    skip_stage5_excel: bool = False,
) -> Dict[str, pd.DataFrame]:
    total_start = time.perf_counter()

    with timed_step("Load Stage3 inputs"):
        s2 = pd.read_csv(stage2_unmapped_csv, dtype="string")
        if "original_query_name" not in s2.columns:
            s2["original_query_name"] = s2["Query_Name"].astype("string")
        mw = pd.read_csv(mw_database_csv, dtype="string")
        mw.columns = [str(c).strip() for c in mw.columns]
        hmdb = pd.read_csv(hmdb_lite_csv, dtype="string") if hmdb_lite_csv else None
        mw_unmapped = None
        if mw_unmapped_db_csv and os.path.exists(mw_unmapped_db_csv):
            mw_unmapped = pd.read_csv(mw_unmapped_db_csv, dtype="string")
            mw_unmapped.columns = [str(c).strip() for c in mw_unmapped.columns]
        else:
            mw_unmapped = mw.copy()
            if mw_unmapped_db_csv:
                print(f"[WARN] MW unmapped DB not found, using mw_database for combined recovery: {mw_unmapped_db_csv}")

    with timed_step("Normalize required columns"):
        s2 = _ensure_cols(
            s2,
            [
                "Query_Name",
                "original_query_name",
                "PubChem_CID",
                "KEGG_ID",
                "HMDB_ID",
                "ChEBI_ID",
                "Standardized_Name",
                "Molecular_Formula",
                "Exact_Mass",
                "InChIKey",
                "SMILES",
                "Synonyms",
                "matched_name",
                "Super_Class",
                "Main_Class",
                "Sub_Class",
            ],
        )

    with timed_step("MW name-based CID/KEGG mapping"):
        mw_name_col = _first_existing_column(mw.columns, ["Metabolite.Name", "Metabolite_Name", "Metabolite Name", "compound_name"])
        mw_pubchem_col = _first_existing_column(mw.columns, ["PubChem.Compound_ID", "PubChem_CID", "PubChem CID", "PubChem.Compound.Id"])
        mw_kegg_col = _first_existing_column(mw.columns, ["Kegg.Id", "KEGG_ID", "Kegg ID", "KEGG Id"])

        if mw_name_col:
            mw_local = mw.copy()
            mw_local["__norm_name"] = mw_local[mw_name_col].map(_normalize_name)
            pubchem_map, kegg_map = _build_mw_name_maps(mw_local, mw_name_col, mw_pubchem_col, mw_kegg_col)

            q_norm = s2["Query_Name"].map(_normalize_name)
            s2["PubChem_CID"] = _naify(q_norm.map(pubchem_map)).combine_first(_naify(s2["PubChem_CID"]))
            s2["KEGG_ID"] = _naify(q_norm.map(kegg_map)).combine_first(_naify(s2["KEGG_ID"]))

            missing_after_exact = s2["PubChem_CID"].isna() & s2["KEGG_ID"].isna()
            if missing_after_exact.any():
                relaxed_pubchem: Dict[str, str] = {}
                relaxed_kegg: Dict[str, str] = {}
                for _, r in mw_local.dropna(subset=["__norm_name"]).iterrows():
                    for key in _relaxed_name_keys(r["__norm_name"]):
                        if key not in relaxed_pubchem and mw_pubchem_col:
                            v = _clean_cid(r.get(mw_pubchem_col))
                            if v:
                                relaxed_pubchem[key] = v
                        if key not in relaxed_kegg and mw_kegg_col:
                            vv = _clean_scalar(r.get(mw_kegg_col))
                            if vv:
                                relaxed_kegg[key] = vv

                def _resolve_relaxed_hits(q: object) -> Tuple[object, object]:
                    hit_pc = pd.NA
                    hit_kg = pd.NA
                    for key in _relaxed_name_keys(q):
                        if key in relaxed_pubchem:
                            hit_pc = relaxed_pubchem[key]
                            break
                    for key in _relaxed_name_keys(q):
                        if key in relaxed_kegg:
                            hit_kg = relaxed_kegg[key]
                            break
                    return hit_pc, hit_kg

                relaxed_hits = _parallel_map_ordered(
                    s2.loc[missing_after_exact, "Query_Name"].tolist(),
                    _resolve_relaxed_hits,
                    max_workers=pubchem_workers,
                )
                fill_pc = [x[0] for x in relaxed_hits]
                fill_kg = [x[1] for x in relaxed_hits]

                s2.loc[missing_after_exact, "PubChem_CID"] = pd.Series(fill_pc, index=s2.index[missing_after_exact], dtype="string")
                s2.loc[missing_after_exact, "KEGG_ID"] = _naify(pd.Series(fill_kg, index=s2.index[missing_after_exact], dtype="string"))

        s2["PubChem_CID"] = s2["PubChem_CID"].map(_clean_cid).astype("string")
        s2["KEGG_ID"] = _naify(s2["KEGG_ID"])

    with timed_step("KEGG->PubChem fallback via MW DB"):
        mw_pubchem_col = _first_existing_column(mw.columns, ["PubChem.Compound_ID", "PubChem_CID", "PubChem CID", "PubChem.Compound.Id"])
        mw_kegg_col = _first_existing_column(mw.columns, ["Kegg.Id", "KEGG_ID", "Kegg ID", "KEGG Id"])
        if mw_kegg_col and mw_pubchem_col:
            kegg_to_cid = _first_mapping_from_split_key(mw, mw_kegg_col, mw_pubchem_col)
            missing = s2["PubChem_CID"].isna()
            fill_vals = []
            for k in s2.loc[missing, "KEGG_ID"].tolist():
                hit = next((kid for kid in _clean_kegg_ids(k) if kid in kegg_to_cid), None)
                fill_vals.append(kegg_to_cid.get(hit, pd.NA) if hit else pd.NA)
            s2.loc[missing, "PubChem_CID"] = pd.Series(fill_vals, index=s2.index[missing], dtype="string")

    with timed_step("KEGG->PubChem fallback via HMDB"):
        if hmdb is not None and "hmdb_metabolite_kegg_id" in hmdb.columns and "hmdb_metabolite_pubchem_compound_id" in hmdb.columns:
            hmdb_map = _first_mapping_from_split_key(
                hmdb.rename(
                    columns={
                        "hmdb_metabolite_kegg_id": "Kegg.Id",
                        "hmdb_metabolite_pubchem_compound_id": "PubChem.Compound_ID",
                    }
                ),
                "Kegg.Id",
                "PubChem.Compound_ID",
            )
            missing = s2["PubChem_CID"].isna()
            fill_vals = []
            for k in s2.loc[missing, "KEGG_ID"].tolist():
                hit = next((kid for kid in _clean_kegg_ids(k) if kid in hmdb_map), None)
                fill_vals.append(hmdb_map.get(hit, pd.NA) if hit else pd.NA)
            s2.loc[missing, "PubChem_CID"] = pd.Series(fill_vals, index=s2.index[missing], dtype="string")

        s2["PubChem_CID"] = s2["PubChem_CID"].map(_clean_cid).astype("string")

    with timed_step("KEGG->PubChem online fallback (KEGG REST)"):
        missing = s2["PubChem_CID"].isna() & s2["KEGG_ID"].notna()
        if missing.any():
            def _resolve_kegg_pubchem(k: object) -> object:
                found = None
                for kid in _clean_kegg_ids(k):
                    found = _extract_pubchem_from_kegg(kid)
                    if found:
                        break
                return found if found else pd.NA

            fill_vals = _parallel_map_ordered(
                s2.loc[missing, "KEGG_ID"].tolist(),
                _resolve_kegg_pubchem,
                max_workers=pubchem_workers,
            )
            s2.loc[missing, "PubChem_CID"] = pd.Series(fill_vals, index=s2.index[missing], dtype="string")

        s2["PubChem_CID"] = s2["PubChem_CID"].map(_clean_cid).astype("string")

    with timed_step("SQLite synonym name->CID fallback"):
        if sqlite_name_fallback and sqlite_db and (_sqlite_table_exists(sqlite_db, "cid_syn_aggr") or _sqlite_table_exists(sqlite_db, "cid_syn")):
            missing = s2["PubChem_CID"].isna()
            def _resolve_sqlite_synonym(q: object) -> object:
                cids_found: List[str] = []
                for qn in _smart_split_metabolite(q):
                    cids_found.extend(fetch_cids_by_name_sqlite(sqlite_db, qn, limit=20))
                return next((c for c in cids_found if _clean_cid(c)), pd.NA)

            fill_vals = _parallel_map_ordered(
                s2.loc[missing, "Query_Name"].tolist(),
                _resolve_sqlite_synonym,
                max_workers=pubchem_workers,
            )
            s2.loc[missing, "PubChem_CID"] = pd.Series(fill_vals, index=s2.index[missing], dtype="string")
            s2["PubChem_CID"] = s2["PubChem_CID"].map(_clean_cid).astype("string")

    with timed_step("Final unresolved-metabolite CID recovery without intermediate outputs"):
        if mw_unmapped is not None:
            for c in ["matched_name"]:
                if c not in s2.columns:
                    s2[c] = pd.NA

            s2["Query_Name_toCheck"] = s2["Query_Name"].astype("string")
            s2["Query_Name_toCheck"] = s2["Query_Name_toCheck"].str.replace("*", "", regex=False)
            s2["Query_Name_toCheck"] = s2["Query_Name_toCheck"].str.replace("&", "", regex=False)
            s2["Query_Name_toCheck"] = s2["Query_Name_toCheck"].str.replace(r"\.[0-9]+$", "", regex=True)
            s2["Query_Name_toCheck"] = s2["Query_Name_toCheck"].str.replace(r"\s*\((\d+)\)$|\s*\[(\d+)\]$", "", regex=True)
            s2["Query_Name_toCheck"] = s2["Query_Name_toCheck"].map(_lower_name).astype("string")

            mw_name_col = _first_existing_column(
                mw_unmapped.columns.tolist(),
                ["Metabolite.Name", "Metabolite_Name", "Metabolite Name", "Query_Name", "name", "title"],
            )
            mw_cid_col = _first_existing_column(
                mw_unmapped.columns.tolist(),
                ["PubChem_CID", "PubChem.Compound_ID", "CID", "cid"],
            )
            mw_refmet_std_col = _first_existing_column(
                mw_unmapped.columns.tolist(),
                [
                    "RefMet.Name.Standardized.Name.",
                    "RefMet_Name_Standardized_Name",
                    "RefMet.Standardized.Name",
                    "RefMet.Name.Standardized.Name",
                    "RefMet_Standardized_Name",
                ],
            )

            if mw_name_col and mw_refmet_std_col:
                refmet_work = mw_unmapped[[mw_name_col, mw_refmet_std_col]].copy()
                refmet_work["__query_key"] = refmet_work[mw_name_col].map(_normalize_name)
                refmet_work = refmet_work.dropna(subset=["__query_key"])
                refmet_work = refmet_work.drop_duplicates(subset=["__query_key"], keep="first")
                refmet_std_map = dict(zip(refmet_work["__query_key"], refmet_work[mw_refmet_std_col].astype("string")))
                s2["refmet_standardized_name"] = s2["Query_Name"].map(_normalize_name).map(refmet_std_map)
            else:
                s2["refmet_standardized_name"] = pd.NA

            s2["refmet_standardized_name"] = s2["refmet_standardized_name"].astype("string").str.replace("*", "", regex=False)
            s2["refmet_standardized_name"] = s2["refmet_standardized_name"].astype("string").str.replace("&", "", regex=False)
            s2["refmet_standardized_name"] = _naify(s2["refmet_standardized_name"])
            s2["refmet_standardized_name"] = s2["refmet_standardized_name"].map(_lower_name).astype("string")

            refmet_bridge_path = REFMET_BRIDGE_R if HAS_RPY2 else None
            refmet_missing = s2["refmet_standardized_name"].isna()
            if refmet_missing.any() and refmet_bridge_path:
                try:
                    bridge_queries = (
                        s2.loc[refmet_missing, "Query_Name"]
                        .map(_refmet_bridge_query)
                        .astype("string")
                        .tolist()
                    )
                    bridge_queries = [q for q in bridge_queries if _clean_scalar(q)]
                    bridge_df = map_refmet_queries(bridge_queries, refmet_bridge_path)
                    if not bridge_df.empty:
                        bridge_cols = [str(c) for c in bridge_df.columns]
                        input_col = _first_existing_column(bridge_cols, ["Input.name", "Input_name", "input_name", "query_name"])
                        std_col = _first_existing_column(bridge_cols, ["Standardized.name", "Standardized_name", "standardized_name", "RefMet.Name.Standardized.Name."])
                        if input_col and std_col:
                            bridge_work = bridge_df[[input_col, std_col]].copy()
                            bridge_work["__query_key"] = bridge_work[input_col].map(_refmet_bridge_query)
                            bridge_work[std_col] = bridge_work[std_col].astype("string").str.replace("*", "", regex=False)
                            bridge_work[std_col] = bridge_work[std_col].astype("string").str.replace("&", "", regex=False)
                            bridge_work[std_col] = bridge_work[std_col].map(_lower_name)
                            bridge_work = bridge_work.dropna(subset=["__query_key", std_col])
                            bridge_work = bridge_work.drop_duplicates(subset=["__query_key"], keep="first")
                            bridge_map = dict(zip(bridge_work["__query_key"], bridge_work[std_col].astype("string")))
                            bridge_fill = s2.loc[refmet_missing, "Query_Name"].map(_refmet_bridge_query).map(bridge_map).astype("string")
                            s2.loc[refmet_missing, "refmet_standardized_name"] = s2.loc[refmet_missing, "refmet_standardized_name"].astype("string").fillna(bridge_fill)
                except Exception as e:
                    print(f"[WARN] RefMet bridge fallback failed, continuing without it: {e}")

            refmet_debug = s2.loc[
                s2["refmet_standardized_name"].notna(),
                ["Query_Name", "refmet_standardized_name"],
            ].drop_duplicates().head(50)
            print(f"[DEBUG] RefMet standardized-name rows: {len(refmet_debug)} shown (up to 50)")
            if not refmet_debug.empty:
                print(refmet_debug.to_string(index=False))

            if sqlite_name_fallback and sqlite_db and (_sqlite_table_exists(sqlite_db, "cid_syn_aggr") or _sqlite_table_exists(sqlite_db, "cid_syn")):
                missing = s2["PubChem_CID"].isna() & s2["refmet_standardized_name"].notna()
                def _resolve_sqlite_refmet_synonym(q: object) -> object:
                    cids_found: List[str] = []
                    for candidate in _candidate_lookup_names(q):
                        cids_found.extend(fetch_cids_by_name_sqlite(sqlite_db, candidate, limit=20))
                    return next((c for c in cids_found if _clean_cid(c)), pd.NA)

                refmet_fill = _parallel_map_ordered(
                    s2.loc[missing, "refmet_standardized_name"].tolist(),
                    _resolve_sqlite_refmet_synonym,
                    max_workers=pubchem_workers,
                )
                s2.loc[missing, "PubChem_CID"] = pd.Series(refmet_fill, index=s2.index[missing], dtype="string")
                s2["PubChem_CID"] = s2["PubChem_CID"].map(_clean_cid).astype("string")

            if mw_name_col and mw_cid_col:
                mw_map_lc = _build_name_to_cid_map_lower(mw_unmapped, mw_name_col, mw_cid_col)
                mw_map_aggr = _build_name_to_cid_map_aggressive(mw_unmapped, mw_name_col, mw_cid_col)

                missing = s2["PubChem_CID"].isna()
                cid_fill = pd.Series(
                    _parallel_map_ordered(
                        s2.loc[missing, "Query_Name"].tolist(),
                        lambda x: _resolve_candidate_cid(x, [mw_map_lc], [mw_map_aggr]),
                        max_workers=pubchem_workers,
                    ),
                    index=s2.index[missing],
                    dtype="string",
                )
                s2.loc[missing, "PubChem_CID"] = cid_fill.astype("string")

                missing = s2["PubChem_CID"].isna()
                rcid_fill = pd.Series(
                    _parallel_map_ordered(
                        s2.loc[missing, "refmet_standardized_name"].tolist(),
                        lambda x: _resolve_candidate_cid(x, [mw_map_lc], [mw_map_aggr]),
                        max_workers=pubchem_workers,
                    ),
                    index=s2.index[missing],
                    dtype="string",
                )
                s2.loc[missing, "PubChem_CID"] = rcid_fill.astype("string")
                s2["PubChem_CID"] = s2["PubChem_CID"].map(_clean_cid).astype("string")

            if sqlite_db:
                try:
                    needed_cids = {c for c in s2["PubChem_CID"].astype("string").dropna().tolist() if _clean_cid(c)}
                    needed_names_norm = {
                        _lower_name(candidate)
                        for value in s2.loc[s2["PubChem_CID"].isna(), "Query_Name_toCheck"].tolist()
                        for candidate in _candidate_lookup_names(value)
                        if _lower_name(candidate)
                    } | {
                        _lower_name(candidate)
                        for value in s2.loc[s2["PubChem_CID"].isna(), "refmet_standardized_name"].tolist()
                        for candidate in _candidate_lookup_names(value)
                        if _lower_name(candidate)
                    }

                    names_idx = load_pubchem_names_index(
                        sqlite_db=sqlite_db,
                        needed_cids=needed_cids,
                        needed_names_norm=needed_names_norm,
                    )
                    names_idx["PubChem_CID"] = names_idx["PubChem_CID"].map(_clean_cid).astype("string")

                    title_map_lc = _build_name_to_cid_map_lower(names_idx, "title", "PubChem_CID")
                    iupac_map_lc = _build_name_to_cid_map_lower(names_idx, "iupac", "PubChem_CID")
                    syn_map_lc = _build_synonym_to_cid_map_lower(names_idx, synonym_col="synonyms", cid_col="PubChem_CID")

                    title_map_aggr = _build_name_to_cid_map_aggressive(names_idx, "title", "PubChem_CID")
                    iupac_map_aggr = _build_name_to_cid_map_aggressive(names_idx, "iupac", "PubChem_CID")
                    syn_map_aggr = _build_synonym_to_cid_map_aggressive(names_idx, synonym_col="synonyms", cid_col="PubChem_CID")

                    if sqlite_db:
                        all_names = {
                            _lower_name(candidate)
                            for value in s2["Query_Name_toCheck"].dropna().astype("string").tolist()
                            for candidate in _candidate_lookup_names(value)
                            if _lower_name(candidate)
                        } | {
                            _lower_name(candidate)
                            for value in s2["refmet_standardized_name"].dropna().astype("string").tolist()
                            for candidate in _candidate_lookup_names(value)
                            if _lower_name(candidate)
                        }
                        syn_map_lc = {
                            **_fetch_synonym_name_to_cid_sqlite(sqlite_db, all_names),
                            **syn_map_lc,
                        }
                        syn_map_aggr = {
                            **_fetch_synonym_name_to_cid_sqlite_aggressive(
                                sqlite_db,
                                {_aggressive_name_key(name) for name in all_names if _aggressive_name_key(name)},
                            ),
                            **syn_map_aggr,
                        }

                    missing = s2["PubChem_CID"].isna()
                    qname_series = s2.loc[missing, "Query_Name_toCheck"].astype("string")
                    cid_fill = pd.Series(
                        _parallel_map_ordered(
                            qname_series.tolist(),
                            lambda x: _resolve_candidate_cid(
                                x,
                                [title_map_lc, iupac_map_lc, syn_map_lc],
                                [title_map_aggr, iupac_map_aggr, syn_map_aggr],
                            ),
                            max_workers=pubchem_workers,
                        ),
                        index=qname_series.index,
                        dtype="string",
                    )
                    s2.loc[missing, "PubChem_CID"] = cid_fill.astype("string")
                    hit_idx = cid_fill[cid_fill.notna()].index
                    if len(hit_idx) > 0:
                        s2.loc[hit_idx, "matched_name"] = s2.loc[hit_idx, "matched_name"].astype("string").fillna(
                            qname_series.loc[hit_idx].astype("string")
                        )

                    missing = s2["PubChem_CID"].isna()
                    refmet_series = s2.loc[missing, "refmet_standardized_name"].astype("string")
                    rcid_fill = pd.Series(
                        _parallel_map_ordered(
                            refmet_series.tolist(),
                            lambda x: _resolve_candidate_cid(
                                x,
                                [title_map_lc, iupac_map_lc, syn_map_lc],
                                [title_map_aggr, iupac_map_aggr, syn_map_aggr],
                            ),
                            max_workers=pubchem_workers,
                        ),
                        index=refmet_series.index,
                        dtype="string",
                    )
                    s2.loc[missing, "PubChem_CID"] = rcid_fill.astype("string")
                    hit_idx = rcid_fill[rcid_fill.notna()].index
                    if len(hit_idx) > 0:
                        s2.loc[hit_idx, "matched_name"] = s2.loc[hit_idx, "matched_name"].astype("string").fillna(
                            refmet_series.loc[hit_idx].astype("string")
                        )

                    missing = s2["PubChem_CID"].isna()
                    if missing.any():
                        syn_candidate_fill = pd.Series(
                            _parallel_map_ordered(
                                s2.loc[missing, "Query_Name_toCheck"].tolist(),
                                lambda x: next(iter(_fetch_synonym_candidates_sqlite(sqlite_db, x, limit=20)), pd.NA),
                                max_workers=pubchem_workers,
                            ),
                            index=s2.index[missing],
                            dtype="string",
                        )
                        s2.loc[missing, "PubChem_CID"] = syn_candidate_fill.astype("string")
                        hit_idx = syn_candidate_fill[syn_candidate_fill.notna()].index
                        if len(hit_idx) > 0:
                            s2.loc[hit_idx, "matched_name"] = s2.loc[hit_idx, "matched_name"].astype("string").fillna(
                                s2.loc[hit_idx, "Query_Name_toCheck"].astype("string")
                            )

                    missing = s2["PubChem_CID"].isna()
                    if missing.any():
                        refmet_syn_fill = pd.Series(
                            _parallel_map_ordered(
                                s2.loc[missing, "refmet_standardized_name"].tolist(),
                                lambda x: next(iter(_fetch_synonym_candidates_sqlite(sqlite_db, x, limit=20)), pd.NA),
                                max_workers=pubchem_workers,
                            ),
                            index=s2.index[missing],
                            dtype="string",
                        )
                        s2.loc[missing, "PubChem_CID"] = refmet_syn_fill.astype("string")
                        hit_idx = refmet_syn_fill[refmet_syn_fill.notna()].index
                        if len(hit_idx) > 0:
                            s2.loc[hit_idx, "matched_name"] = s2.loc[hit_idx, "matched_name"].astype("string").fillna(
                                s2.loc[hit_idx, "refmet_standardized_name"].astype("string")
                            )

                    s2["PubChem_CID"] = s2["PubChem_CID"].map(_clean_cid).astype("string")

                    title_by_cid = (
                        names_idx[["PubChem_CID", "title"]]
                        .dropna(subset=["PubChem_CID", "title"])
                        .drop_duplicates(subset=["PubChem_CID"], keep="first")
                        .set_index("PubChem_CID")["title"]
                        .to_dict()
                    )
                    iupac_by_cid = (
                        names_idx[["PubChem_CID", "iupac"]]
                        .dropna(subset=["PubChem_CID", "iupac"])
                        .drop_duplicates(subset=["PubChem_CID"], keep="first")
                        .set_index("PubChem_CID")["iupac"]
                        .to_dict()
                    )
                    s2["Standardized_Name"] = _naify(s2["Standardized_Name"])
                    s2["Standardized_Name"] = s2["Standardized_Name"].combine_first(s2["PubChem_CID"].map(title_by_cid))
                    s2["Standardized_Name"] = _naify(s2["Standardized_Name"]).combine_first(s2["PubChem_CID"].map(iupac_by_cid))
                except Exception as e:
                    print(f"[WARN] PubChem names index lookup failed, continuing without it: {e}")

            s2["PubChem_CID"] = s2["PubChem_CID"].map(_clean_cid).astype("string")
            drop_cols = [c for c in ["Query_Name_toCheck", "refmet_standardized_name"] if c in s2.columns]
            if drop_cols:
                s2 = s2.drop(columns=drop_cols)

    with timed_step("Offline SQLite enrichment"):
        cids = {c for c in s2["PubChem_CID"].dropna().tolist() if _clean_cid(c)}
        if sqlite_db:
            smiles_map = _fetch_sqlite_map(sqlite_db, "cid_smiles", "cid", "smiles", cids)
            inchi_map = _fetch_sqlite_map(sqlite_db, "cid_inchi", "cid", "inchikey", cids) if _sqlite_table_exists(sqlite_db, "cid_inchi") else {}
            title_map = _fetch_sqlite_map(sqlite_db, "cid_title", "cid", "title", cids) if _sqlite_table_exists(sqlite_db, "cid_title") else {}
            iupac_map = _fetch_sqlite_map(sqlite_db, "cid_iupac", "cid", "iupac_name", cids) if _sqlite_table_exists(sqlite_db, "cid_iupac") else {}
            synonyms_map = _fetch_synonyms_sqlite(sqlite_db, cids)

            s2["SMILES"] = _naify(s2["SMILES"]).combine_first(s2["PubChem_CID"].map(smiles_map))
            s2["InChIKey"] = _naify(s2["InChIKey"]).combine_first(s2["PubChem_CID"].map(inchi_map))
            s2["Standardized_Name"] = _naify(s2["Standardized_Name"]).combine_first(s2["PubChem_CID"].map(title_map))
            s2["Standardized_Name"] = _naify(s2["Standardized_Name"]).combine_first(s2["PubChem_CID"].map(iupac_map))
            s2["Synonyms"] = _naify(s2["Synonyms"]).combine_first(s2["PubChem_CID"].map(synonyms_map))

    with timed_step("PubChem API enrichment"):
        api = fetch_pubchem_props_api(s2["PubChem_CID"].dropna().tolist(), workers=pubchem_workers)
        if not api.empty:
            s2 = s2.merge(api, on="PubChem_CID", how="left", suffixes=("", "_api"))
            for col in ["Standardized_Name", "Molecular_Formula", "Exact_Mass", "InChIKey", "SMILES"]:
                api_col = f"{col}_api"
                if api_col in s2.columns:
                    s2[col] = _naify(s2[col]).combine_first(_naify(s2[api_col]))
                    s2 = s2.drop(columns=[api_col])

    with timed_step("PubChem PUG-View ID backfill"):
        ann = fetch_pugview_ids_for_cids(
            s2["PubChem_CID"].astype("string").dropna().unique().tolist(),
            sleep=sleep,
            max_workers=pubchem_workers,
        )
        if not ann.empty:
            s2 = s2.merge(ann, on="PubChem_CID", how="left", suffixes=("", "_pv"))
            for c in ["HMDB_ID", "KEGG_ID", "ChEBI_ID"]:
                pv = f"{c}_pv"
                if pv in s2.columns:
                    s2[c] = _naify(s2[c]).combine_first(_naify(s2[pv]))
                    s2 = s2.drop(columns=[pv])

    if not skip_classyfire:
        with timed_step("ClassyFire annotation"):
            cf = annotate_classyfire(
                s2["InChIKey"].astype("string").tolist(),
                r_classyfire_bridge,
                print_output=print_classyfire_output,
            )
            if not cf.empty:
                s2["InChIKey"] = s2["InChIKey"].astype("string").map(_clean_inchikey).astype("string")
                cf["InChIKey"] = cf["InChIKey"].astype("string").map(_clean_inchikey).astype("string")
                s2 = s2.merge(cf, on="InChIKey", how="left", suffixes=("", "_cf"))
                for c in ["Super_Class", "Main_Class", "Sub_Class"]:
                    cc = f"{c}_cf"
                    if cc in s2.columns:
                        s2[c] = _naify(s2[c]).combine_first(_naify(s2[cc]))
                        s2 = s2.drop(columns=[cc])

    s2["ChEBI_ID"] = s2["ChEBI_ID"].astype("string").map(_normalize_chebi).astype("string")
    if "Database_Source" not in s2.columns:
        s2["Database_Source"] = pd.NA
    s2["Database_Source"] = _naify(s2["Database_Source"])
    s2.loc[s2["Database_Source"].isna() & s2["PubChem_CID"].astype("string").map(_clean_cid).notna(), "Database_Source"] = "PubChem"
    s2.loc[s2["Database_Source"].isna() & s2["HMDB_ID"].notna(), "Database_Source"] = "HMDB"
    s2.loc[s2["Database_Source"].isna() & s2["KEGG_ID"].notna(), "Database_Source"] = "KEGG"

    with timed_step("Split final mapped vs unmapped"):
        has_primary_id = _has_primary_mapping_identifier(s2)
        mapped = s2[has_primary_id].copy()
        unmapped = s2[~has_primary_id].copy()
        mapped["matched_name"] = _naify(mapped["matched_name"]).combine_first(mapped["Query_Name"].astype("string"))
        unmapped["Database_Source"] = pd.NA

    mapped = _finalize_public_stage_output(mapped)
    unmapped = _finalize_public_stage_output(unmapped)

    with timed_step("Write combined Stage3 output CSVs"):
        mapped.to_csv(f"{out_prefix}_mapped.csv", index=False)
        unmapped.to_csv(f"{out_prefix}_unmapped.csv", index=False)

    excel_sheets = {
        "mapped": mapped,
        "unmapped": unmapped,
    }
    if not skip_stage5_excel:
        excel_sheets.update(build_stage5_annotation_sheets(mapped))

    with timed_step("Write Stage3 Excel workbook"):
        write_excel_workbook(excel_sheets, f"{out_prefix}.xlsx")

    print(f"[TIMER] TOTAL run_stage3: {time.perf_counter() - total_start:.2f}s")
    return {"mapped": mapped, "unmapped": unmapped}


def main():
    parser = argparse.ArgumentParser(description="HuMANet combined Stage3/Stage4 Python pipeline")
    parser.add_argument("--stage2_unmapped", required=True, help="Stage2 unmapped CSV")
    parser.add_argument("--out_prefix", default="stage3", help="Final combined output prefix")
    parser.add_argument("--pubchem_workers", type=int, default=8, help="Worker count for PubChem property API")
    parser.add_argument("--sleep", type=float, default=0.2, help="Sleep seconds between PubChem PUG-View requests")
    parser.add_argument("--skip_classyfire", action="store_true", help="Disable ClassyFire bridge")
    parser.add_argument("--print_classyfire_output", action="store_true", help="Print ClassyFire rows for debugging")
    parser.add_argument("--disable_sqlite_name_fallback", action="store_true", help="Disable sqlite name->CID fallback if you need stricter parity")
    parser.add_argument("--skip_stage5_excel", action="store_true", help="Skip adding Stage5 sheets to the Excel workbook")
    args = parser.parse_args()

    bridge = CLASSYFIRE_BRIDGE_R if HAS_RPY2 else None

    run_stage3(
        stage2_unmapped_csv=args.stage2_unmapped,
        mw_database_csv=MW_DATABASE_CSV,
        out_prefix=args.out_prefix,
        hmdb_lite_csv=HMDB_LITE_CSV,
        sqlite_db=PUBCHEM_OFFLINE_SQLITE,
        pubchem_workers=args.pubchem_workers,
        r_classyfire_bridge=bridge,
        skip_classyfire=args.skip_classyfire,
        print_classyfire_output=args.print_classyfire_output,
        sleep=args.sleep,
        sqlite_name_fallback=(not args.disable_sqlite_name_fallback),
        mw_unmapped_db_csv=MW_UNMAPPED_DATABASE_CSV,
        skip_stage5_excel=args.skip_stage5_excel,
    )


if __name__ == "__main__":
    main()
