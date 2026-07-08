#!/usr/bin/env python3
"""
inspect command  (cli/inspect_libraries.py)
============================================

Build a reviewable "library sheet" from the SRAscraper metadata pickle so a
paired 5' GEX + VDJ (TCR/BCR) dataset can be run through `cellranger multi`.
One row per library (SRR), grouped into patients, tagged Gene Expression / VDJ-T.

Run it AFTER SRAscraper and BEFORE create-config. Exposed two ways:

    ClusterCatcher inspect --metadata-dir <dir> --fastq-base <dir> --out sheet.csv
    python cli/inspect_libraries.py --metadata-dir <dir> --fastq-base <dir> --out sheet.csv

What it does
------------
- Label backfill for runs where the grouping column is NaN (one GSM sequenced
  across multiple runs). Cascade with per-row provenance:
    1) parse experiment_title / experiment_desc  ("GSM12345: <label>")
    2) mode of non-null label within a GSM-level key (library_name, etc.)
    3) run_alias stem  ("GSM12345_r2" -> GSM12345) mapped to a known label
- multi_id derived by a selectable rule (default: patient token = first field of
  the label), so a patient's GEX and VDJ-T collapse into ONE multi group. Keying
  on the full label would split them and defeat `cellranger multi`.
- feature_type / chemistry recomputed from the *filled* label.
- HPV read from a dedicated hpv column (e.g. "hpv state"), with a per-patient
  consensus (GEX-preferred) and a conflict flag.
- Composition + orphan checks: every group should be GEX-only (fine, e.g. 3'
  samples) or GEX+VDJ-T (paired). A VDJ-only group is flagged.

Nothing is downloaded, aligned, or modified. Only a draft CSV is written.
The draft is meant to be reviewed (and edited if needed) before create-config.
"""
import argparse
import os
import re
import pickle
import sys
from types import SimpleNamespace

import click
import numpy as np
import pandas as pd

LIBTYPE_TOKENS = (r"(gex|gene[\s_-]*expression|scrna|rna[\s_-]*seq|transcriptome|"
                  r"tcr|vdj|v\(d\)j|bcr|immune[\s_-]*profiling|5[\s'’-]*prime|3[\s'’-]*prime|"
                  r"5gex|3gex|5p|3p|sc5p|sc3p|five[\s_-]*prime|three[\s_-]*prime)")
NAME_HINTS = ["subject", "patient", "donor", "individual", "participant", "case",
              "title", "alias", "sample_name"]
GSM_KEY_HINTS = ["library_name", "experiment_accession", "sample_accession",
                 "biosample", "geo_accession", "gsm", "experiment_alias"]


# ---------------------------------------------------------------- loaders
def load_pickle_any(path):
    with open(path, "rb") as fh:
        obj = pickle.load(fh)
    frames = []
    if isinstance(obj, dict):
        vals = list(obj.values())
        if vals and hasattr(vals[0], "columns"):            # {GSE: DataFrame}
            for k, v in obj.items():
                if hasattr(v, "columns") and len(v):
                    d = v.copy(); d["__series__"] = k; frames.append(d)
        elif vals and isinstance(vals[0], dict):            # {id: {...}}
            frames.append(pd.DataFrame.from_dict(obj, orient="index").reset_index(names="__key__"))
        else:
            frames.append(pd.DataFrame([obj]))
    elif hasattr(obj, "columns"):
        frames.append(obj.copy())
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def find_run_col(df):
    for c in ["run_accession", "run", "Run", "run_alias", "srr", "SRR"]:
        if c in df.columns:
            return c
    for c in df.columns:
        if df[c].astype(str).str.match(r"^[SED]RR\d+").mean() > 0.8:
            return c
    return None


# ---------------------------------------------------------------- profiling
def profile(df):
    n = len(df)
    print(f"\n{'column':28s} {'dtype':10s} {'n_uniq':>7s} {'n_null':>7s}  examples")
    print("-" * 100)
    for c in df.columns:
        u = df[c].nunique(dropna=True)
        nn = df[c].isna().sum()
        ex = ", ".join(str(x)[:22] for x in df[c].dropna().unique()[:4])
        print(f"{str(c)[:28]:28s} {str(df[c].dtype)[:10]:10s} {u:7d} {nn:7d}  {ex[:60]}")
    print(f"\nrows (libraries/SRRs) = {n}")


def rank_patient_columns(df, run_col):
    n = len(df)
    scored = []
    for c in df.columns:
        if c == run_col or str(c).startswith("__"):
            continue
        u = df[c].nunique(dropna=True)
        if u <= 1 or u > n:
            continue
        spread = 1.0 - abs((u / n) - 0.6)
        hint = 1.0 if any(h in str(c).lower() for h in NAME_HINTS) else 0.0
        scored.append((c, u, round(spread + hint, 3)))
    return sorted(scored, key=lambda x: -x[2])


def detect_gsm_keys(df, run_col):
    """GSM-level keys for backfill: low-null, groups runs (n_unique < n_rows)."""
    n = len(df)
    keys = []
    for c in df.columns:
        if c == run_col or str(c).startswith("__"):
            continue
        u = df[c].nunique(dropna=True)
        nn = df[c].isna().sum()
        if 1 < u < n and nn < 0.5 * n:
            pri = next((i for i, h in enumerate(GSM_KEY_HINTS) if h in str(c).lower()), 99)
            keys.append((c, u, nn, pri))
    keys.sort(key=lambda x: (x[3], x[2], -x[1]))     # name priority, then fewest nulls
    return [(c, u, nn) for c, u, nn, _ in keys]


# ---------------------------------------------------------------- backfill
def strip_gsm_prefix(s):
    return re.sub(r"^\s*GSM\d+\s*[:\-]\s*", "", str(s)).strip()


def backfill_label(df, src, gsm_keys):
    s = df[src].copy()
    prov = ["original" if pd.notna(v) else None for v in s]

    # 1) parse from experiment_title / experiment_desc
    for tcol in ["experiment_title", "experiment_desc"]:
        if tcol in df.columns and s.isna().any():
            for i in s.index:
                if pd.isna(s.iloc[i]):
                    cand = strip_gsm_prefix(df[tcol].iloc[i])
                    if cand and cand.lower() not in ("nan", "none", ""):
                        s.iloc[i] = cand; prov[i] = f"parsed:{tcol}"

    # 2) mode of non-null label within each GSM-level key group
    for key, _, _ in gsm_keys:
        if not s.isna().any():
            break
        grp = df[key]
        mode_map = {}
        for g, sub in s.groupby(grp):
            nn = sub.dropna()
            if not nn.empty:
                mode_map[g] = nn.mode().iloc[0]
        for i in s.index:
            if pd.isna(s.iloc[i]):
                m = mode_map.get(grp.iloc[i])
                if m is not None and pd.notna(m):
                    s.iloc[i] = m; prov[i] = f"mode:{key}"

    # 3) run_alias stem -> known label
    if s.isna().any() and "run_alias" in df.columns:
        stem = df["run_alias"].astype(str).str.replace(r"_r\d+$", "", regex=True)
        label_by_stem = {}
        for i in s.index:
            if pd.notna(s.iloc[i]):
                label_by_stem.setdefault(stem.iloc[i], s.iloc[i])
        for i in s.index:
            if pd.isna(s.iloc[i]):
                m = label_by_stem.get(stem.iloc[i])
                if m:
                    s.iloc[i] = m; prov[i] = "alias_stem"

    return s, prov


# ---------------------------------------------------------------- derivations
def sanitize(x):
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(x)).strip("_")[:60]


def strip_libtokens(x):
    x = re.sub(LIBTYPE_TOKENS, " ", str(x), flags=re.I)
    return re.sub(r"[\s_\-.,#]+", "_", x).strip("_")


def derive_id(label, rule):
    lab = str(label)
    if rule == "token":
        tok = re.split(r"[,_;|/\s]+", lab.strip())
        return sanitize(tok[0]) if tok and tok[0] else sanitize(lab)
    if rule == "stem":
        return sanitize(strip_libtokens(lab))
    return sanitize(lab)


def infer_feature_type(text):
    t = str(text).lower()
    if re.search(r"\bbcr\b|vdj[\s_-]*b|immunoglobulin", t): return "VDJ-B"
    if re.search(r"\btcr\b|vdj|v\(d\)j|t[\s_-]*cell[\s_-]*receptor|immune[\s_-]*profiling", t): return "VDJ-T"
    return "Gene Expression"


def infer_chem(text):
    t = str(text).lower()
    if re.search(r"5[\s'’-]*prime|sc5p|5gex|\b5p\b|five[\s_-]*prime", t): return "5p"
    if re.search(r"3[\s'’-]*prime|sc3p|3gex|\b3p\b|three[\s_-]*prime", t): return "3p"
    return "auto"


def detect_hpv_col(df):
    for c in df.columns:
        cl = str(c).lower()
        if "hpv" in cl or "p16" in cl:
            vals = {str(v).strip().lower() for v in df[c].dropna().unique()}
            if vals & {"+", "-", "positive", "negative", "pos", "neg"}:
                return c
    return None


def hpv_from_val(v):
    t = str(v).strip().lower()
    if t in ("+", "positive", "pos") or t.endswith("positive"): return "HPV+"
    if t in ("-", "negative", "neg") or t.endswith("negative"): return "HPV-"
    return "?"


def hpv_from_label(label):
    t = str(label).lower()
    if "positive" in t: return "HPV+"
    if "negative" in t: return "HPV-"
    return "?"


# ---------------------------------------------------------------- core
def run_inspection(args):
    """Core logic shared by the argparse and click front-ends."""
    pkl = args.pickle
    if not pkl and args.metadata_dir:
        for cand in ["dictionary_file.pkl", "file_metadata.pkl"]:
            p = os.path.join(args.metadata_dir, cand)
            if os.path.exists(p):
                pkl = p; break
    if not pkl or not os.path.exists(pkl):
        sys.exit(f"ERROR: could not find a metadata pickle (pass --pickle). looked in {args.metadata_dir}")
    print(f"metadata pickle: {pkl}")

    df = load_pickle_any(pkl)
    if df.empty:
        sys.exit("ERROR: pickle produced an empty table; inspect it manually.")
    run_col = find_run_col(df)
    print(f"run/SRR column: {run_col}")
    profile(df)

    # ---- choose label (patient) column
    ranked = rank_patient_columns(df, run_col)
    print("\ncandidate label columns (higher score = better):")
    for c, u, s in ranked[:8]:
        print(f"   {str(c):28s} n_unique={u:4d}  score={s}")
    label_col = args.patient_column or (ranked[0][0] if ranked else run_col)

    # ---- GSM-level keys for backfill
    gsm_keys = detect_gsm_keys(df, run_col)
    if args.group_key_column:
        gsm_keys = [(args.group_key_column,
                     df[args.group_key_column].nunique(),
                     df[args.group_key_column].isna().sum())] + \
                   [k for k in gsm_keys if k[0] != args.group_key_column]
    print(f"\nlabel column: '{label_col}'")
    print("GSM-level backfill keys (in order tried):", ", ".join(k[0] for k in gsm_keys[:5]) or "none")

    # ---- backfill label
    label, prov = backfill_label(df, label_col, gsm_keys)
    filled = [(str(df[run_col].iloc[i]), prov[i], str(label.iloc[i])[:40])
              for i in range(len(df)) if prov[i] not in ("original", None)]
    unresolved = [str(df[run_col].iloc[i]) for i in range(len(df)) if prov[i] is None]

    # ---- preview id groupings so the grain is visible
    print("\nmulti_id grouping preview (n groups by rule):")
    for rule in ["token", "stem", "value"]:
        ng = label.map(lambda x: derive_id(x, rule)).nunique()
        star = "  <-- selected" if rule == args.id_from else ""
        print(f"   id-from={rule:6s} -> {ng} groups{star}")

    # ---- feature / chemistry sources (default to filled label)
    feat_src = df[args.feature_column] if (args.feature_column and args.feature_column in df.columns) else label
    chem_src = df[args.chemistry_column] if (args.chemistry_column and args.chemistry_column in df.columns) else label

    # ---- HPV
    hpv_col = args.hpv_column or detect_hpv_col(df)
    print(f"HPV source: {'column ' + repr(hpv_col) if hpv_col else 'label token (no HPV column found)'}")

    series_col = "__series__" if "__series__" in df.columns else None
    rows = []
    for i in range(len(df)):
        srr = str(df[run_col].iloc[i])
        series = str(df[series_col].iloc[i]) if series_col else ""
        fq = os.path.join(args.fastq_base, series, srr) if series else os.path.join(args.fastq_base, srr)
        hpv = hpv_from_val(df[hpv_col].iloc[i]) if hpv_col else hpv_from_label(label.iloc[i])
        # source_title kept clean: experiment_title-derived labels can carry
        # trailing "; Homo sapiens; OTHER" fields, so trim at the first ';'.
        title_clean = str(label.iloc[i]).split(";")[0].strip()[:60]
        rows.append({
            "multi_id": derive_id(label.iloc[i], args.id_from),
            "fastq_id": srr,
            "feature_type": infer_feature_type(feat_src.iloc[i]),
            "chemistry": infer_chem(chem_src.iloc[i]),
            "fastqs": fq,
            "hpv": hpv,
            "source_title": title_clean,
            "_prov": prov[i] or "UNRESOLVED",
        })
    sheet = pd.DataFrame(rows)

    # ---- per-patient HPV consensus (GEX-preferred) + conflict flags
    hpv_conflicts = []
    for mid, sub in sheet.groupby("multi_id"):
        setv = {v for v in sub["hpv"] if v in ("HPV+", "HPV-")}
        gex = sub[sub.feature_type == "Gene Expression"]
        if len(gex) and gex["hpv"].iloc[0] in ("HPV+", "HPV-"):
            cons = gex["hpv"].iloc[0]
        elif setv:
            cons = sorted(setv)[0]
        else:
            cons = "?"
        if len(setv) > 1:
            hpv_conflicts.append((mid, "/".join(sorted(setv)), f"-> {cons} (GEX-preferred)"))
        sheet.loc[sheet.multi_id == mid, "hpv"] = cons

    sheet = sheet.sort_values(["multi_id", "feature_type"]).reset_index(drop=True)
    sheet.drop(columns=["_prov"]).to_csv(args.out, index=False)

    # ---- report
    print(f"\ndraft library sheet written: {args.out}  "
          f"({len(sheet)} libraries, {sheet['multi_id'].nunique()} groups)")
    print(sheet.drop(columns=["_prov"]).to_string(index=False))

    print("\n--- BACKFILL PROVENANCE ---")
    if filled:
        for srr, src, lab in filled:
            print(f"   {srr}  filled via {src:22s} -> {lab}")
    else:
        print("   (no rows needed backfill)")
    if unresolved:
        print(f"   !! UNRESOLVED (no label recoverable): {', '.join(unresolved)}")

    print("\n--- GROUP COMPOSITION ---")
    for mid, sub in sheet.groupby("multi_id"):
        comp = "+".join(sorted(set(sub["feature_type"])))
        chem = "/".join(sorted(set(sub["chemistry"])))
        n = len(sub)
        flag = ""
        if "Gene Expression" not in comp:
            flag = "  <-- WARNING: VDJ with no GEX (orphan)"
        elif "VDJ-T" in comp and "5p" not in chem and "auto" not in chem:
            flag = "  <-- WARNING: VDJ paired with non-5' chemistry"
        print(f"   {mid:16s} {comp:28s} n={n} chem={chem}{flag}")

    if hpv_conflicts:
        print("\n--- HPV LABEL CONFLICTS (resolved, review) ---")
        for mid, seen, res in hpv_conflicts:
            print(f"   {mid:16s} title said {seen}  {res}")

    print("\nReview the sheet. If grain or feature_type is wrong, rerun with")
    print("  --patient-column / --group-key-column / --id-from / --feature-column / --hpv-column.")
    print("A paired 5' patient = one 'Gene Expression' row + one 'VDJ-T' row sharing multi_id.")
    return sheet


# ---------------------------------------------------------------- argparse front-end
def main():
    ap = argparse.ArgumentParser(description="Build a cellranger multi library sheet from SRAscraper metadata.")
    ap.add_argument("--pickle", help="Path to the metadata pickle (auto-located if omitted)")
    ap.add_argument("--metadata-dir", help="SRAscraper metadata dir (dictionary_file.pkl / file_metadata.pkl)")
    ap.add_argument("--fastq-base", required=True, help="Base dir: {fastq_base}/{SERIES}/{SRR}/*.fastq.gz")
    ap.add_argument("--out", default="./library_sheet_draft.csv")
    ap.add_argument("--patient-column", help="Column holding the label (default: auto)")
    ap.add_argument("--group-key-column", help="GSM-level key for backfill (default: auto)")
    ap.add_argument("--id-from", choices=["token", "stem", "value"], default="token",
                    help="How multi_id derives from the label (default: token = first field)")
    ap.add_argument("--feature-column", help="Column used to infer feature_type (default: label)")
    ap.add_argument("--chemistry-column", help="Column used to infer chemistry (default: label)")
    ap.add_argument("--hpv-column", help="Column holding HPV status (default: auto-detect)")
    run_inspection(ap.parse_args())


# ---------------------------------------------------------------- click front-end
@click.command('inspect')
@click.option('--pickle', 'pickle_path', default=None,
              help='Path to the SRAscraper metadata pickle (auto-located from --metadata-dir if omitted)')
@click.option('--metadata-dir', default=None,
              help='SRAscraper metadata dir (looks for dictionary_file.pkl / file_metadata.pkl)')
@click.option('--fastq-base', required=True,
              help='Base FASTQ dir; rows point at {fastq_base}/{SERIES}/{SRR}/')
@click.option('--out', default='./library_sheet_draft.csv',
              help='Draft library sheet output path (default: ./library_sheet_draft.csv)')
@click.option('--patient-column', default=None, help='Force the label/patient column (default: auto-rank)')
@click.option('--group-key-column', default=None, help='Force the GSM-level backfill key (default: auto-detect)')
@click.option('--id-from', type=click.Choice(['token', 'stem', 'value']), default='token',
              help="How multi_id derives from the label: token=first field (default), stem=strip lib tokens, value=whole")
@click.option('--feature-column', default=None, help='Column used to infer feature_type (default: label)')
@click.option('--chemistry-column', default=None, help='Column used to infer chemistry (default: label)')
@click.option('--hpv-column', default=None, help='Column holding HPV status (default: auto-detect e.g. "hpv state")')
def inspect_libraries(pickle_path, metadata_dir, fastq_base, out, patient_column,
                      group_key_column, id_from, feature_column, chemistry_column, hpv_column):
    """
    Build a reviewable cellranger multi library sheet from SRAscraper metadata.

    Run after SRAscraper and before create-config for paired 5' GEX + VDJ
    (TCR/BCR) datasets. Profiles the metadata pickle, groups libraries into
    patients, tags Gene Expression vs VDJ-T, backfills labels for multi-run
    GSMs, and writes a draft CSV for you to review before create-config.

    \b
    Example:
      ClusterCatcher inspect \\
          --metadata-dir /path/to/GSE######/metadata \\
          --fastq-base   /path/to/GSE######/fastq \\
          --out          ./library_sheet_draft.csv
    """
    args = SimpleNamespace(
        pickle=pickle_path,
        metadata_dir=metadata_dir,
        fastq_base=fastq_base,
        out=out,
        patient_column=patient_column,
        group_key_column=group_key_column,
        id_from=id_from,
        feature_column=feature_column,
        chemistry_column=chemistry_column,
        hpv_column=hpv_column,
    )
    run_inspection(args)


if __name__ == '__main__':
    main()
