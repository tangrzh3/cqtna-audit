#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 114 -- refresh the C1 resource survey without the eQTL Catalogue API.

DYNAMIC_EQTL_SURVEY.md rests on Step 39, run 2026-08-10, and says the table must
be re-checked before any C1 work starts. The API has returned HTTP 500 on every
data endpoint all session (the api-docs page and the FTP tree are up, so it is the
API backend, not the site). The metadata that endpoint serves is also published as
a static table in the eQTL-Catalogue-resources repository, which is reachable, so
the survey can be refreshed from there instead.

This reproduces Step 39's two passes over the current release:
  1. datasets whose cell type looks like a T cell AND whose condition looks
     stimulated or activated
  2. every T-cell dataset whose condition is not simply naive, in case wording
     differs

and re-applies the Step 39b circularity test: a candidate is flagged if its design
and per-timepoint sample sizes match the Soskic exposure used in this paper.

Outputs: 114a_stim_datasets.tsv, 114b_circularity.tsv, 114_console.log
"""
import io
import sys
import urllib.request

import pandas as pd

MR = r"D:/R_ex/MR"
UA = {"User-Agent": "Mozilla/5.0 (dynamic-eqtl-survey)"}
BASE = ("https://raw.githubusercontent.com/eQTL-Catalogue/"
        "eQTL-Catalogue-resources/master/data_tables/")

# Soskic CD4 activation design used as this paper's exposure
SOSKIC_N = {("naive", "0h"): 99, ("naive", "16h"): 99, ("naive", "40h"): 89,
            ("naive", "5d"): 85, ("memory", "0h"): 100, ("memory", "16h"): 95,
            ("memory", "40h"): 89, ("memory", "5d"): 90}
SOSKIC_SIZES = set(SOSKIC_N.values())

T_WORDS = ("t-cell", "t_cell", "tcell", "cd4", "cd8", "th1", "th2", "th17",
           "treg", "tfh", "lymphocyte")
STIM_WORDS = ("stim", "activ", "anti-cd3", "cd3-cd28", "cd3_cd28", "pha",
              "pma", "iono", "il2", "tcr")


class Tee:
    def __init__(self, p):
        self.f = io.open(p, "w", encoding="utf-8")

    def write(self, s):
        enc = sys.__stdout__.encoding or "utf-8"
        sys.__stdout__.write(s.encode(enc, "replace").decode(enc, "replace"))
        self.f.write(s)

    def flush(self):
        sys.__stdout__.flush()
        self.f.flush()


sys.stdout = Tee(f"{MR}/114_console.log")


def grab(name):
    u = BASE + name
    try:
        raw = urllib.request.urlopen(
            urllib.request.Request(u, headers=UA), timeout=90
        ).read().decode("utf-8", "replace")
        return pd.read_csv(io.StringIO(raw), sep="\t", dtype=str)
    except Exception as ex:
        print(f"  !! {name}: {type(ex).__name__} {getattr(ex,'code','')}")
        return None


print("=" * 78)
print("Step 114 -- C1 resource survey, refreshed from the static metadata table")
print("=" * 78)

md = grab("dataset_metadata_r8.tsv")
if md is None:
    sys.exit("could not fetch dataset metadata")
print(f"dataset_metadata_r8.tsv: {len(md)} rows, {md.shape[1]} columns")
print(f"columns: {list(md.columns)}\n")

low = {c.lower(): c for c in md.columns}


def col(*cands, required=True):
    for c in cands:
        if c in low:
            return low[c]
    if required:
        raise KeyError(f"none of {cands} in {list(md.columns)}")
    return None


C_STUDY = col("study_label", "study")
C_DS = col("dataset_id", "dataset")
C_QM = col("quant_method", "quantification_method", required=False)
C_SG = col("sample_group")
C_TIS = col("tissue_label", "tissue", required=False)
C_CON = col("condition_label", "condition", required=False)
C_N = col("sample_size", "n_samples", required=False)

d = md.copy()
if C_QM:
    d = d[d[C_QM].astype(str).str.lower() == "ge"]
    print(f"gene-expression datasets: {len(d)}")

d["_blob"] = (d[C_SG].fillna("") + " " +
              (d[C_CON].fillna("") if C_CON else "") + " " +
              (d[C_TIS].fillna("") if C_TIS else "")).str.lower()
d["_isT"] = d._blob.apply(lambda b: any(w in b for w in T_WORDS))
d["_isStim"] = d._blob.apply(lambda b: any(w in b for w in STIM_WORDS))

cand = d[d._isT & d._isStim]
print("\n" + "=" * 78)
print(f"pass 1 -- stimulated / activated T-cell datasets: {len(cand)}")
print("=" * 78)
cols = [c for c in (C_STUDY, C_DS, C_SG, C_CON, C_N) if c]
for _, r in cand.iterrows():
    n = r[C_N] if C_N else "?"
    print(f"  {str(r[C_STUDY])[:20]:<22}{str(r[C_DS]):<12}"
          f"{str(r[C_SG])[:34]:<36}{str(r[C_CON])[:22] if C_CON else '':<24}n={n}")

extra = d[d._isT & ~d._isStim]
if C_CON is not None:
    extra = extra[~extra[C_CON].astype(str).str.lower().isin(("naive", "", "nan"))]
print("\n" + "=" * 78)
print(f"pass 2 -- other T-cell datasets, condition not simply naive: {len(extra)}")
print("=" * 78)
for _, r in extra.iterrows():
    n = r[C_N] if C_N else "?"
    print(f"  {str(r[C_STUDY])[:20]:<22}{str(r[C_DS]):<12}"
          f"{str(r[C_SG])[:34]:<36}{str(r[C_CON])[:22] if C_CON else '':<24}n={n}")

allc = pd.concat([cand, extra]).drop_duplicates(subset=[C_DS])
allc[cols].to_csv(f"{MR}/114a_stim_datasets.tsv", sep="\t", index=False)

# ---- circularity against the Soskic exposure ------------------------------
print("\n" + "=" * 78)
print("circularity test (Step 39b rule): does the design match our exposure?")
print("=" * 78)
rows = []
for study, sub in allc.groupby(C_STUDY):
    ns = set()
    if C_N:
        for v in sub[C_N]:
            try:
                ns.add(int(float(v)))
            except (TypeError, ValueError):
                pass
    overlap = ns & SOSKIC_SIZES
    cd4_tp = sub[C_SG].astype(str).str.lower().str.contains("cd4").sum()
    flag = len(overlap) >= 2 and cd4_tp >= 2
    rows.append(dict(study=study, n_datasets=len(sub), sizes=sorted(ns),
                     matching_soskic_sizes=sorted(overlap), circular_flag=flag))
    mark = "  <-- CIRCULAR RISK" if flag else ""
    print(f"  {str(study)[:24]:<26}n={sorted(ns)}"
          f"  matches Soskic {sorted(overlap)}{mark}")

pd.DataFrame(rows).to_csv(f"{MR}/114b_circularity.tsv", sep="\t", index=False)
print("\nSoskic per-timepoint sizes for reference:", sorted(SOSKIC_SIZES))
print("\nwrote 114a_stim_datasets.tsv, 114b_circularity.tsv")
print("\nA flag is a prompt to check the study by hand, not a verdict; two or")
print("more coinciding sample sizes on CD4 arms is suggestive, not proof.")
