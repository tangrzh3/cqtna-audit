#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 117 -- is there any large multi-condition stimulation eQTL resource?

The Cano-Gamez suggestion exposed a gap: studies with several cytokine-polarising
conditions have single-digit donor counts, and the stimulated resources that do
carry eQTLs have one or few conditions. This screens every gene-expression dataset
in the catalogue against the three criteria that would close it:

    1. it is an eQTL dataset (all catalogue entries are)
    2. donors >= 100, the rough floor for cis-eQTL mapping
    3. its study carries more than one non-resting condition

Metadata only; no summary statistics are downloaded.

Outputs: 117a_multicondition_studies.tsv, 117_console.log
"""
import io
import sys
import urllib.request

import pandas as pd

MR = r"D:/R_ex/MR"
UA = {"User-Agent": "Mozilla/5.0 (condition-scan)"}
META = ("https://raw.githubusercontent.com/eQTL-Catalogue/"
        "eQTL-Catalogue-resources/master/data_tables/dataset_metadata_r8.tsv")
MIN_N = 100

RESTING = {"naive", "", "nan", "none", "control", "unstimulated", "untreated"}


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


sys.stdout = Tee(f"{MR}/117_console.log")

raw = urllib.request.urlopen(
    urllib.request.Request(META, headers=UA), timeout=90
).read().decode("utf-8", "replace")
md = pd.read_csv(io.StringIO(raw), sep="\t", dtype=str)
d = md[md.quant_method == "ge"].copy()
d["n"] = pd.to_numeric(d.sample_size, errors="coerce")
d["cond"] = d.condition_label.fillna("").astype(str)
d["is_stim"] = ~d.cond.str.lower().isin(RESTING)

print("=" * 78)
print("Step 117 -- multi-condition eQTL resources with usable donor numbers")
print("=" * 78)
print(f"gene-expression datasets: {len(d)}   studies: {d.study_label.nunique()}")

rows = []
for study, sub in d.groupby("study_label"):
    stim = sub[sub.is_stim]
    conds = sorted(set(stim.cond))
    big = stim[stim.n >= MIN_N]
    rows.append(dict(study=study,
                     n_datasets=len(sub),
                     n_stim_datasets=len(stim),
                     n_conditions=len(conds),
                     max_n=int(sub.n.max()) if sub.n.notna().any() else None,
                     max_n_stim=int(stim.n.max()) if stim.n.notna().any() else None,
                     n_stim_datasets_ge100=len(big),
                     conditions="; ".join(conds)[:160]))

s = pd.DataFrame(rows).sort_values(["n_conditions", "max_n_stim"],
                                   ascending=False)
s.to_csv(f"{MR}/117a_multicondition_studies.tsv", sep="\t", index=False)

print("\n" + "=" * 78)
print(f"studies with >1 non-resting condition (any donor count)")
print("=" * 78)
multi = s[s.n_conditions > 1]
for _, r in multi.iterrows():
    flag = "  <== MEETS ALL THREE" if r.n_stim_datasets_ge100 > 0 else ""
    print(f"  {str(r.study)[:26]:<28}conditions={r.n_conditions:<3}"
          f"max n(stim)={str(r.max_n_stim):<6}"
          f"stim datasets n>={MIN_N}: {r.n_stim_datasets_ge100}{flag}")
    print(f"      {r.conditions}")

print("\n" + "=" * 78)
print(f"criterion 3 alone: studies with >1 condition = {len(multi)}")
hit = multi[multi.n_stim_datasets_ge100 > 0]
print(f"all three criteria (>1 condition AND a stimulated dataset at n>={MIN_N}): "
      f"{len(hit)}")
print("=" * 78)
if len(hit):
    for _, r in hit.iterrows():
        print(f"  {r.study}: {r.n_conditions} conditions, "
              f"{r.n_stim_datasets_ge100} stimulated datasets at n>={MIN_N}, "
              f"max n {r.max_n_stim}")
else:
    print("  none -- the gap the Cano-Gamez suggestion pointed at is real")

print("\nwrote 117a_multicondition_studies.tsv")
