#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 39b  CIRCULARITY CHECK before anything else.

Cytoimmgen in the eQTL Catalogue exposes combined_CD4_Naive_STIM_16H (n=99),
combined_CD4_Memory_STIM_16H (n=95), plus 40H and 5D arms -- the same design and,
for 16 h, the same sample sizes as the Soskic CD4 activation data used as our
EXPOSURE (Naive 0h=99, 16h=99, 40h=89, 5d=85; Memory 0h=100, 16h=95, 40h=89,
5d=90).

If Cytoimmgen is that same experiment, then "replicating" our activation-window
finding in it would be circular and must not be done. This script gathers the
study metadata needed to decide, and separately identifies which stimulated
T-cell datasets come from genuinely independent studies AND quantify TPI1.
"""
import json
import time
import urllib.request
import urllib.parse
import pandas as pd

API = "https://www.ebi.ac.uk/eqtl/api/v2"
UA = {"User-Agent": "Mozilla/5.0"}
MR = r"D:/R_ex/MR"
TPI1, SPSB2 = "ENSG00000111669", "ENSG00000111671"

SOSKIC_N = {("Naive", "0h"): 99, ("Naive", "16h"): 99, ("Naive", "40h"): 89,
            ("Naive", "5d"): 85, ("Memory", "0h"): 100, ("Memory", "16h"): 95,
            ("Memory", "40h"): 89, ("Memory", "5d"): 90}


def get(path, **p):
    url = f"{API}/{path}" + ("?" + urllib.parse.urlencode(p) if p else "")
    for a in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=120) as r:
                return json.load(r), None
        except Exception as ex:
            code = getattr(ex, "code", None)
            if a == 4:
                return None, code
            time.sleep(2 * (a + 1))
    return None, None


# ---- study-level metadata -------------------------------------------------
print("=" * 88)
print("study metadata")
print("=" * 88)
studies, _ = get("studies/", size=1000)
if studies:
    for s in studies:
        if s.get("study_label") in ("Cytoimmgen", "Schmiedel_2018", "Nathan_2022",
                                    "Soskic_2022", "Gilchrist_2021"):
            print(f"  {s}")

ds, _ = get("datasets/", size=1000)
cyto = [d for d in ds if d["study_label"] == "Cytoimmgen" and d["quant_method"] == "ge"]
print(f"\nCytoimmgen gene-expression datasets: {len(cyto)}")

print("\n" + "=" * 88)
print("sample-size comparison with the Soskic exposure data")
print("=" * 88)
for d in sorted(cyto, key=lambda x: x["sample_group"]):
    sg = d["sample_group"]
    if not sg.startswith("combined_CD4"):
        continue
    lin = "Naive" if "Naive" in sg else "Memory"
    tp = ("16h" if "16H" in sg else "40h" if "40H" in sg
          else "5d" if "5D" in sg else "0h")
    exp = SOSKIC_N.get((lin, tp))
    flag = "  <-- SAME N" if exp == int(d["sample_size"]) else ""
    print(f"  {sg:<34}n={d['sample_size']:<5}Soskic {lin} {tp} n={exp}{flag}")

# ---- gene coverage in the independent candidates --------------------------
print("\n" + "=" * 88)
print("gene coverage in candidate stimulated T-cell datasets")
print("=" * 88)
T_WORDS = ("t-cell", "cd4", "cd8", "th1", "th2", "th17", "treg", "tfh", "tn",
           "tcm", "tem", "temra")
STIM = ("stim", "activ", "anti-cd3")
cands = [d for d in ds if d["quant_method"] == "ge"
         and any(w in f"{d['sample_group']} {d['tissue_label']}".lower() for w in T_WORDS)
         and any(w in f"{d['sample_group']} {d['condition_label']}".lower() for w in STIM)]
by_study = {}
for d in cands:
    by_study.setdefault(d["study_label"], []).append(d)

rows = []
for study, lst in by_study.items():
    # test at most 4 datasets per study to keep the query count sane
    for d in lst[:4]:
        rec = dict(study=study, dataset=d["dataset_id"], cell=d["sample_group"],
                   condition=d["condition_label"], n=int(d["sample_size"]))
        for tag, ensg in (("TPI1", TPI1), ("SPSB2", SPSB2)):
            a, code = get(f"datasets/{d['dataset_id']}/associations",
                          molecular_trait_id=ensg, size=1)
            time.sleep(0.35)
            rec[tag] = bool(a)
            rec[f"{tag}_tpm"] = a[0].get("median_tpm") if a else None
        rows.append(rec)
        print(f"  {study:<16}{d['sample_group']:<34}{d['condition_label']:<20}"
              f"TPI1={rec['TPI1']!s:<6}SPSB2={rec['SPSB2']}", flush=True)

t = pd.DataFrame(rows)
t.to_csv(f"{MR}/39b_stim_dataset_coverage.tsv", sep="\t", index=False)
print("\n" + "=" * 88)
print("summary")
print("=" * 88)
print(t.groupby("study").agg(datasets=("dataset", "size"),
                             TPI1_quantified=("TPI1", "sum"),
                             SPSB2_quantified=("SPSB2", "sum")).to_string())
print("\nwritten: 39b")
