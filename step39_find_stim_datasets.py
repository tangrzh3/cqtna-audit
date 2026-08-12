#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 39  Is there any other eQTL dataset with a stimulated CD4 T-cell condition
         in which TPI1 is actually quantified?

DICE could not test the activation-window hypothesis: TPI1 is absent from both
of its stimulated datasets (Step 38d), despite median TPM 98-194 where it is
measured. Before concluding that no independent replication is possible, sweep
the whole eQTL Catalogue.

Two passes:
  1. every dataset whose tissue/condition looks like a stimulated or activated
     T cell (or any CD4 population), regardless of study
  2. for each candidate, ask whether TPI1 and SPSB2 are quantified at all --
     the DICE failure mode -- before counting it as usable
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


def get(path, **p):
    url = f"{API}/{path}?" + urllib.parse.urlencode(p)
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


ds, _ = get("datasets/", size=1000)
print(f"catalogue datasets: {len(ds)}\n")

# ---- pass 1: candidate conditions ---------------------------------------
T_WORDS = ("t-cell", "t_cell", "tcell", "cd4", "cd8", "th1", "th2", "th17",
           "treg", "tfh", "lymphocyte")
STIM_WORDS = ("stim", "activ", "anti-cd3", "cd3-cd28", "cd3_cd28", "pha",
              "pma", "iono", "il2", "tcr")

cands = []
for d in ds:
    if d["quant_method"] != "ge":
        continue
    blob = (f"{d['sample_group']} {d['condition_label']} "
            f"{d['tissue_label']}").lower()
    is_t = any(w in blob for w in T_WORDS)
    is_stim = any(w in blob for w in STIM_WORDS)
    if is_t and is_stim:
        cands.append(d)

print("=" * 96)
print(f"stimulated / activated T-cell datasets (gene expression): {len(cands)}")
print("=" * 96)
for d in cands:
    print(f"  {d['study_label']:<22}{d['dataset_id']:<11}{d['sample_group']:<34}"
          f"{d['condition_label']:<24}n={d['sample_size']}")

# also list every non-naive condition, in case wording differs
print("\n" + "=" * 96)
print("all T-cell datasets whose condition is not simply 'naive'")
print("=" * 96)
extra = [d for d in ds if d["quant_method"] == "ge"
         and any(w in f"{d['sample_group']} {d['tissue_label']}".lower()
                 for w in T_WORDS)
         and d["condition_label"].lower() not in ("naive", "")]
for d in extra:
    print(f"  {d['study_label']:<22}{d['dataset_id']:<11}{d['sample_group']:<34}"
          f"{d['condition_label']:<24}n={d['sample_size']}")

# ---- pass 2: is TPI1 actually quantified there? -------------------------
check = {d["dataset_id"]: d for d in cands + extra}
print("\n" + "=" * 96)
print("gene coverage check (the DICE failure mode)")
print("=" * 96)
rows = []
for did, d in check.items():
    rec = dict(study=d["study_label"], dataset=did, cell=d["sample_group"],
               condition=d["condition_label"], n=int(d["sample_size"]))
    for tag, ensg in (("TPI1", TPI1), ("SPSB2", SPSB2)):
        a, code = get(f"datasets/{did}/associations",
                      molecular_trait_id=ensg, size=1)
        time.sleep(0.35)
        rec[tag] = bool(a)
        rec[f"{tag}_tpm"] = a[0].get("median_tpm") if a else None
    rows.append(rec)
    print(f"  {d['study_label']:<20}{d['sample_group']:<32}"
          f"{d['condition_label']:<22}TPI1={rec['TPI1']!s:<6}SPSB2={rec['SPSB2']}",
          flush=True)

t = pd.DataFrame(rows)
t.to_csv(f"{MR}/39a_stimulated_Tcell_datasets.tsv", sep="\t", index=False)
usable = t[t.TPI1]
print("\n" + "=" * 96)
print(f"USABLE (TPI1 quantified): {len(usable)} of {len(t)}")
print("=" * 96)
if len(usable):
    print(usable[["study", "dataset", "cell", "condition", "n", "TPI1_tpm"]]
          .to_string(index=False))
else:
    print("  none -- no independent dataset can test the activation-window "
          "hypothesis for TPI1")
print("\nwritten: 39a")
