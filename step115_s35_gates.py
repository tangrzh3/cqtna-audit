#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 115 -- S35 §3.2 gate 1: is each new resource independent of our exposure?

S35 requires a hand check of study design, donor source and sample sizes before
any cell is run, and says explicitly that step114's automatic flag is a prompt
rather than a verdict. This script gathers the evidence for that judgement: the
catalogue's own study metadata, and each study's abstract and methods-level
description of its donors and design.

It reads no eQTL statistic. Instrument extraction is step 116 and only proceeds
for resources that pass this gate.

Outputs: 115a_study_evidence.txt, 115b_gate1.tsv, 115_console.log
"""
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

MR = r"D:/R_ex/MR"
UA = {"User-Agent": "Mozilla/5.0 (s35-gate-check)"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
META = ("https://raw.githubusercontent.com/eQTL-Catalogue/"
        "eQTL-Catalogue-resources/master/data_tables/dataset_metadata_r8.tsv")

CELLS = {           # resource -> dataset ids registered in S35
    "Nathan_2022": ["QTD000666"],
    "Randolph_2021": ["QTD000588"],
    "Schmiedel_2018": ["QTD000484"],
}
# Soskic design, for the comparison the gate is about
SOSKIC = ("naive and memory CD4 at 0 h, 16 h, 40 h and 5 d; 85-100 donors per "
          "profile; anti-CD3/CD28; Nature Genetics 2022")


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


sys.stdout = Tee(f"{MR}/115_console.log")


def get(url, as_json=False, tries=4):
    for a in range(tries):
        try:
            r = urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=90)
            return json.load(r) if as_json else r.read().decode("utf-8",
                                                                "replace")
        except Exception as ex:
            if a == tries - 1:
                print(f"    !! {type(ex).__name__} {getattr(ex,'code','')}")
                return None
            time.sleep(2 * (a + 1))


def abstract(pmid):
    q = urllib.parse.urlencode(dict(db="pubmed", id=pmid, rettype="abstract",
                                    retmode="text"))
    time.sleep(0.3)
    return get(f"{EUTILS}/efetch.fcgi?{q}")


print("=" * 78)
print("Step 115 -- S35 gate 1: circularity hand check")
print("=" * 78)
print(f"our exposure: {SOSKIC}\n")

raw = get(META)
md = pd.read_csv(io.StringIO(raw), sep="\t", dtype=str)

out = io.open(f"{MR}/115a_study_evidence.txt", "w", encoding="utf-8")
rows = []

for study, dsets in CELLS.items():
    sub = md[md.study_label == study]
    reg = sub[sub.dataset_id.isin(dsets)]
    pmid = sorted(set(sub.pmid.dropna()))
    print("=" * 78)
    print(f"{study}")
    print("=" * 78)
    print(f"  datasets in catalogue: {len(sub)}   registered for S35: {dsets}")
    print(f"  pmid: {pmid}")
    for _, r in reg.iterrows():
        print(f"    {r.dataset_id}  {r.sample_group:<26} "
              f"cond={r.condition_label:<20} n={r.sample_size} "
              f"type={r.study_type}")

    sizes = sorted({int(float(v)) for v in sub.sample_size.dropna()
                    if str(v).replace('.', '').isdigit()})
    print(f"  all sample sizes in this study: {sizes}")

    out.write(f"\n\n{'='*78}\n### {study}  pmid={pmid}\n{'='*78}\n")
    ab = ""
    if pmid:
        ab = abstract(pmid[0]) or ""
        out.write(ab)

    # surface the sentences that speak to donors and design
    DON = re.compile(r"\b(donor|individual|participant|volunteer|subject|"
                     r"sample[sd]? from|cohort|recruit)\w*\b", re.I)
    sents = [s.strip() for s in re.split(r"(?<=[.;])\s+", ab)
             if DON.search(s) and 40 < len(s) < 500]
    print("  --- donor / cohort statements from the abstract:")
    if not sents:
        print("      (none found in abstract; full text needed for a firm call)")
    for s in sents[:4]:
        print(f"      * {s[:260]}")

    rows.append(dict(study=study, datasets=";".join(dsets), pmid=";".join(pmid),
                     n_datasets_in_study=len(sub), sizes=str(sizes),
                     donor_sentences=len(sents)))
    print()

pd.DataFrame(rows).to_csv(f"{MR}/115b_gate1.tsv", sep="\t", index=False)
out.close()

print("=" * 78)
print("wrote 115a_study_evidence.txt, 115b_gate1.tsv")
print("The verdict per resource is a hand judgement to be entered in S35 §10.1.")
print("No eQTL statistic has been read.")
print("=" * 78)
