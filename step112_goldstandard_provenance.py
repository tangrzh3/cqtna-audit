#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 112 -- S34 §3.1 criterion 6 on the class-B candidates.

Step 110 disqualified FLAMES, cS2G and PEGASUS as truth sets: they are effector-
gene PREDICTORS that take eQTL evidence as input, which is the very thing this
paper audits. What may still qualify is the gold-standard set each was EVALUATED
against, which is a different object.

The question this script serves is S34 §3.1 criterion 6: is that gold standard's
gene assignment decided independently of eQTL evidence? If every candidate gold
standard also rests on eQTL, the check lands on §5 branch F and is not run.

Only the passages describing how each gold standard was BUILT are pulled. No
locus-to-gene table is downloaded, per S34 §3.2.

Outputs: 112a_goldstandard_passages.txt, 112_console.log
"""
import os
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
UA = {"User-Agent": "Mozilla/5.0 (attribution-benchmark-survey)"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CANDIDATES = [
    ("cS2G", "35668300", "10.1038/s41588-022-01087-y"),
    ("FLAMES", "39930082", "10.1038/s41588-025-02084-7"),
]

# terms that mark a passage describing how the gold standard was assembled
MARK = re.compile(
    r"gold[- ]standard|silver[- ]standard|benchmark set|ground truth|"
    r"curated set|positive control gene|training set|known causal gene",
    re.I)
# terms that would disqualify under criterion 6
EQTL = re.compile(r"\beQTL|expression quantitative|colocali[sz]|"
                  r"TWAS|SMR\b|S2G|SNP-to-gene", re.I)


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


sys.stdout = Tee(f"{MR}/112_console.log")


def get(url, parse=None, tries=4):
    for a in range(tries):
        try:
            r = urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=90)
            return json.load(r) if parse == "json" else r.read().decode(
                "utf-8", "replace")
        except Exception as ex:
            if a == tries - 1:
                print(f"    !! {type(ex).__name__} {getattr(ex,'code','')}")
                return None
            time.sleep(2 * (a + 1))


def pmcid_for(pmid):
    q = urllib.parse.urlencode(dict(dbfrom="pubmed", db="pmc", id=pmid,
                                    retmode="json"))
    r = get(f"{EUTILS}/elink.fcgi?{q}", parse="json")
    try:
        return r["linksets"][0]["linksetdbs"][0]["links"][0]
    except Exception:
        return None


def pmc_text(pmcid):
    q = urllib.parse.urlencode(dict(db="pmc", id=pmcid, rettype="full",
                                    retmode="xml"))
    xml = get(f"{EUTILS}/efetch.fcgi?{q}")
    if not xml:
        return None
    txt = re.sub(r"<[^>]+>", " ", xml)
    txt = re.sub(r"&[a-z]+;", " ", txt)
    return re.sub(r"[ \t]+", " ", txt)


print("=" * 78)
print("Step 112 -- how was each candidate gold standard built?")
print("=" * 78)
print("S34 §3.1 criterion 6: the assignment must not rest on eQTL evidence.\n")

out = io.open(f"{MR}/112a_goldstandard_passages.txt", "w", encoding="utf-8")

for name, pmid, doi in CANDIDATES:
    print("=" * 78)
    print(f"{name}   pmid {pmid}   doi {doi}")
    print("=" * 78)
    pmcid = pmcid_for(pmid)
    print(f"  PMC id: {pmcid or 'NOT IN PMC -- full text unavailable here'}")
    if not pmcid:
        print("  -> criterion 6 cannot be judged from open full text.\n")
        out.write(f"\n\n### {name}: not in PMC\n")
        continue

    txt = pmc_text(pmcid)
    if not txt:
        print("  -> fetch failed\n")
        continue
    print(f"  full text: {len(txt):,} chars")

    sents = re.split(r"(?<=[.;])\s+", txt)
    hits = [s.strip() for s in sents if MARK.search(s) and 60 < len(s) < 900]
    seen, keep = set(), []
    for s in hits:
        k = s[:80]
        if k in seen:
            continue
        seen.add(k)
        keep.append(s)

    print(f"  passages mentioning a gold/benchmark set: {len(keep)}")
    out.write(f"\n\n{'='*78}\n### {name}  pmid:{pmid} doi:{doi}\n{'='*78}\n")
    for s in keep:
        out.write(s + "\n\n")

    flagged = [s for s in keep if EQTL.search(s)]
    print(f"  of which ALSO mention eQTL/S2G/coloc/TWAS: {len(flagged)}")
    print("\n  --- first passages describing the set ---")
    for s in keep[:6]:
        mark = " [eQTL-linked]" if EQTL.search(s) else ""
        print(f"    * {s[:300]}{mark}")
    print()
    time.sleep(0.5)

out.close()
print("=" * 78)
print("wrote 112a_goldstandard_passages.txt")
print("Judgement against criterion 6 is recorded by hand in S34 §9.1.")
print("No locus-to-gene table has been downloaded.")
print("=" * 78)
