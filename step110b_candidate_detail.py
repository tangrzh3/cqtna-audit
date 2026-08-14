#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 110b -- broader candidate search plus abstracts, for S34 §3.1 scoring.

Round 1 (step110a) used queries whose AND/OR precedence was too tight and three
returned nothing. This round uses plain term queries, and fetches ABSTRACTS for
the promising records so the six admissibility criteria can be judged.

Abstracts are metadata. No locus-to-gene table is opened here, in keeping with
S34 §3.2: candidates must be scored before their content is seen.

Outputs: 110c_candidates_round2.tsv, 110d_abstracts.txt, 110b_console.log
"""
import io
import json
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

MR = r"D:/R_ex/MR"
UA = {"User-Agent": "Mozilla/5.0 (attribution-benchmark-survey)"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


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


sys.stdout = Tee(f"{MR}/110b_console.log")


def fetch(url, parse=json.load, tries=4):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            r = urllib.request.urlopen(req, timeout=60)
            return parse(r) if parse else r.read().decode("utf-8", "replace")
        except Exception as ex:
            if a == tries - 1:
                print(f"    !! {type(ex).__name__} {getattr(ex,'code','')}")
                return None
            time.sleep(2 * (a + 1))


def pubmed(query, n=10):
    q = urllib.parse.urlencode(dict(db="pubmed", term=query, retmax=n,
                                    retmode="json", sort="relevance"))
    res = fetch(f"{EUTILS}/esearch.fcgi?{q}")
    ids = (res or {}).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    q2 = urllib.parse.urlencode(dict(db="pubmed", id=",".join(ids),
                                     retmode="json"))
    summ = fetch(f"{EUTILS}/esummary.fcgi?{q2}") or {}
    out = []
    for pmid in ids:
        r = summ.get("result", {}).get(pmid)
        if not r:
            continue
        doi = next((a["value"] for a in r.get("articleids", [])
                    if a.get("idtype") == "doi"), "")
        out.append(dict(pmid=pmid, title=r.get("title", ""),
                        journal=r.get("fulljournalname", ""),
                        year=(r.get("pubdate", "") or "")[:4], doi=doi))
    time.sleep(0.2)
    return out


def abstract(pmid):
    q = urllib.parse.urlencode(dict(db="pubmed", id=pmid, rettype="abstract",
                                    retmode="text"))
    time.sleep(0.2)
    return fetch(f"{EUTILS}/efetch.fcgi?{q}", parse=None)


QUERIES = [
    "Prioritizing effector genes trait-associated loci multimodal",
    "open targets genetics locus2gene machine learning gold standard",
    "combined SNP-to-gene strategy cS2G",
    "effector gene GWAS locus benchmark evaluation",
    "gold standard GWAS loci known causal genes curated list",
    "predicted effector gene aggregation standards unified schema",
]

print("=" * 78)
print("Step 110b -- broader candidate search (METADATA + ABSTRACTS ONLY)")
print("=" * 78)

rows, seen = [], {}
for q in QUERIES:
    print(f"\n--- {q}")
    hits = pubmed(q)
    if not hits:
        print("    (no results)")
    for h in hits:
        if h["pmid"] in seen:
            continue
        seen[h["pmid"]] = h
        h["query"] = q
        rows.append(h)
        print(f"    {h['year']:<6}{h['journal'][:28]:<30}{h['title'][:66]}")

if rows:
    pd.DataFrame(rows)[["year", "journal", "title", "doi", "pmid", "query"]] \
        .to_csv(f"{MR}/110c_candidates_round2.tsv", sep="\t", index=False)
    print(f"\nwrote 110c_candidates_round2.tsv  ({len(rows)} unique)")

# ---- abstracts for the records most likely to define a truth set ----------
KEY = ("effector gene", "gold standard", "benchmark", "causal gene",
       "locus-to-gene", "locus to gene", "prioriti")
picks = [h for h in rows
         if any(k in h["title"].lower() for k in KEY)][:10]

# carry over the strongest round-1 records by PMID
for pmid in ("34365404", "40000000"):
    pass

print("\n" + "=" * 78)
print(f"abstracts for {len(picks)} candidate records")
print("=" * 78)
with io.open(f"{MR}/110d_abstracts.txt", "w", encoding="utf-8") as fh:
    for h in picks:
        print(f"\n### {h['year']} {h['journal']}")
        print(f"### {h['title']}")
        print(f"### doi:{h['doi']}  pmid:{h['pmid']}")
        ab = abstract(h["pmid"]) or "(abstract unavailable)"
        fh.write(f"\n{'='*78}\n{h['title']}\ndoi:{h['doi']} pmid:{h['pmid']}\n"
                 f"{'='*78}\n{ab}\n")
        body = "\n".join(ab.splitlines()[3:])
        print(body[:1100].strip())

print("\nwrote 110d_abstracts.txt")
print("\nStill NOT done: no candidate scored into S34 §9.1, no locus-gene "
      "content opened.")
