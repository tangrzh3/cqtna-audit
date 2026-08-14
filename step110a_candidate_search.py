#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 110a -- S34 §3.2 step 1: find candidate truth sets and verify one reference.

S34 §3.2 requires candidates to be scored against the six admissibility criteria
BEFORE any candidate's locus-to-gene pairs are opened. This script therefore
collects METADATA ONLY -- titles, years, DOIs, journals -- and downloads no
locus-gene table. Scoring is a judgement recorded by hand in S34 §9.1.

It also resolves the one missing bibliographic entry flagged in REFERENCES.md:
BioGRID ORCS is cited in the text ("628 of 1,471 human CRISPR screens") but has
no verified entry, so the claim currently carries no reference number.

The eQTL Catalogue API is returning HTTP 500 on every endpoint in this session,
so the C1 resource survey cannot be refreshed here. NCBI E-utilities and CrossRef
are reachable.

Outputs: 110a_benchmark_candidates.tsv
         110b_orcs_reference.tsv
         110a_console.log
"""
import io
import json
import sys
import time
import urllib.parse
import urllib.request

MR = r"D:/R_ex/MR"
UA = {"User-Agent": "Mozilla/5.0 (attribution-benchmark-survey)"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CROSSREF = "https://api.crossref.org/works"


class Tee:
    def __init__(self, path):
        self.f = io.open(path, "w", encoding="utf-8")

    def write(self, s):
        sys.__stdout__.write(s)
        self.f.write(s)

    def flush(self):
        sys.__stdout__.flush()
        self.f.flush()


sys.stdout = Tee(f"{MR}/110a_console.log")


def fetch(url, tries=4):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.load(urllib.request.urlopen(req, timeout=60))
        except Exception as ex:
            if a == tries - 1:
                print(f"    !! failed: {type(ex).__name__} "
                      f"{getattr(ex, 'code', '')}")
                return None
            time.sleep(2 * (a + 1))
    return None


def pubmed(query, n=8):
    """esearch + esummary. Returns list of dicts."""
    q = urllib.parse.urlencode(
        dict(db="pubmed", term=query, retmax=n, retmode="json", sort="relevance"))
    res = fetch(f"{EUTILS}/esearch.fcgi?{q}")
    if not res:
        return []
    ids = res.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    q2 = urllib.parse.urlencode(
        dict(db="pubmed", id=",".join(ids), retmode="json"))
    summ = fetch(f"{EUTILS}/esummary.fcgi?{q2}")
    if not summ:
        return []
    out = []
    for pmid in ids:
        r = summ.get("result", {}).get(pmid)
        if not r:
            continue
        doi = ""
        for aid in r.get("articleids", []):
            if aid.get("idtype") == "doi":
                doi = aid.get("value", "")
        out.append(dict(pmid=pmid, title=r.get("title", ""),
                        journal=r.get("fulljournalname", ""),
                        year=(r.get("pubdate", "") or "")[:4], doi=doi))
        time.sleep(0.15)
    return out


# --------------------------------------------------------------------------
# 1. candidate truth sets -- metadata only, no locus-gene content
# --------------------------------------------------------------------------
QUERIES = [
    ("L2G / Open Targets",
     "locus-to-gene OR (locus AND gene AND prioritisation) AND GWAS AND "
     "(gold standard OR benchmark)"),
    ("gold-standard gene sets",
     "GWAS gold standard causal gene set benchmark locus assignment"),
    ("ProGeM-type frameworks",
     "framework prioritising causal genes GWAS loci gold standard evaluation"),
    ("experimentally anchored",
     "GWAS locus causal gene experimentally validated benchmark set human"),
    ("drug-target anchored",
     "drug target genetic support GWAS locus causal gene gold standard"),
]

print("=" * 78)
print("Step 110a -- S34 §3.2: candidate truth sets (METADATA ONLY)")
print("=" * 78)
print("No locus-to-gene table is downloaded here. Admissibility scoring against")
print("S34 §3.1 is a hand judgement recorded in §9.1 of that document.\n")

rows, seen = [], set()
for label, q in QUERIES:
    print(f"--- {label}")
    print(f"    query: {q}")
    hits = pubmed(q)
    if not hits:
        print("    (no results)")
    for h in hits:
        key = h["pmid"]
        if key in seen:
            continue
        seen.add(key)
        h["query_label"] = label
        rows.append(h)
        print(f"    {h['year']:<6}{h['journal'][:30]:<32}{h['title'][:70]}")
    print()

if rows:
    import pandas as pd
    df = pd.DataFrame(rows)[["query_label", "year", "journal", "title",
                             "doi", "pmid"]]
    df.to_csv(f"{MR}/110a_benchmark_candidates.tsv", sep="\t", index=False)
    print(f"wrote 110a_benchmark_candidates.tsv  ({len(df)} unique records)")

# --------------------------------------------------------------------------
# 2. BioGRID ORCS -- the missing reference
# --------------------------------------------------------------------------
print("\n" + "=" * 78)
print("BioGRID ORCS: resolving the reference cited but not listed")
print("=" * 78)

orcs = []
for q in ["BioGRID ORCS CRISPR screen database",
          "BioGRID Open Repository CRISPR Screens"]:
    print(f"--- pubmed: {q}")
    for h in pubmed(q, n=6):
        print(f"    {h['year']:<6}{h['journal'][:30]:<32}{h['title'][:70]}")
        print(f"           doi:{h['doi']}  pmid:{h['pmid']}")
        orcs.append(h)

# CrossRef cross-check, filtering the traps REFERENCES.md documents
print("\n--- crossref bibliographic query (filtering preprints / recommendations)")
q = urllib.parse.urlencode({"query.bibliographic":
                            "BioGRID ORCS open repository CRISPR screens",
                            "rows": 8})
cr = fetch(f"{CROSSREF}?{q}")
if cr:
    for it in cr.get("message", {}).get("items", []):
        typ = it.get("type", "")
        cont = " ".join(it.get("container-title", []) or [])
        if typ == "posted-content" or "Faculty Opinions" in cont:
            continue
        title = (it.get("title") or [""])[0]
        yr = ""
        for k in ("published-print", "published-online", "issued"):
            if it.get(k, {}).get("date-parts"):
                yr = str(it[k]["date-parts"][0][0])
                break
        print(f"    {yr:<6}{cont[:30]:<32}{title[:70]}")
        print(f"           doi:{it.get('DOI','')}  type:{typ}")

if orcs:
    import pandas as pd
    pd.DataFrame(orcs).to_csv(f"{MR}/110b_orcs_reference.tsv",
                              sep="\t", index=False)
    print(f"\nwrote 110b_orcs_reference.tsv")

print("\n" + "=" * 78)
print("NOT done here: no candidate has been scored against S34 §3.1, and no")
print("locus-to-gene content has been opened. Scoring is the next step and must")
print("be written into S34 §9.1 before any hit rate is computed.")
print("=" * 78)
