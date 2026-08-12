"""Step 91 -- systematic audit of eQTL-instrumented MR target-nomination studies.

WHY, AND WHAT THIS FIXES
Supplementary S16 compares this paper against three studies, chosen because we
had them to hand, and only against their main text. Reviewer 2 objected that
three papers cannot show whether the problems we describe are general. This
script replaces that with a search-defined sample.

THE TRAP WE ARE AVOIDING
Concluding "the field does not do check X" from the absence of X in abstracts
would be exactly the absence-of-evidence error Reviewer 1 caught us making about
naevi and the five cancers. Abstracts are ~250 words; nobody reports locus
annotation there. So:

  * abstracts are used ONLY to decide eligibility;
  * every practice claim is scored on FULL TEXT, and only for studies whose full
    text we could actually obtain (PMC open access);
  * the denominator reported is always the full-text set, never the search set;
  * "not found" is reported as not found in the text we could read, and the
    regexes are printed so a reader can see what would have counted.

SCORING CRITERIA -- fixed here before any full text was read.
  C1 known-locus attribution: does the paper compare its significant signal
     against previously reported loci for its own outcome trait?
  C2 colocalisation: performed at all?
  C3 outcome-power / list-stability: does it report how its candidate list
     depends on the outcome GWAS used, or on its power?
  C4 single-instrument disclosure: does it state how many instruments per
     exposure, or acknowledge that a one-SNP MR P value is the outcome P value?
  C5 SMR/HEIDI used as the arbiter without colocalisation?

Output: 91a_search.tsv, 91b_fulltext_scores.tsv, 91c_summary.tsv
"""
import io
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

MR = r"D:/R_ex/MR"
SCRATCH = os.path.join(MR, "litaudit")
os.makedirs(SCRATCH, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ---------------------------------------------------------------- criteria
CRITERIA = {
    "C1_known_locus": [
        r"known (?:risk )?loc[iu]s", r"previously (?:reported|identified|known) loc",
        r"established (?:risk )?loc", r"reported GWAS loc", r"known susceptibility loc",
        r"annotat\w+ (?:against|to) (?:known|previously)",
    ],
    "C2_coloc": [
        r"\bcoloc", r"colocali[sz]", r"\bSuSiE\b", r"\beCAVIAR\b", r"\bmoloc\b",
        r"posterior probability of (?:a )?shared", r"\bPP\.?H4\b", r"\bH4\b",
    ],
    "C3_power_stability": [
        r"(?:power|sample size) of the outcome", r"outcome GWAS power",
        r"different outcome (?:GWAS|dataset)", r"sensitivity to the outcome",
        r"down-?sampl", r"replicat\w+ (?:the )?(?:candidate )?list",
        r"stability of the (?:candidate )?list",
    ],
    "C4_instrument_count": [
        r"single (?:genetic )?(?:instrument|variant|SNP)", r"one instrument",
        r"number of instruments", r"instruments? per (?:gene|exposure|probe)",
        r"only one SNP", r"a single cis-",
    ],
    "C5_smr_heidi": [r"\bSMR\b", r"\bHEIDI\b"],
}
ELIGIBLE = [r"mendelian randomi", r"\beQTL\b|expression quantitative trait"]
NOMINATION = [r"drug target", r"therapeutic target", r"target identification",
              r"causal gene", r"candidate gene", r"nominat", r"prioriti[sz]"]


def get(url, tries=3):
    for k in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=90).read()
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(2)


def efetch_abstracts(pmids):
    out = {}
    for i in range(0, len(pmids), 150):
        chunk = pmids[i:i + 150]
        url = (f"{EUTILS}/efetch.fcgi?db=pubmed&retmode=xml&id=" + ",".join(chunk))
        xml = get(url).decode("utf-8", "replace")
        for art in xml.split("<PubmedArticle>")[1:]:
            pm = re.search(r"<PMID[^>]*>(\d+)</PMID>", art)
            ti = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", art, re.S)
            ab = " ".join(re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", art, re.S))
            jr = re.search(r"<Title>(.*?)</Title>", art, re.S)
            yr = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", art, re.S)
            pmc = re.search(r'<ArticleId IdType="pmc">(PMC\d+)</ArticleId>', art)
            if pm:
                clean = lambda t: re.sub(r"<[^>]+>", "", t).strip() if t else ""
                out[pm.group(1)] = dict(
                    pmid=pm.group(1), title=clean(ti.group(1) if ti else ""),
                    abstract=clean(ab), journal=clean(jr.group(1) if jr else ""),
                    year=yr.group(1) if yr else "", pmc=pmc.group(1) if pmc else "")
        time.sleep(0.4)
        print(f"  abstracts {min(i+150, len(pmids))}/{len(pmids)}", end="\r")
    print()
    return out


def fetch_pmc(pmcid):
    path = os.path.join(SCRATCH, f"{pmcid}.txt")
    if os.path.exists(path):
        return io.open(path, encoding="utf-8").read()
    url = f"{EUTILS}/efetch.fcgi?db=pmc&retmode=xml&id={pmcid}"
    try:
        xml = get(url).decode("utf-8", "replace")
    except Exception:
        return ""
    body = re.sub(r"<[^>]+>", " ", xml)
    body = re.sub(r"\s+", " ", body)
    io.open(path, "w", encoding="utf-8").write(body)
    time.sleep(0.4)
    return body


def score(text):
    t = text.lower()
    return {k: any(re.search(p, t) for p in pats) for k, pats in CRITERIA.items()}


def main():
    ids = json.load(io.open(os.path.join(MR, "..", "..", "litaudit_ids.json"),
                            encoding="utf-8")) if False else None
    src = sys.argv[1]
    ids = json.load(io.open(src, encoding="utf-8"))["union"]
    print(f"search set: {len(ids)} PMIDs")
    meta = efetch_abstracts(ids)
    print(f"fetched metadata for {len(meta)}")

    rows = []
    for pm, m in meta.items():
        blob = (m["title"] + " " + m["abstract"]).lower()
        elig = all(any(re.search(p, blob) for p in [pat]) for pat in ELIGIBLE)
        nom = any(re.search(p, blob) for p in NOMINATION)
        rows.append(dict(**m, eligible=bool(elig and nom)))
    import csv
    with io.open(f"{MR}/91a_search.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            r = dict(r); r["abstract"] = r["abstract"][:400]
            w.writerow(r)
    elig = [r for r in rows if r["eligible"]]
    withpmc = [r for r in elig if r["pmc"]]
    print(f"eligible (eQTL + MR + nomination language): {len(elig)}")
    print(f"  of which with a PMC id: {len(withpmc)}")

    scored = []
    for i, r in enumerate(withpmc):
        txt = fetch_pmc(r["pmc"])
        if len(txt) < 5000:
            continue
        sc = score(txt)
        scored.append(dict(pmid=r["pmid"], pmc=r["pmc"], year=r["year"],
                           journal=r["journal"], title=r["title"][:120],
                           chars=len(txt), **{k: int(v) for k, v in sc.items()}))
        print(f"  full text {i+1}/{len(withpmc)}  {r['pmc']}  {len(txt):,} chars", end="\r")
    print()
    with io.open(f"{MR}/91b_fulltext_scores.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(scored[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(scored)

    n = len(scored)
    summ = []
    for k in CRITERIA:
        c = sum(r[k] for r in scored)
        summ.append(dict(criterion=k, n_with=c, n_total=n, pct=round(100 * c / n, 1)))
    with io.open(f"{MR}/91c_summary.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["criterion", "n_with", "n_total", "pct"], delimiter="\t")
        w.writeheader()
        w.writerows(summ)
    print(f"\n=== full-text scored: {n} studies ===")
    for s in summ:
        print(f"  {s['criterion']:<22} {s['n_with']:3d}/{s['n_total']}  {s['pct']:5.1f}%")
    print("\nwrote 91a-91c")


if __name__ == "__main__":
    main()
