"""Step 91b -- rescore the literature audit after a bug, and under a second,
deliberately broader pattern set.

TWO CHANGES, KEPT SEPARATE ON PURPOSE

(1) BUG. step91 lower-cased each full text and then matched patterns that
    contained upper-case literals (\\bSMR\\b, \\bSuSiE\\b, \\bPP\\.?H4\\b,
    \\beCAVIAR\\b). Those could never match, which is why C5 came back 0/152 --
    an impossible number for this literature. Fixed by matching
    case-insensitively against the original text. This is a bug fix: it was
    wrong before and is right now.

(2) BROADENING. Having seen the first pass, we widened C1 and C3 to catch
    phrasings the pre-fixed patterns would have missed ("previously associated
    loci", "already known", "replication in an independent outcome", ...).
    This happened AFTER seeing results, so it is post-hoc and is reported as a
    separate column rather than replacing the strict count. Where the two
    disagree we report the WIDER (more generous to the literature) number in
    the manuscript, because the claim being made is that these checks are rare;
    a generous denominator makes that claim harder, not easier, to sustain.

Output: 91d_rescored.tsv, 91e_summary_both.tsv
"""
import csv
import glob
import io
import os
import re

MR = r"D:/R_ex/MR"
SCRATCH = os.path.join(MR, "litaudit")

STRICT = {
    "C1_known_locus": [
        r"known (?:risk )?loc[iu]s", r"previously (?:reported|identified|known) loc",
        r"established (?:risk )?loc", r"reported GWAS loc", r"known susceptibility loc",
        r"annotat\w+ (?:against|to) (?:known|previously)",
    ],
    "C2_coloc": [
        r"\bcoloc", r"colocali[sz]", r"\bSuSiE\b", r"\beCAVIAR\b", r"\bmoloc\b",
        r"posterior probability of (?:a )?shared", r"\bPP\.?H4\b",
    ],
    "C3_power_stability": [
        r"(?:power|sample size) of the outcome", r"outcome GWAS power",
        r"different outcome (?:GWAS|dataset)", r"sensitivity to the outcome",
        r"down-?sampl", r"stability of the (?:candidate )?list",
    ],
    "C4_instrument_count": [
        r"single (?:genetic )?(?:instrument|variant|SNP)", r"one instrument",
        r"number of instruments", r"instruments? per (?:gene|exposure|probe)",
        r"only one SNP", r"a single cis-",
    ],
    "C5_smr_heidi": [r"\bSMR\b", r"\bHEIDI\b"],
}

WIDE = {
    "C1_known_locus": STRICT["C1_known_locus"] + [
        r"previously associated loc", r"already (?:known|reported|identified)",
        r"novel loc", r"known gene[s]? for", r"GWAS[- ]significant loc",
        r"reported (?:in|by) (?:previous|prior) GWAS", r"catalog(?:ue)? of known",
        r"distance to (?:the )?(?:nearest )?known",
    ],
    "C2_coloc": STRICT["C2_coloc"] + [r"\bH4\b", r"shared causal variant",
                                      r"\bHyPrColoc\b", r"\bCOLOC\b"],
    "C3_power_stability": STRICT["C3_power_stability"] + [
        r"alternative outcome", r"second(?:ary)? outcome GWAS", r"replication (?:GWAS|cohort) for the outcome",
        r"statistical power", r"power calculation", r"power analysis",
        r"varying the outcome", r"another (?:melanoma |cancer )?GWAS",
    ],
    "C4_instrument_count": STRICT["C4_instrument_count"] + [
        r"number of (?:genetic )?variants? (?:used )?as instruments",
        r"one SNP per", r"lead (?:cis-)?eQTL", r"top (?:cis-)?eQTL",
    ],
    "C5_smr_heidi": STRICT["C5_smr_heidi"] + [r"summary[- ]data[- ]based mendelian"],
}


def score(text, pats):
    return {k: int(any(re.search(p, text, re.I) for p in v)) for k, v in pats.items()}


def main():
    prev = {r["pmc"]: r for r in csv.DictReader(
        io.open(f"{MR}/91b_fulltext_scores.tsv", encoding="utf-8"), delimiter="\t")}
    rows = []
    for pmc, meta in prev.items():
        path = os.path.join(SCRATCH, f"{pmc}.txt")
        if not os.path.exists(path):
            continue
        t = io.open(path, encoding="utf-8").read()
        s_strict = score(t, STRICT)
        s_wide = score(t, WIDE)
        rows.append(dict(pmid=meta["pmid"], pmc=pmc, year=meta["year"],
                         journal=meta["journal"], title=meta["title"],
                         **{f"strict_{k}": v for k, v in s_strict.items()},
                         **{f"wide_{k}": v for k, v in s_wide.items()}))
    with io.open(f"{MR}/91d_rescored.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    summ = []
    for k in STRICT:
        cs = sum(r[f"strict_{k}"] for r in rows)
        cw = sum(r[f"wide_{k}"] for r in rows)
        summ.append(dict(criterion=k, n=n, strict_n=cs, strict_pct=round(100*cs/n, 1),
                         wide_n=cw, wide_pct=round(100*cw/n, 1)))
    with io.open(f"{MR}/91e_summary_both.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(summ)

    print(f"=== {n} full texts, case-insensitive matching ===")
    print(f"{'criterion':<22}{'strict':>14}{'wide (post-hoc)':>20}")
    for s in summ:
        print(f"  {s['criterion']:<20} {s['strict_n']:3d}/{n} {s['strict_pct']:5.1f}%   "
              f"{s['wide_n']:3d}/{n} {s['wide_pct']:5.1f}%")
    # journals and years, for the descriptive table
    yrs = {}
    for r in rows:
        yrs[r["year"]] = yrs.get(r["year"], 0) + 1
    print("\nby year:", dict(sorted(yrs.items())))
    print("\nwrote 91d, 91e")


if __name__ == "__main__":
    main()
