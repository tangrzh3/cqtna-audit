#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 158 -- execute S53: is instrumentability associated with readability?

Executes manuscript/PREREG_visibility_overlap.md (S53), committed 2026-09-06 in
its own commit before any crosstab was formed. S53 fixes the universe, the four
variables, the primary endpoint, the calibration and the stopping rules; this
script implements them and adds nothing.

ORDER MATTERS AND IS ENFORCED HERE
S53 section 4 makes the within-stratum permutation a kill criterion, not a
footnote. If the permuted CMH statistic does not centre on zero the
stratification is wrong, and section 3's numbers are not interpreted. That test
runs before the result is printed as a result.

THE PRIMARY ENDPOINT IS THE CONDITIONAL ONE
Both variables are functions of expression in opposite directions -- eQTL
detection power rises with expression, probe panels drop the top of the
distribution -- so the marginal odds ratio is expected to be negative for a
reason that is not about the platforms being independently biased. S53 section
5 says so in advance, and says that a null conditional association is the
finding rather than a failure.

  python step158_visibility_overlap.py <dir-with-downloaded-tables> [repo]

Outputs: 158a_visibility_2x2.tsv, 158b_visibility_by_decile.tsv, 158c_console.log
"""
import glob
import gzip
import io
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

DATA = sys.argv[1] if len(sys.argv) > 1 else "."
MR = (sys.argv[2] if len(sys.argv) > 2
      else os.environ.get("CQTNA_DIR") or os.path.dirname(os.path.abspath(__file__)))

N_PERM = 1000
SEED = 158
MIN_CELL = 10          # S53 section 6
MIN_STRATA = 5         # S53 section 6
CALIB_TOL = 0.05       # S53 section 4

LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def cmh(tabs):
    """Cochran-Mantel-Haenszel odds ratio and test over 2x2 strata."""
    num = den = 0.0
    e = v = obs = 0.0
    for a, b, c, d in tabs:
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
        r1, r2 = a + b, c + d
        c1, c2 = a + c, b + d
        obs += a
        e += r1 * c1 / n
        if n > 1:
            v += r1 * r2 * c1 * c2 / (n * n * (n - 1))
    orr = num / den if den > 0 else np.nan
    chi = (abs(obs - e) - 0.5) ** 2 / v if v > 0 else np.nan
    p = stats.chi2.sf(chi, 1) if np.isfinite(chi) else np.nan
    return orr, chi, p


def strata_tables(df, col):
    out = []
    for _, g in df.groupby("decile"):
        a = int(((g.instrumentable == 1) & (g[col] == 1)).sum())
        b = int(((g.instrumentable == 1) & (g[col] == 0)).sum())
        c = int(((g.instrumentable == 0) & (g[col] == 1)).sum())
        d = int(((g.instrumentable == 0) & (g[col] == 0)).sum())
        out.append((a, b, c, d))
    return out


def main():
    # ------------------------------------------------- expression (S53 sect. 2)
    files = sorted(glob.glob(os.path.join(DATA, "bulk", "Donor*.csv.gz")))
    if not files:
        say("no bulk RNA-seq files -> cannot build the conditioning variable.")
        return 2
    frames = []
    for f in files:
        with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:
            b = pd.read_csv(fh)
        cpm = [c for c in b.columns if c.endswith("_cpm")]
        frames.append(pd.DataFrame({"gene_id": b.gene_id.astype(str).str.split(".").str[0],
                                    "cpm": b[cpm].mean(axis=1)}))
    expr = (pd.concat(frames).groupby("gene_id").cpm.mean().reset_index())
    say("expression: %d donor files, %d genes with a CPM"
        % (len(files), len(expr)))

    # ------------------------------------------------------- universe + labels
    kd = pd.read_csv(os.path.join(DATA, "kd.csv"))
    kd["gene_id"] = kd.perturbed_gene_id.astype(str).str.split(".").str[0]
    universe = pd.DataFrame({"gene_id": sorted(kd.gene_id.unique())})
    say("universe: %d genes in the perturbation library" % len(universe))

    with gzip.open(os.path.join(DATA, "downstream.csv.gz"), "rt",
                   encoding="utf-8", errors="replace") as fh:
        dn = pd.read_csv(fh)
    readable_ids = set(dn.downstream_gene_ids.astype(str).str.split(".").str[0])
    say("readable set: %d gene ids appear as measurable downstream genes"
        % len(readable_ids))

    inst = pd.read_csv(os.path.join(MR, "CD4_dynamic_top_eqtl_instruments_all.csv"))
    inst_ids = set(inst.gene_id.astype(str).str.split(".").str[0])
    say("instrumented set: %d genes with a genome-wide significant top cis-eQTL"
        % len(inst_ids))

    vk = kd.groupby("gene_id").signif_knockdown.max()
    df = universe.merge(expr, on="gene_id", how="inner")
    df["instrumentable"] = df.gene_id.isin(inst_ids).astype(int)
    df["readable"] = df.gene_id.isin(readable_ids).astype(int)
    df["verifiable_kd"] = df.gene_id.map(vk).fillna(False).astype(int)
    say("analysed: %d genes with both a CPM and library membership" % len(df))

    # ------------------------------------------------ marginal (S53 sect. 3.1)
    say()
    say("=" * 74)
    say("1. marginal, shown for contrast only -- NOT the endpoint")
    say("=" * 74)
    rows = []
    for col in ("readable", "verifiable_kd"):
        a = int(((df.instrumentable == 1) & (df[col] == 1)).sum())
        b = int(((df.instrumentable == 1) & (df[col] == 0)).sum())
        c = int(((df.instrumentable == 0) & (df[col] == 1)).sum())
        d = int(((df.instrumentable == 0) & (df[col] == 0)).sum())
        say("  %s: instrumentable&yes %d, instrumentable&no %d,"
            " not-instrumentable&yes %d, not&no %d" % (col, a, b, c, d))
        if min(a, b, c, d) < MIN_CELL:
            say("    a cell has < %d genes -> no odds ratio (S53 section 6)"
                % MIN_CELL)
            rows.append(dict(variable=col, a=a, b=b, c=c, d=d,
                             odds_ratio=np.nan, p=np.nan))
            continue
        orr, p = stats.fisher_exact([[a, b], [c, d]])
        say("    OR %.3f   Fisher P %.3g" % (orr, p))
        rows.append(dict(variable=col, a=a, b=b, c=c, d=d, odds_ratio=orr, p=p))
    pd.DataFrame(rows).to_csv(os.path.join(MR, "158a_visibility_2x2.tsv"),
                              sep="\t", index=False)

    # ---------------------------------------------- conditioning (sect. 3.2/4)
    df["decile"] = pd.qcut(df.cpm.rank(method="first"), 10, labels=False)
    usable = sum(1 for _, g in df.groupby("decile")
                 if g.readable.nunique() > 1 and g.instrumentable.nunique() > 1)
    say()
    say("=" * 74)
    say("2. calibration BEFORE the conditional result (S53 section 4)")
    say("=" * 74)
    say("  expression deciles carrying both labels: %d of 10" % usable)
    if usable < MIN_STRATA:
        say("  fewer than %d -> marginal only, no conditional claim." % MIN_STRATA)
        return 0

    rng = np.random.default_rng(SEED)
    obs_or, _, _ = cmh(strata_tables(df, "readable"))
    null = []
    for _ in range(N_PERM):
        p = df.copy()
        p["readable"] = (p.groupby("decile").readable
                          .transform(lambda s: rng.permutation(s.values)))
        o, _, _ = cmh(strata_tables(p, "readable"))
        null.append(np.log(o) if o and np.isfinite(o) and o > 0 else np.nan)
    null = np.array([x for x in null if np.isfinite(x)])
    med = float(np.median(null))
    say("  permuted log OR: median %+.4f, sd %.4f, n = %d"
        % (med, float(np.std(null)), len(null)))
    if abs(med) > CALIB_TOL:
        say("  |median| > %.2f -> CALIBRATION FAILED. The stratification is"
            % CALIB_TOL)
        say("  wrong, and nothing in section 3 is interpreted (S53 section 4).")
        return 1
    say("  calibration passes: the permuted statistic centres on zero.")

    say()
    say("=" * 74)
    say("3. PRIMARY ENDPOINT -- conditional on expression (S53 section 3.2)")
    say("=" * 74)
    out = []
    for col in ("readable", "verifiable_kd"):
        tabs = strata_tables(df, col)
        orr, chi, p = cmh(tabs)
        emp = (np.mean(np.abs(null) >= abs(np.log(orr)))
               if col == "readable" and orr and orr > 0 else np.nan)
        say("  %-14s CMH OR %.3f   chi2 %.2f   P %.3g%s"
            % (col, orr, chi, p,
               "   permutation P %.3f" % emp if np.isfinite(emp) else ""))
        for i, (a, b, c, d) in enumerate(tabs):
            out.append(dict(variable=col, decile=i, a=a, b=b, c=c, d=d))
    pd.DataFrame(out).to_csv(os.path.join(MR, "158b_visibility_by_decile.tsv"),
                             sep="\t", index=False)

    say()
    say("Read against S53 section 5, written before this ran: a negative")
    say("marginal association that largely disappears under conditioning IS")
    say("the predicted result, and a null conditional association is a finding")
    say("and not a failure -- it would mean the two platforms are not")
    say("independently biased but reading the same axis. S53 section 7 caps")
    say("what this may be used for.")
    return 0


if __name__ == "__main__":
    rc = main()
    io.open(os.path.join(MR, "158c_console.log"), "w",
            encoding="utf-8", newline="\n").write("\n".join(LOG) + "\n")
    raise SystemExit(rc)
