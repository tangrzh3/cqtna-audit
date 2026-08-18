#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""CQTNA -- Context-specific QTL Target Nomination Audit.

Runs the automatable part of the eight diagnostics on a candidate list produced
by any cis-eQTL-instrumented MR pipeline, and emits an evidence tier per gene.

    python cqtna.py --config my_audit.yaml
    python cqtna.py --demo                 # runs on the audit's own data

What it does NOT do is as important as what it does. Three of the eight
diagnostics cannot be automated from a candidate list, and the report says so by
name rather than omitting them -- see CHECKLIST in cqtna_manual.md. A run that
reports five green modules is not an audited nomination; it is five of eight.

Modules
  A  locus attribution        fold enrichment of significant loci on the
                              outcome's own known loci, one-sided Fisher, with a
                              mismatched-list negative control when supplied
  B  inference unit           the FDR<t list recomputed over record, variant,
                              gene and independent locus
  C  list stability           the same pipeline against a second outcome GWAS:
                              Jaccard, and the known/novel split of what moves
  D  instrument attrition     instrumentable -> analysable -> associated
  E  compartment attribution  per-cell-type expression ratio for nominated genes
  F  peak distance            base pairs between eQTL and GWAS peaks

Everything is computed by independent locus where a count carries an argument,
because a significant locus does not name a gene.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from math import exp, lgamma

import numpy as np
import pandas as pd

VERSION = "0.1"
DEFAULTS = dict(fdr_threshold=0.05, locus_window_kb=1000, known_window_kb=1000,
                target_cell_type=None, top_n=None)


# ----------------------------------------------------------------- statistics
def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


def simes(p):
    p = np.sort(np.asarray(p, float))
    n = len(p)
    return float(np.min(p * n / np.arange(1, n + 1)))


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, min(r1, c1) + 1)), 1.0)


def md_table(df):
    """markdown table without pulling in tabulate -- this tool is meant to be
    run by other people, so its dependencies stay at numpy and pandas"""
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |",
           "|" + "|".join("---" for _ in cols) + "|"]
    for r in df.itertuples(index=False):
        out.append("| " + " | ".join("" if v is None or (isinstance(v, float)
                                                          and np.isnan(v))
                                     else str(v) for v in r) + " |")
    return "\n".join(out)


def assign_loci(chrs, poss, window_kb):
    order = sorted(range(len(chrs)), key=lambda i: (str(chrs[i]), int(poss[i])))
    out = [None] * len(chrs)
    lc, lp, lid = None, None, 0
    for i in order:
        c, p = str(chrs[i]), int(poss[i])
        if c != lc or p - lp > window_kb * 1000:
            lid += 1
        out[i] = lid
        lc, lp = c, p
    return out


def known_flagger(known, window_kb):
    KN = {c: np.sort(s.pos.values)
          for c, s in known.astype({"chr": str}).groupby("chr")}

    def f(ch, pos):
        arr = KN.get(str(ch))
        if arr is None or not len(arr):
            return False
        i = np.searchsorted(arr, pos)
        return any(0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= window_kb * 1000
                   for j in (i - 1, i))
    return f


def nearest_known_bp(known):
    """distance in bp to the nearest known lead SNP; inf if that chromosome
    carries none. Same neighbour search as known_flagger, so the two can never
    disagree about whether a variant is inside a given window."""
    KN = {c: np.sort(s.pos.values)
          for c, s in known.astype({"chr": str}).groupby("chr")}

    def f(ch, pos):
        arr = KN.get(str(ch))
        if arr is None or not len(arr):
            return float("inf")
        i = np.searchsorted(arr, pos)
        return min((abs(int(arr[j]) - int(pos))
                    for j in (i - 1, i) if 0 <= j < len(arr)), default=float("inf"))
    return f


# ----------------------------------------------------------------- input
REQUIRED_MR = ["record_id", "gene", "chr", "pos", "p"]
REQUIRED_KNOWN = ["chr", "pos"]


def load_mr(path, cfg):
    d = pd.read_csv(path, sep=None, engine="python")
    missing = [c for c in REQUIRED_MR if c not in d.columns]
    if missing:
        sys.exit(f"{path}: missing required column(s) {missing}. "
                 f"Required: {REQUIRED_MR}")
    d = d[d.p.notna()].copy()
    d["chr"] = d["chr"].astype(str)
    d["pos"] = d["pos"].astype(int)
    d["locus"] = assign_loci(d.chr.tolist(), d.pos.tolist(), cfg["locus_window_kb"])
    d["fdr"] = bh(d.p.values)
    return d


def load_known(path):
    k = pd.read_csv(path, sep=None, engine="python")
    missing = [c for c in REQUIRED_KNOWN if c not in k.columns]
    if missing:
        sys.exit(f"{path}: missing required column(s) {missing}")
    k["chr"] = k["chr"].astype(str)
    k["pos"] = k["pos"].astype(int)
    return k


# ----------------------------------------------------------------- modules
def module_a(d, known, cfg, label="known-locus list"):
    """locus attribution, by independent locus

    ⚠ "known" is a property of the LOCUS, computed once over every record and
    inherited by the significant loci and by the gene labels. The first version
    took the denominator over all records, the numerator over significant
    records only, and labelled genes from per-record flags -- three conventions
    that can disagree, so a gene could be listed as novel while its locus was
    counted as known. They happen to agree on the demo data, which is why this
    needed a test rather than an inspection.
    """
    f = known_flagger(known, cfg["known_window_kb"])
    d = d.copy()
    rec = pd.Series([f(c, p) for c, p in zip(d.chr, d.pos)], index=d.index)
    locus_known = rec.groupby(d.locus).any()          # 唯一口径
    d["known"] = d.locus.map(locus_known)             # 广播回记录

    BT, BK = len(locus_known), int(locus_known.sum())
    sig = d[d.fdr < cfg["fdr_threshold"]]
    sig_loci = sig.locus.unique()
    ST, SK = len(sig_loci), int(locus_known[sig_loci].sum())
    fold = ((SK / ST) / (BK / BT)) if ST and BK else float("nan")
    p = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)) if ST else float("nan")
    return dict(reference=label, background_known=BK, background_loci=BT,
                background_pct=round(100 * BK / BT, 2) if BT else None,
                significant_loci=ST, significant_known=SK,
                pct_known=round(100 * SK / ST, 1) if ST else None,
                fold=round(fold, 2) if ST and BK else None,
                fisher_p_one_sided=p,
                known_genes=sorted(set(sig.gene[sig.known])),
                novel_genes=sorted(set(sig.gene[~sig.known])))


def module_b(d, cfg):
    """the list under four inference units"""
    t = cfg["fdr_threshold"]
    out = []
    rec_genes = set(d.gene[d.fdr < t])
    out.append(dict(unit="record", n_tests=len(d),
                    n_significant=int((d.fdr < t).sum()),
                    n_genes=len(rec_genes), n_loci=d.locus[d.fdr < t].nunique()))
    for unit, key in (("variant", ["chr", "pos"]), ("gene", ["gene"]),
                      ("locus", ["locus"])):
        g = d.groupby(key).p.apply(lambda s: simes(s.values)).reset_index(name="p")
        g["fdr"] = bh(g.p.values)
        hit = g[g.fdr < t]
        sub = d.merge(hit[key], on=key, how="inner")
        out.append(dict(unit=unit, n_tests=len(g), n_significant=len(hit),
                        n_genes=sub.gene.nunique(), n_loci=sub.locus.nunique()))
    tab = pd.DataFrame(out)
    return tab


def module_c(d, d2, cfg):
    """list stability against a second outcome GWAS

    ⚠ Locus identity must be recomputed on the UNION of both tables' coordinates.
    Each table's own `locus` column is an integer counter assigned within that
    table, so the counters are not comparable: on the demo data both tables carry
    an id 206, one at 16:87.7 Mb and one at 16:89.7 Mb, and a naive set
    intersection called them the same locus. Clustering the union gives both
    tables one shared partition, and the key carries the chromosome so two
    chromosomes can never match.
    """
    t = cfg["fdr_threshold"]
    n1 = len(d)
    chrs = d.chr.tolist() + d2.chr.tolist()
    poss = d.pos.tolist() + d2.pos.tolist()
    uni = assign_loci(chrs, poss, cfg["locus_window_kb"])
    # 位点键写成 chr:start-end，既可比又能看出簇的跨度
    ext = {}
    for c, p, g in zip(chrs, poss, uni):
        lo, hi = ext.get(g, (p, p))
        ext[g] = (min(lo, p), max(hi, p))
    key = [f"{c}:{ext[g][0]}-{ext[g][1]}" for c, g in zip(chrs, uni)]
    k1, k2 = key[:n1], key[n1:]

    s1 = (d.fdr < t).values
    s2 = (d2.fdr < t).values
    a = set(d.gene[s1])
    b = set(d2.gene[s2])
    la = {k for k, keep in zip(k1, s1) if keep}
    lb = {k for k, keep in zip(k2, s2) if keep}
    v1 = set(zip(d.chr, d.pos, d.gene))
    v2 = set(zip(d2.chr, d2.pos, d2.gene))
    return dict(genes_outcome1=len(a), genes_outcome2=len(b),
                genes_shared=len(a & b),
                jaccard_gene=round(len(a & b) / len(a | b), 3) if (a | b) else None,
                loci_outcome1=len(la), loci_outcome2=len(lb),
                loci_shared=len(la & lb),
                jaccard_locus=round(len(la & lb) / len(la | lb), 3) if (la | lb) else None,
                lost=sorted(a - b), gained=sorted(b - a),
                shared_locus_keys=sorted(la & lb),
                coverage=dict(records_outcome1=n1, records_outcome2=len(d2),
                              gene_variant_pairs_shared=len(v1 & v2),
                              gene_variant_pairs_only1=len(v1 - v2),
                              gene_variant_pairs_only2=len(v2 - v1),
                              genes_never_tested_in_2=sorted(a - set(d2.gene)),
                              genes_never_tested_in_1=sorted(b - set(d.gene))))


def module_d(path):
    """instrument attrition at three levels"""
    t = pd.read_csv(path, sep=None, engine="python")
    for c in ("instrumentable", "analysable", "associated"):
        if c not in t.columns:
            sys.exit(f"{path}: missing column {c}")
    n = len(t)
    return dict(pathway_genes=n,
                instrumentable=int(t.instrumentable.sum()),
                analysable=int(t.analysable.sum()),
                associated=int(t.associated.sum()),
                note="report all three; collapsing them attributes an "
                     "outcome-side limit to the biology")


def module_e(path, genes, cfg):
    """compartment attribution: per-cell-type ratio for nominated genes"""
    x = pd.read_csv(path, sep=None, engine="python")
    for c in ("gene", "cell_type", "mean_expression"):
        if c not in x.columns:
            sys.exit(f"{path}: missing column {c}")
    tgt = cfg.get("target_cell_type")
    rows = []
    for g in genes:
        s = x[x.gene == g]
        if s.empty or tgt is None or tgt not in set(s.cell_type):
            continue
        v = s.set_index("cell_type").mean_expression
        ref = v[tgt]
        top = v.drop(index=tgt).idxmax()
        rows.append(dict(gene=g, target_cell_type=tgt,
                         target_expression=round(float(ref), 4),
                         highest_other=top,
                         highest_other_expression=round(float(v[top]), 4),
                         ratio_other_over_target=round(float(v[top] / ref), 2)
                         if ref else None))
    return pd.DataFrame(rows)


SWEEP_KB = (100, 250, 500, 1000)


def module_g(d, known, cfg):
    """window sensitivity, and the continuous distance behind the binary flag

    Two different conventions get called "the 1 Mb window" and they carry
    different weight, so they are swept separately:

      known_window_kb  what counts as landing on a known locus -- this one moves
                       the fold directly
      locus_window_kb  how independent loci are defined -- this one moves the
                       denominator

    A threshold chosen to flatter a result weakens when tightened. Reporting the
    sweep is what distinguishes a convention from a tuned parameter.

    The distance for a SIGNIFICANT locus is taken over that locus's significant
    records only, because that is what the binary flag is computed from. Taking
    it over every record in the locus lets a non-significant variant supply the
    distance and produces a number that contradicts the flag.
    """
    dist = nearest_known_bp(known)
    d = d.copy()
    d["dist_bp"] = [dist(c, p) for c, p in zip(d.chr, d.pos)]
    t = cfg["fdr_threshold"]

    def enrich(loci, known_kb):
        kn = d.dist_bp <= known_kb * 1000
        tmp = pd.DataFrame(dict(locus=loci, known=kn.values, fdr=d.fdr.values))
        bg = tmp.groupby("locus").agg(known=("known", "any"))
        sg = tmp[tmp.fdr < t].groupby("locus").agg(known=("known", "any"))
        BT, BK, ST, SK = len(bg), int(bg.known.sum()), len(sg), int(sg.known.sum())
        fold = ((SK / ST) / (BK / BT)) if ST and BK else None
        p = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)) if ST else None
        return dict(background_loci=BT, background_known=BK,
                    significant_loci=ST, significant_known=SK,
                    fold=round(fold, 2) if fold is not None else None,
                    fisher_p_one_sided=p)

    known_sweep = [dict(known_window_kb=kb,
                        locus_window_kb=cfg["locus_window_kb"],
                        **enrich(d.locus.values, kb)) for kb in SWEEP_KB]

    locus_sweep = []
    for kb in SWEEP_KB:
        loci = assign_loci(d.chr.tolist(), d.pos.tolist(), kb)
        locus_sweep.append(dict(locus_window_kb=kb,
                                known_window_kb=cfg["known_window_kb"],
                                **enrich(loci, cfg["known_window_kb"])))

    # ⚠ 距离必须与二分类同源。位点的 known 状态由**该位点全部记录**决定，
    # 所以主报的距离也取全部记录的最小值；另报只用显著记录算的那一列，
    # 两者不同时说明该位点是靠一条不显著的记录才贴近已知 lead SNP。
    sig = d[d.fdr < t]
    bg_all = (d.groupby("locus").dist_bp.min()
                .replace(float("inf"), np.nan).dropna())
    sig_keys = sig.locus.unique()
    sig_d = bg_all[bg_all.index.isin(sig_keys)].sort_values()
    sig_rec = (sig.groupby("locus").dist_bp.min()
                  .replace(float("inf"), np.nan).dropna().sort_values())
    bg_d = bg_all
    # 阈值落在数据的空隙里，还是正踩在数据上？后者说明结论依赖这个阈值
    gap = None
    if len(sig_d) > 1:
        v = sig_d.values
        below = v[v <= cfg["known_window_kb"] * 1000]
        above = v[v > cfg["known_window_kb"] * 1000]
        if len(below) and len(above):
            gap = dict(nearest_below_bp=int(below[-1]),
                       nearest_above_bp=int(above[0]),
                       threshold_bp=cfg["known_window_kb"] * 1000)
    return dict(known_sweep=known_sweep, locus_sweep=locus_sweep,
                significant_distances_bp=[int(v) for v in sig_d.values],
                significant_distances_sig_records_bp=[int(v) for v in sig_rec.values],
                significant_median_bp=int(sig_d.median()) if len(sig_d) else None,
                background_median_bp=int(bg_d.median()) if len(bg_d) else None,
                threshold_gap=gap)


def module_f(path):
    """distance between eQTL and GWAS peaks"""
    x = pd.read_csv(path, sep=None, engine="python")
    for c in ("gene", "eqtl_pos", "gwas_pos"):
        if c not in x.columns:
            sys.exit(f"{path}: missing column {c}")
    x["distance_bp"] = (x.eqtl_pos.astype(int) - x.gwas_pos.astype(int)).abs()
    return x[["gene", "eqtl_pos", "gwas_pos", "distance_bp"]]


# ----------------------------------------------------------------- tiers
def assign_tiers(a, c, e, cfg):
    """screened / unresolved / state-informative / target-supported

    The tier says what KIND of evidence a gene has, and stability across outcome
    GWAS is reported beside it rather than folded into it. Those are different
    facts: a gene on a known locus that also fails to replicate is still a gene
    on a known locus, and collapsing the two produced rows reading "unresolved"
    over a basis saying "already known".

    Deliberately conservative: nothing reaches target-supported from this tool
    alone, because the three diagnostics it cannot run are exactly the ones that
    would license that word.
    """
    rows = {}
    for g in a["known_genes"]:
        rows[g] = dict(tier="screened",
                       basis="significant, but on a locus already known for "
                             "this outcome",
                       stability="not tested")
    for g in a["novel_genes"]:
        rows[g] = dict(tier="unresolved",
                       basis="significant on a locus not previously reported "
                             "for this outcome",
                       stability="not tested")
    if c is not None:
        lost, gained = set(c["lost"]), set(c["gained"])
        for g, r in rows.items():
            r["stability"] = ("lost under the second outcome GWAS" if g in lost
                              else "retained under the second outcome GWAS")
        for g in gained:
            if g not in rows:
                rows[g] = dict(tier="unresolved",
                               basis="significant only under the second outcome "
                                     "GWAS",
                               stability="gained under the second outcome GWAS")
    if e is not None and len(e):
        for r in e.itertuples():
            if r.gene in rows and r.ratio_other_over_target and \
                    r.ratio_other_over_target > 2:
                rows[r.gene]["tier"] = "state-informative"
                rows[r.gene]["basis"] = (
                    f"expression dominated by {r.highest_other} "
                    f"({r.ratio_other_over_target}x the target cell type), so "
                    f"tissue-level data cannot validate a target-cell-specific "
                    f"mechanism")
    return rows


# ----------------------------------------------------------------- report
MANUAL = [
    ("(ii) colocalisation with explicit multiple-signal modelling and an LD "
     "reference matched to the outcome cohort",
     "needs regional summary statistics and an LD panel, not a candidate list"),
    ("(vii) cell-level matching on lineage composition when splitting cells by "
     "a score",
     "needs the single-cell data and the split itself; a cluster-level control "
     "can pass while the cell-level one fails"),
    ("(viii) code-by-code verification that an endpoint definition is unchanged "
     "before comparing candidate lists across releases",
     "needs the phenotype definitions, which release notes do not reliably "
     "summarise"),
]


def report(res, cfg, out_dir):
    L = []
    w = L.append
    w(f"# CQTNA report (v{VERSION})\n")
    w(f"Settings: FDR < {cfg['fdr_threshold']}, locus window "
      f"{cfg['locus_window_kb']} kb, known-locus window "
      f"{cfg['known_window_kb']} kb.\n")

    a = res["A"]
    w("## A. Locus attribution\n")
    w(f"- background: {a['background_known']}/{a['background_loci']} loci "
      f"({a['background_pct']}%) carry a known lead SNP for this outcome")
    w(f"- significant: {a['significant_known']}/{a['significant_loci']} "
      f"({a['pct_known']}%)")
    w(f"- **enrichment {a['fold']}x, one-sided Fisher P = "
      f"{a['fisher_p_one_sided']:.3g}**")
    if res.get("A_nc"):
        nc = res["A_nc"]
        w(f"- mismatched-list control: {nc['fold']}x, P = "
          f"{nc['fisher_p_one_sided']:.3g} "
          f"{'(clean)' if (nc['fold'] or 0) <= 1 or nc['fisher_p_one_sided'] > .05 else '**(POSITIVE -- attribution may be a density artefact)**'}")
    else:
        w("- ⚠ no mismatched-list control supplied; enrichment specific to this "
          "outcome's genetics is not established")
    w("")

    w("## B. Inference unit\n")
    w(md_table(res["B"]))
    w("\n⚠ A shortest list is not evidence that FDR is controlled under "
      "dependence. State which unit the conclusions are in.\n")

    if res.get("C"):
        c = res["C"]
        w("## C. List stability across outcome GWAS\n")
        w(f"- genes: {c['genes_outcome1']} vs {c['genes_outcome2']}, "
          f"{c['genes_shared']} shared, Jaccard {c['jaccard_gene']}")
        w(f"- loci: {c['loci_outcome1']} vs {c['loci_outcome2']}, "
          f"Jaccard {c['jaccard_locus']}")
        w(f"- lost: {', '.join(c['lost']) or 'none'}")
        w(f"- gained: {', '.join(c['gained']) or 'none'}\n")
    else:
        w("## C. List stability across outcome GWAS\n\n"
          "⚠ **not run** -- no second outcome supplied. This is the diagnostic "
          "that most often changes a conclusion.\n")

    if res.get("D"):
        d = res["D"]
        w("## D. Instrument attrition\n")
        w(f"- of {d['pathway_genes']} pathway genes: "
          f"**{d['instrumentable']} instrumentable, {d['analysable']} analysable "
          f"against this outcome, {d['associated']} nominally associated**\n")
    # A module that was configured but produced nothing is reported as such.
    # Silently dropping it would be the failure this tool exists to detect.
    for key, title, why in (
            ("E", "E. Compartment attribution",
             "no significant gene appears in the expression table, or no "
             "`target_cell_type` was set"),
            ("F", "F. eQTL-to-GWAS peak distance", "no peak table supplied")):
        w(f"## {title}\n")
        v = res.get(key)
        if v is not None and len(v):
            w(md_table(v))
            w("")
        elif v is not None:
            w(f"⚠ **configured but empty** — {why}.\n")
        else:
            w("⚠ **not run** — no input supplied.\n")

    if res.get("G"):
        g = res["G"]
        w("## G. Window sensitivity and the distances behind the flag\n")
        def fmt(rows, cols):
            t = pd.DataFrame(rows)[cols].copy()
            t["fisher_p_one_sided"] = [None if v is None else f"{v:.3g}"
                                       for v in t.fisher_p_one_sided]
            return md_table(t)

        w(fmt(g["known_sweep"],
              ["known_window_kb", "significant_known", "significant_loci",
               "background_known", "background_loci", "fold",
               "fisher_p_one_sided"]))
        w("")
        w(fmt(g["locus_sweep"],
              ["locus_window_kb", "background_loci", "significant_loci", "fold",
               "fisher_p_one_sided"]))
        w("")
        if g["significant_distances_bp"]:
            w("- distances of significant loci to the nearest known lead SNP (kb): "
              + ", ".join(f"{v/1000:,.0f}" for v in g["significant_distances_bp"]))
            w(f"- median {g['significant_median_bp']/1000:,.0f} kb, against "
              f"{g['background_median_bp']/1000:,.0f} kb for all testable loci")
        gap = g["threshold_gap"]
        if gap:
            w(f"- the {gap['threshold_bp']/1000:,.0f} kb threshold falls between "
              f"{gap['nearest_below_bp']/1000:,.0f} kb and "
              f"{gap['nearest_above_bp']/1000:,.0f} kb — "
              + ("**a gap, so the cut is not near your data**"
                 if gap["nearest_above_bp"] >= 4 * max(gap["nearest_below_bp"], 1)
                 else "**close to your data, so the result depends on this choice**"))
        w("\n⚠ A window chosen to flatter a result weakens when tightened. If the "
          "fold *rises* as the window narrows, the reported value is the "
          "conservative one.\n")

    w("## Evidence tiers\n")
    t = res["tiers"]
    if t:
        w("| gene | tier | basis | stability across outcome GWAS |")
        w("|---|---|---|---|")
        for g, r in sorted(t.items(), key=lambda kv: (kv[1]["tier"], kv[0])):
            w(f"| {g} | **{r['tier']}** | {r['basis']} | {r['stability']} |")
    else:
        w("no significant genes at this threshold")
    w("\n⚠ **No gene can reach `target-supported` from this tool.** The three "
      "diagnostics it cannot run are the ones that would license that word.\n")

    w("## Not run here -- these require manual work\n")
    for name, why in MANUAL:
        w(f"- **{name}**  \n  {why}")
    w("\nA nomination audited on five of eight diagnostics is audited on five "
      "of eight.\n")

    path = os.path.join(out_dir, "cqtna_report.md")
    open(path, "w", encoding="utf-8").write("\n".join(L))
    return path


# ----------------------------------------------------------------- driver
def run(cfg, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    d = load_mr(cfg["mr_results"], cfg)
    known = load_known(cfg["known_loci"])
    res = {"A": module_a(d, known, cfg)}
    if cfg.get("known_loci_mismatched"):
        res["A_nc"] = module_a(d, load_known(cfg["known_loci_mismatched"]), cfg,
                               label="mismatched (negative control)")
    res["B"] = module_b(d, cfg)
    res["C"] = module_c(d, load_mr(cfg["mr_results_alt"], cfg), cfg) \
        if cfg.get("mr_results_alt") else None
    res["D"] = module_d(cfg["instrument_table"]) if cfg.get("instrument_table") else None
    sig_genes = sorted(set(d.gene[d.fdr < cfg["fdr_threshold"]]))
    res["E"] = module_e(cfg["expression_table"], sig_genes, cfg) \
        if cfg.get("expression_table") else None
    res["F"] = module_f(cfg["peak_table"]) if cfg.get("peak_table") else None
    res["G"] = module_g(d, known, cfg)
    res["tiers"] = assign_tiers(res["A"], res["C"], res["E"], cfg)

    res["B"].to_csv(os.path.join(out_dir, "cqtna_units.tsv"), sep="\t", index=False)
    with open(os.path.join(out_dir, "cqtna_results.json"), "w", encoding="utf-8") as fh:
        json.dump({k: (v.to_dict("records") if isinstance(v, pd.DataFrame) else v)
                   for k, v in res.items() if k != "tiers"}, fh, indent=1,
                  default=str)
    p = report(res, cfg, out_dir)
    print(f"wrote {p}")
    return res


def main():
    ap = argparse.ArgumentParser(description="CQTNA " + VERSION)
    ap.add_argument("--config", help="JSON config; see cqtna_template/")
    ap.add_argument("--demo", action="store_true",
                    help="run on the audit's own melanoma data")
    ap.add_argument("--out", default="cqtna_out")
    a = ap.parse_args()
    if a.demo:
        here = os.path.dirname(os.path.abspath(__file__))
        cfg = dict(DEFAULTS)
        cfg.update(json.load(open(os.path.join(here, "demo", "demo_config.json"),
                                  encoding="utf-8")))
        for k in ("mr_results", "known_loci", "known_loci_mismatched",
                  "mr_results_alt", "instrument_table", "expression_table", "peak_table"):
            if cfg.get(k):
                cfg[k] = os.path.join(here, "demo", cfg[k])
        run(cfg, a.out)
    elif a.config:
        cfg = dict(DEFAULTS)
        cfg.update(json.load(open(a.config, encoding="utf-8")))
        # paths in a config are resolved relative to the config file, so a config
        # can be moved with its data and still run
        base = os.path.dirname(os.path.abspath(a.config))
        for k in ("mr_results", "mr_results_alt", "known_loci",
                  "known_loci_mismatched", "instrument_table",
                  "expression_table", "peak_table"):
            v = cfg.get(k)
            if v and not os.path.isabs(v):
                cfg[k] = os.path.join(base, v)
        run(cfg, a.out)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
