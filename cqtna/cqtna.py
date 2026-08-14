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
    """locus attribution, by independent locus"""
    f = known_flagger(known, cfg["known_window_kb"])
    d = d.copy()
    d["known"] = [f(c, p) for c, p in zip(d.chr, d.pos)]
    bg = d.groupby("locus").agg(known=("known", "any"))
    BT, BK = len(bg), int(bg.known.sum())
    sig = d[d.fdr < cfg["fdr_threshold"]]
    sg = sig.groupby("locus").agg(known=("known", "any"))
    ST, SK = len(sg), int(sg.known.sum())
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
    """list stability against a second outcome GWAS"""
    t = cfg["fdr_threshold"]
    a = set(d.gene[d.fdr < t])
    b = set(d2.gene[d2.fdr < t])
    la = set(d.locus[d.fdr < t])
    lb = set(d2.locus[d2.fdr < t])
    return dict(genes_outcome1=len(a), genes_outcome2=len(b),
                genes_shared=len(a & b),
                jaccard_gene=round(len(a & b) / len(a | b), 3) if (a | b) else None,
                loci_outcome1=len(la), loci_outcome2=len(lb),
                jaccard_locus=round(len(la & lb) / len(la | lb), 3) if (la | lb) else None,
                lost=sorted(a - b), gained=sorted(b - a))


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
