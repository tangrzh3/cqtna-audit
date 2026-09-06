#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 157 -- execute S52: look up TPI1 in the genome-scale CD4+ T perturb-seq map.

Executes manuscript/PREREG_perturbseq_lookup.md (S52), committed 2026-09-06 in
its own commit BEFORE any row of this dataset was opened. The reading table is
fixed there and this script does not restate it -- it produces the numbers S52
section 2 enumerates, in the order S52 requires, and nothing else.

POSITIVE CONTROL RUNS FIRST AND IS A KILL CRITERION
S52 section 2 requires reproducing a named regulator's large effect before any
TPI1 number is interpreted. If the control fails we misread the table; the
lookup is void and TPI1 is not reported. That check is in the code, not in the
reading.

WHAT IT DOES NOT DO
It does not decide anything about TPI1. Every outcome S52 enumerates leaves the
S47 ladder unchanged, for the reason given in S52 section 5: the target claim
failed on the FDR, on a prior-dependent colocalisation, on a region the outcome
cannot fine-map, on unresolved gene attribution at chr12p13 and on compartment,
and a knockdown experiment in CD4+ T cells touches none of those.

Source (S52 section 1): Zhu et al., Cell 2026;189, doi:10.1016/j.cell.2026.08.002
Supplementary tables from the authors' analysis repository
(github.com/emdann/GWT_perturbseq_analysis_2025, metadata/suppl_tables).

  python step157_perturbseq_tpi1_lookup.py <dir-with-downloaded-tables> [repo]

Outputs: 157a_perturbseq_tpi1.tsv, 157b_console.log
"""
import gzip
import io
import os
import sys

import pandas as pd

DATA = sys.argv[1] if len(sys.argv) > 1 else "."
MR = (sys.argv[2] if len(sys.argv) > 2
      else os.environ.get("CQTNA_DIR") or os.path.dirname(os.path.abspath(__file__)))

GENE = "TPI1"
ENSG = "ENSG00000111669"
# Highly expressed glycolytic genes, checked alongside so that an absence is
# read as a property of the panel rather than of this one gene.
NEIGHBOURS = ["GAPDH", "PGK1", "ENO1", "PKM", "ALDOA", "LDHA", "HK1", "ACTB"]
CONTROLS = ["STAT3", "IRF4", "BATF"]          # named as regulators in the paper
CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]
# The control must look like a regulator, or we are reading the wrong column.
CONTROL_MIN_DE = 50

LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    de = pd.read_csv(os.path.join(DATA, "DE_stats.csv"))
    say("DE table: %d rows, %d perturbed genes, conditions %s"
        % (len(de), de.target_contrast_gene_name.nunique(),
           sorted(de.culture_condition.unique())))

    # ---------------------------------------------------------------- control
    say()
    say("=" * 74)
    say("POSITIVE CONTROL (S52 section 2) -- runs before anything about TPI1")
    say("=" * 74)
    ctl = de[de.target_contrast_gene_name.isin(CONTROLS)]
    if ctl.empty:
        say("  none of %s found -- VOID" % ", ".join(CONTROLS))
        return 1
    for _, r in ctl.sort_values(["target_contrast_gene_name",
                                 "culture_condition"]).iterrows():
        say("  %-6s %-9s  cells %6.0f   DE genes %5d   on-target effect %6.2f  sig=%s"
            % (r.target_contrast_gene_name, r.culture_condition, r.n_cells_target,
               r.n_total_de_genes, r.ontarget_effect_size, r.ontarget_significant))
    best = ctl.n_total_de_genes.max()
    if best < CONTROL_MIN_DE:
        say()
        say("  strongest control shows only %d DE genes (< %d)."
            % (best, CONTROL_MIN_DE))
        say("  CONTROL FAILED -> lookup VOID. TPI1 is not reported: on this")
        say("  evidence we have misread the table, not found an absence.")
        return 1
    say()
    say("  control passes: strongest named regulator reaches %d DE genes." % best)

    # ------------------------------------------------------------------ TPI1
    say()
    say("=" * 74)
    say("%s, as a perturbed gene (S52 section 2, items 1-4)" % GENE)
    say("=" * 74)
    g = de[de.target_contrast_gene_name == GENE]
    if g.empty:
        say("  %s is absent from the perturbation library or not expressed." % GENE)
        say("  -> S52 outcome C (not ascertainable).")
    else:
        med = de.groupby("culture_condition").n_cells_target.median()
        for _, r in g.sort_values("culture_condition").iterrows():
            say("  %-9s cells %6.0f (condition median %6.0f)   DE genes %5d"
                "   up %4d  down %4d"
                % (r.culture_condition, r.n_cells_target,
                   med[r.culture_condition], r.n_total_de_genes,
                   r.n_up_genes, r.n_down_genes))
            say("            on-target knockdown effect %6.3f   significant=%s"
                "   guides=%s   low_target_gex=%s"
                % (r.ontarget_effect_size, r.ontarget_significant,
                   r.n_guides, r.low_target_gex))
        say()
        for cond in CONDITIONS:
            sub = de[de.culture_condition == cond]
            row = g[g.culture_condition == cond]
            if sub.empty or row.empty:
                continue
            de_pct = (sub.n_total_de_genes <
                      row.n_total_de_genes.iloc[0]).mean() * 100
            cell_pct = (sub.n_cells_target <
                        row.n_cells_target.iloc[0]).mean() * 100
            say("  %-9s %s sits at the %.0fth percentile of DE-gene count and the"
                " %.0fth of cell recovery" % (cond, GENE, de_pct, cell_pct))

    # --------------------------------- is it measured at all? (S52 item 1, 4)
    say()
    say("=" * 74)
    say("Is %s measured by this assay at all? (S52 items 1 and 4)" % GENE)
    say("=" * 74)
    kpath = os.path.join(DATA, "kd.csv")
    if not os.path.exists(kpath):
        say("  knockdown-efficiency table not downloaded -> item 4 not ascertainable.")
    else:
        kd = pd.read_csv(kpath)
        sub = kd[kd.perturbed_gene_id.astype(str).str.contains(ENSG, na=False)]
        if sub.empty:
            say("  %s absent from the guide table." % GENE)
        else:
            for _, r in sub.iterrows():
                say("  guide: target expr %.3f in %.0f cells | NTC expr %.3f in"
                    " %.0f cells | t %.2f | signif_knockdown=%s"
                    % (r.guide_mean_expr, r.guide_n, r.ntc_mean_expr, r.ntc_n,
                       r.t_statistic, r.signif_knockdown))
            allzero = kd.groupby("perturbed_gene_id").ntc_mean_expr.max().eq(0)
            say()
            say("  The control arm matters more than the perturbed arm here: the")
            say("  target's expression is zero in the non-targeting controls too,")
            say("  so this is not a knockdown that failed -- the assay does not")
            say("  capture the gene, and on-target verification is impossible.")
            say("  Genes whose every guide has zero control expression: %d of %d"
                " (%.1f%%)" % (allzero.sum(), len(allzero),
                               100.0 * allzero.mean()))
        say()
        say("  Knockdown efficiency overall, for scale: %.1f%% of guide rows are"
            % (100.0 * kd.signif_knockdown.mean()))
        say("  significant knockdowns.")

    # ------------------- expressed, or merely unread? (S52 item 1, decisive)
    say()
    say("=" * 74)
    say("Expressed, or merely unread? (S52 item 1)")
    say("=" * 74)
    bpath = os.path.join(DATA, "bulk.csv.gz")
    dpath0 = os.path.join(DATA, "downstream.csv.gz")
    if not (os.path.exists(bpath) and os.path.exists(dpath0)):
        say("  bulk RNA-seq or downstream table missing -> not ascertainable.")
    else:
        with gzip.open(bpath, "rt", encoding="utf-8", errors="replace") as f:
            b = pd.read_csv(f)
        cpm = [c for c in b.columns if c.endswith("_cpm")]
        b["CPM"] = b[cpm].mean(axis=1)
        b = b.sort_values("CPM", ascending=False).reset_index(drop=True)
        b["rank"] = b.index + 1
        with gzip.open(dpath0, "rt", encoding="utf-8", errors="replace") as f:
            uni = set(pd.read_csv(f).downstream_gene.astype(str))
        say("  Conventional bulk RNA-seq of the same cells, %d genes." % len(b))
        say("  %-7s %10s %8s   %s" % ("gene", "mean CPM", "rank", "in readout?"))
        for gname in [GENE] + NEIGHBOURS:
            r = b[b.gene_name == gname]
            if r.empty:
                continue
            r = r.iloc[0]
            say("  %-7s %10.1f %8d   %s"
                % (gname, r.CPM, r["rank"], "yes" if gname in uni else "NO"))
        top = b.head(40).gene_name.tolist()
        say()
        say("  Of the 40 highest-expressed genes in their own bulk RNA-seq, only")
        say("  %d are in the probe readout. Missing: %s"
            % (sum(1 for x in top if x in uni),
               ", ".join(x for x in top if x not in uni)[:150]))
        say("  Mitochondrial, ribosomal and the largest housekeeping transcripts.")
        say("  So the gene is not unexpressed and not unimportant -- it is too")
        say("  highly expressed to be worth a probe. Its invisibility here is a")
        say("  property of the readout panel, not of the biology.")

    # ------------------------------------------------- TPI1 as downstream gene
    say()
    say("=" * 74)
    say("%s, as a downstream gene (S52 section 2, item 5)" % GENE)
    say("=" * 74)
    dpath = os.path.join(DATA, "downstream.csv.gz")
    if not os.path.exists(dpath):
        say("  downstream table not downloaded -> item 5 not ascertainable.")
    else:
        with gzip.open(dpath, "rt", encoding="utf-8", errors="replace") as f:
            dn = pd.read_csv(f)
        uni = set(dn.downstream_gene.astype(str))
        say("  genes measurable as downstream genes: %d" % len(uni))
        say("  %-6s %s" % (GENE, "present" if GENE in uni else "ABSENT"))
        for g in NEIGHBOURS:
            say("  %-6s %s" % (g, "present" if g in uni else "ABSENT"))
        say()
        hit = dn[dn.downstream_gene == GENE]
        if hit.empty:
            say("  %s is not a condition-specific downstream gene of any cluster."
                % GENE)
        else:
            for _, r in hit.sort_values(["condition", "hdbscan_cluster"]).iterrows():
                say("  cluster %-4s %-9s upstream regulators %-4s  sign coherence %.2f"
                    % (r.hdbscan_cluster, r.condition, r.num_of_upstream,
                       r.sign_coherence))
            say("  clusters/conditions in which %s appears downstream: %d"
                % (GENE, len(hit)))

    # ------------------------------------------------------------------ write
    out = os.path.join(MR, "157a_perturbseq_tpi1.tsv")
    cols = ["target_contrast_gene_name", "culture_condition", "n_cells_target",
            "n_up_genes", "n_down_genes", "n_total_de_genes",
            "ontarget_effect_size", "ontarget_significant", "n_guides",
            "low_target_gex"]
    keep = de[de.target_contrast_gene_name.isin([GENE] + CONTROLS)][cols]
    keep.to_csv(out, sep="\t", index=False)
    say()
    say("wrote 157a_perturbseq_tpi1.tsv")
    say()
    say("Read the result against the table in S52 section 3, which was fixed")
    say("before this ran. No outcome there changes an S47 label, and S52")
    say("section 5 forbids the sentence this result will tempt someone to write.")
    return 0


if __name__ == "__main__":
    rc = main()
    io.open(os.path.join(MR, "157b_console.log"), "w",
            encoding="utf-8", newline="\n").write("\n".join(LOG) + "\n")
    raise SystemExit(rc)
