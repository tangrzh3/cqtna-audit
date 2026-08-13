"""Step 94c -- the disease x eQTL-resource grid.

Executes PREREG_generality_grid.md. Every cell asks the same question with the
same machinery: of the FDR<0.05 MR records under this (exposure, outcome) pair,
what share of independent loci carries a known lead SNP FOR THAT DISEASE?

Two things are different from every earlier cross-disease analysis in this paper
and both are pre-registered:
  * each disease is scored against ITS OWN known-locus list, not melanoma's.
    Using melanoma's list on other cancers (step 30) tests specificity; this
    tests attribution, and the two are not interchangeable.
  * a mismatched negative control is included: melanoma's known-locus list
    applied to pancreatic cancer. If that also enriches, the whole grid is an
    artefact of locus density and is void (prereg §5, §6-E).

Known-locus lists arrive as rsIDs; positions come from the outcome file itself,
so no liftover is involved and every disease is scored in its own build.

Output: 94d_grid.tsv, 94e_grid_summary.tsv
"""
import csv
import gzip
import io
import os
from math import lgamma, exp

import numpy as np
import pandas as pd
from scipy.stats import norm

MR = r"D:/R_ex/MR"
LOCUS_KB = 1000
KNOWN_KB = 1000


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    hi = min(r1, c1)
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, hi + 1)), 1.0)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


def assign_loci(df):
    out = {}
    for ch, sub in df.groupby("chr"):
        sub = sub.sort_values("pos")
        lid, prev = 0, None
        for _, r in sub.iterrows():
            if prev is not None and r.pos - prev > LOCUS_KB * 1000:
                lid += 1
            out[(ch, r.pos)] = f"{ch}_{lid}"
            prev = r.pos
    return out


def known_positions(rsids, outcome_df):
    """map known-locus rsIDs to positions using the outcome file's own rsid column"""
    hits = outcome_df[outcome_df.rsid.isin(rsids)]
    return {c: np.sort(s.pos.values) for c, s in hits.groupby("chr")}


def attribute(d, kn):
    def is_known(ch, pos):
        arr = kn.get(str(ch))
        if arr is None or not len(arr):
            return False
        i = np.searchsorted(arr, pos)
        return any(0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= KNOWN_KB * 1000
                   for j in (i - 1, i))
    snp = d.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp)
    d = d.copy()
    d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
    d["known"] = [is_known(c, p) for c, p in zip(d.chr, d.pos)]
    bg = d.groupby("locus").agg(known=("known", "any"))
    BG_T, BG_K = len(bg), int(bg.known.sum())
    sig = d[d.fdr < 0.05]
    sg = sig.groupby("locus").agg(known=("known", "any"))
    S_T, S_K = len(sg), int(sg.known.sum())
    fold = ((S_K / S_T) / (BG_K / BG_T)) if S_T and BG_K else np.nan
    p = fisher_greater(S_K, S_T - S_K, BG_K - S_K,
                       (BG_T - BG_K) - (S_T - S_K)) if S_T else np.nan
    return dict(bg_loci=BG_T, bg_known=BG_K, sig_records=len(sig),
                sig_loci=S_T, sig_known=S_K,
                pct_known=round(100 * S_K / S_T, 1) if S_T else np.nan,
                fold=round(fold, 2) if S_T else np.nan,
                fisher_p=p, nominal_hits=int((d.p_mr < .05).sum()))


def load_known(name):
    """GRCh38-resolved list (step94b2). Falls back to None so the cell is deferred
    rather than scored against a list we could not place on the genome."""
    path = f"{MR}/known_loci_{name}_grch38.csv"
    if not os.path.exists(path) or os.path.getsize(path) < 200:
        return None
    d = pd.read_csv(path)
    return {c: np.sort(s.pos.values) for c, s in d.astype({"chr": str}).groupby("chr")}


def mr_from_extract(inst, ext, label):
    """inst: exposure instruments with chr,pos,beta_exp; ext: FinnGen extract"""
    ext = ext.rename(columns={"#chrom": "chr"})
    ext["chr"] = ext["chr"].astype(str)
    ext["pos"] = ext["pos"].astype(int)
    m = inst.merge(ext[["chr", "pos", "rsids", "ref", "alt", "beta", "sebeta", "af_alt"]],
                   on=["chr", "pos"], how="inner", suffixes=("", "_out"))
    m = m[(m.sebeta > 0) & (m.beta_exp != 0)]
    z = m.beta / m.sebeta
    out = pd.DataFrame(dict(chr=m.chr, pos=m.pos, gene_id=m.gene_id,
                            rsid=m.rsids.fillna(""), p_mr=2 * norm.sf(np.abs(z))))
    out["fdr"] = bh(out.p_mr.values)
    print(f"    {label}: {len(out):,} MR records, FDR<0.05 {int((out.fdr<.05).sum())}")
    return out


def main():
    rows = []

    # ---------------- column 1: Soskic (already computed, read back) ----------
    mel_sos = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
    mel_sos[["chr", "pos"]] = mel_sos.SNP.str.split(":", expand=True)
    mel_sos["pos"] = mel_sos.pos.astype(int)
    mel_sos = mel_sos.rename(columns={"FDR": "fdr", "pval": "p_mr"})
    # melanoma known list is already encoded in the 'category' column
    sig = mel_sos[mel_sos.fdr < 0.05]
    snp = mel_sos.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp)
    mel_sos["locus"] = [loci[(c, p)] for c, p in zip(mel_sos.chr, mel_sos.pos)]
    mel_sos["known"] = mel_sos.category != "潜在新位点"
    bg = mel_sos.groupby("locus").agg(known=("known", "any"))
    sg = mel_sos[mel_sos.fdr < 0.05].groupby("locus").agg(known=("known", "any"))
    rows.append(dict(exposure="Soskic_CD4", disease="melanoma",
                     bg_loci=len(bg), bg_known=int(bg.known.sum()),
                     sig_records=len(sig), sig_loci=len(sg),
                     sig_known=int(sg.known.sum()),
                     pct_known=round(100 * sg.known.sum() / len(sg), 1),
                     fold=round((sg.known.sum() / len(sg)) /
                                (bg.known.sum() / len(bg)), 2),
                     fisher_p=fisher_greater(int(sg.known.sum()),
                                             len(sg) - int(sg.known.sum()),
                                             int(bg.known.sum()) - int(sg.known.sum()),
                                             (len(bg) - int(bg.known.sum())) -
                                             (len(sg) - int(sg.known.sum()))),
                     nominal_hits=int((mel_sos.p_mr < .05).sum()), source="13"))

    for tag, fn in (("HCC_high", "85a_HCC_high_annotated.tsv"),
                    ("HCC_low", "85a_HCC_low_annotated.tsv")):
        d = pd.read_csv(f"{MR}/{fn}", sep="\t")
        d["chr"] = d["chr"].astype(str)
        bg = d.groupby("locus").agg(known=("known", "any"))
        sig = d[d.fdr < 0.05]
        sg = sig.groupby("locus").agg(known=("known", "any"))
        rows.append(dict(exposure="Soskic_CD4", disease=tag,
                         bg_loci=len(bg), bg_known=int(bg.known.sum()),
                         sig_records=len(sig), sig_loci=len(sg),
                         sig_known=int(sg.known.sum()),
                         pct_known=round(100 * sg.known.sum() / len(sg), 1),
                         fold=round((sg.known.sum() / len(sg)) /
                                    (bg.known.sum() / len(bg)), 2),
                         fisher_p=fisher_greater(int(sg.known.sum()),
                                                 len(sg) - int(sg.known.sum()),
                                                 int(bg.known.sum()) - int(sg.known.sum()),
                                                 (len(bg) - int(bg.known.sum())) -
                                                 (len(sg) - int(sg.known.sum()))),
                         nominal_hits=int((d.p_mr < .05).sum()), source="85a"))

    # ---------------- column 2: eQTLGen ---------------------------------------
    eq_mel = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
    eq_mel["chr"] = eq_mel["chr"].astype(str)
    bg = eq_mel.groupby("locus").agg(known=("known", "any"))
    sig = eq_mel[eq_mel.fdr < 0.05]
    sg = sig.groupby("locus").agg(known=("known", "any"))
    rows.append(dict(exposure="eQTLGen_blood", disease="melanoma",
                     bg_loci=len(bg), bg_known=int(bg.known.sum()),
                     sig_records=len(sig), sig_loci=len(sg),
                     sig_known=int(sg.known.sum()),
                     pct_known=round(100 * sg.known.sum() / len(sg), 1),
                     fold=round((sg.known.sum() / len(sg)) /
                                (bg.known.sum() / len(bg)), 2),
                     fisher_p=fisher_greater(int(sg.known.sum()),
                                             len(sg) - int(sg.known.sum()),
                                             int(bg.known.sum()) - int(sg.known.sum()),
                                             (len(bg) - int(bg.known.sum())) -
                                             (len(sg) - int(sg.known.sum()))),
                     nominal_hits=int((eq_mel.p_mr < .05).sum()), source="92c"))

    inst = eq_mel[["chr", "pos", "gene_id", "beta_exp"]].drop_duplicates()
    CAN = {"lung": "C3_BRONCHUS_LUNG_EXALLC", "colorectal": "C3_COLORECTAL_EXALLC",
           "pancreas": "C3_PANCREAS_EXALLC", "breast": "C3_BREAST_EXALLC",
           "prostate": "C3_PROSTATE_EXALLC"}
    mel_known = set(pd.read_csv(f"{MR}/landi2020_known_loci_grch38.csv").rsid) \
        if "rsid" in pd.read_csv(f"{MR}/landi2020_known_loci_grch38.csv", nrows=1).columns else set()

    for name, endpoint in CAN.items():
        path = f"{MR}/cancer_extracts_eqtlgen/{endpoint}.tsv"
        if not os.path.exists(path):
            print(f"  [missing extract] {name}")
            continue
        ext = pd.read_csv(path, sep="\t")
        d = mr_from_extract(inst, ext, f"eQTLGen x {name}")
        kn_ids = load_known(name)
        if kn_ids is None:
            print(f"  [no known-locus list] {name} -- cell deferred")
            continue
        ext2 = ext.rename(columns={"#chrom": "chr"})
        ext2["chr"] = ext2["chr"].astype(str)
        ext2["rsid"] = ext2["rsids"].fillna("").str.split(",").str[0]
        res = attribute(d, kn_ids)
        rows.append(dict(exposure="eQTLGen_blood", disease=name,
                         n_known_list=sum(len(v) for v in kn_ids.values()),
                         source="94a", **res))
        # ---- negative control: melanoma's list on pancreatic cancer
        if name == "pancreas" and mel_known:
            knm = known_positions(mel_known, ext2[["chr", "pos", "rsid"]])
            nc = attribute(d, knm)
            rows.append(dict(exposure="eQTLGen_blood",
                             disease="pancreas_NC_melanoma_list",
                             n_known_list=len(mel_known), source="NC", **nc))

    # ---------------- eQTLGen x HCC -------------------------------------------
    print("  (eQTLGen x HCC cells require the HCC outcome scan; see step94d)")

    g = pd.DataFrame(rows)
    g.to_csv(f"{MR}/94d_grid.tsv", sep="\t", index=False)
    print("\n=== grid ===")
    cols = ["exposure", "disease", "sig_loci", "sig_known", "pct_known", "fold",
            "fisher_p", "bg_known", "bg_loci"]
    print(g[cols].to_string(index=False))

    testable = g[(g.sig_loci >= 2) & (~g.disease.str.contains("_NC_"))]
    frac = (testable.fold > 1).mean() if len(testable) else np.nan
    print(f"\nU1: cells with >=2 significant loci: {len(testable)}; "
          f"fold>1 in {(testable.fold>1).sum()}/{len(testable)} = {100*frac:.0f}% "
          f"(pre-registered threshold 80%)")
    pd.DataFrame([dict(n_testable=len(testable),
                       n_fold_gt1=int((testable.fold > 1).sum()),
                       pct=round(100 * frac, 1) if len(testable) else np.nan)]).to_csv(
        f"{MR}/94e_grid_summary.tsv", sep="\t", index=False)
    print("\nwrote 94d, 94e")


if __name__ == "__main__":
    main()
