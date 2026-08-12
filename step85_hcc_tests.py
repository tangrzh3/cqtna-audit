"""Step 85 -- the three pre-registered tests for the HCC generalisation, plus the
matched-power control that keeps a null result interpretable.

Pre-registered in manuscript/PREREG_hcc_generalisation.md §4-§5. Nothing here is
chosen after seeing the HCC results: locus definition, counting unit (independent
loci, not records), enrichment test, direction, and the Neff matching formula were
all fixed in that document.

  G1  known-locus attribution of the FDR<0.05 signal, counted by independent locus
  G2  overlap of novel-locus gene nominations between HCC and melanoma
  G3  recovery rate known vs novel across the two HCC power levels (DESCRIPTIVE
      only -- the two HCC outcomes are not independent, see prereg §9)
  MP  melanoma down-sampled to HCC-high's effective size, 200 replicates

Output: 85a-85d
"""
import os
from math import lgamma, exp
import numpy as np
import pandas as pd
from scipy.stats import norm

MR = r"D:/R_ex/MR"
LOCUS_KB = 1000          # same 1 Mb single-linkage rule as Step 30e
KNOWN_KB = 1000          # instrument within 1 Mb of a known lead SNP -> known locus
NOVEL_MEL = "潜在新位点"
rng = np.random.default_rng(1)
N_REP = 200

# case/control counts fixed in the prereg
MEL_CASE, MEL_CTRL = 12530, 789099
HI_CASE, HI_CTRL = 3748, 1861536
LO_CASE, LO_CTRL = 947, 378749


def n_eff(nc, nk):
    return 4.0 / (1.0 / nc + 1.0 / nk)


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
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def assign_loci(df):
    """1 Mb single-linkage clustering of instrument positions within a chromosome."""
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


def main():
    known = pd.read_csv(f"{MR}/84a_hcc_known_loci_grch38.csv")
    known["chr"] = known["chr"].astype(str)
    kn_by_chr = {c: np.sort(s.pos.values) for c, s in known.groupby("chr")}
    print(f"known HCC loci with GRCh38 coordinates: {len(known)}")

    def is_known_snp(ch, pos):
        arr = kn_by_chr.get(str(ch))
        if arr is None or len(arr) == 0:
            return False
        i = np.searchsorted(arr, pos)
        for j in (i - 1, i):
            if 0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= KNOWN_KB * 1000:
                return True
        return False

    lines = []
    summary = []
    for tag, fn in (("HCC_high", "84b_mr_HCC_high_GCST90809296.tsv"),
                    ("HCC_low",  "84b_mr_HCC_low_FinnGenR12.tsv")):
        d = pd.read_csv(f"{MR}/{fn}", sep="\t")
        d["chr"] = d["chr"].astype(str)
        snp = d.groupby(["chr", "pos"], as_index=False).size()
        loci = assign_loci(snp)
        d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
        d["known"] = [is_known_snp(c, p) for c, p in zip(d.chr, d.pos)]

        # background: all instrument loci in the strict set
        bg = d.groupby("locus").agg(known=("known", "any"))
        BG_T, BG_K = len(bg), int(bg.known.sum())

        sig = d[d.fdr < 0.05]
        sg = sig.groupby("locus").agg(known=("known", "any"))
        S_T, S_K = len(sg), int(sg.known.sum())

        p = fisher_greater(S_K, S_T - S_K, BG_K - S_K, (BG_T - BG_K) - (S_T - S_K)) if S_T else float("nan")
        fold = ((S_K / S_T) / (BG_K / BG_T)) if S_T and BG_K else float("nan")

        print(f"\n[{tag}]  strict records {len(d)}, instrument loci {BG_T} "
              f"(known background {BG_K}/{BG_T} = {100*BG_K/BG_T:.1f}%)")
        print(f"  FDR<0.05: {len(sig)} records across {S_T} independent loci, "
              f"{S_K} known -> {fold:.2f}x, Fisher one-sided P = {p:.4g}")
        for L, s in sig.groupby("locus"):
            gs = sorted(set(s.gene_id))
            print(f"    {L:<8} known={bool(s.known.any())} "
                  f"pos={s.pos.min()}-{s.pos.max()} genes={','.join(gs)} "
                  f"n_rec={len(s)} minP={s.p_mr.min():.3g}")
            lines.append({"outcome": tag, "locus": L, "known": bool(s.known.any()),
                          "chr": s.chr.iloc[0], "pos_min": int(s.pos.min()),
                          "pos_max": int(s.pos.max()), "genes": ";".join(gs),
                          "n_records": len(s), "min_p": float(s.p_mr.min())})
        summary.append({"outcome": tag, "n_strict_records": len(d),
                        "bg_loci": BG_T, "bg_known": BG_K,
                        "sig_records": len(sig), "sig_loci": S_T, "sig_known": S_K,
                        "fold": fold, "fisher_p": p})
        d.to_csv(f"{MR}/85a_{tag}_annotated.tsv", sep="\t", index=False)

    pd.DataFrame(lines).to_csv(f"{MR}/85b_hcc_significant_loci.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(f"{MR}/85c_G1_enrichment.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- G2
    print("\n=== G2: novel-locus gene nominations, HCC vs melanoma ===")
    hi = pd.read_csv(f"{MR}/85a_HCC_high_annotated.tsv", sep="\t")
    lo = pd.read_csv(f"{MR}/85a_HCC_low_annotated.tsv", sep="\t")
    hcc_novel = set(hi[(hi.fdr < 0.05) & (~hi.known)].gene_id) | \
                set(lo[(lo.fdr < 0.05) & (~lo.known)].gene_id)

    mel = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
    mel_novel = set(mel[(mel.FDR < 0.05) & (mel.category == NOVEL_MEL)].gene_id)
    # the FinnGen round's four novel-locus genes (Step 11-12), by symbol
    mel_novel_symbols = set(mel[(mel.FDR < 0.05) & (mel.category == NOVEL_MEL)].SYMBOL.dropna())
    print(f"  melanoma meta novel-locus genes at FDR<0.05: {len(mel_novel)} {sorted(mel_novel_symbols)}")
    print(f"  HCC novel-locus genes at FDR<0.05: {len(hcc_novel)} {sorted(hcc_novel)}")
    print(f"  INTERSECTION: {sorted(hcc_novel & mel_novel)}  (n={len(hcc_novel & mel_novel)})")

    # ---------------------------------------------------------------- G3
    print("\n=== G3 (DESCRIPTIVE ONLY -- outcomes not independent) ===")
    lo_sig = lo[lo.fdr < 0.05]
    hi_sig = hi[hi.fdr < 0.05]
    for cls, mask_lo, mask_hi in (("known", lo_sig.known, hi_sig.known),
                                  ("novel", ~lo_sig.known, ~hi_sig.known)):
        g_lo = set(lo_sig[mask_lo].gene_id)
        g_hi = set(hi_sig[mask_hi].gene_id)
        rec = len(g_lo & g_hi) / len(g_lo) if g_lo else float("nan")
        print(f"  {cls:<6} low-power genes {len(g_lo)}, high-power genes {len(g_hi)}, "
              f"recovered {len(g_lo & g_hi)} -> {rec:.1%}" if g_lo else
              f"  {cls:<6} low-power genes 0 -- recovery undefined")

    # ---------------------------------------------------------------- MP
    print("\n=== MP: melanoma down-sampled to HCC-high effective size ===")
    ne_mel, ne_hi, ne_lo = n_eff(MEL_CASE, MEL_CTRL), n_eff(HI_CASE, HI_CTRL), n_eff(LO_CASE, LO_CTRL)
    print(f"  Neff  melanoma-meta {ne_mel:,.0f} | HCC-high {ne_hi:,.0f} "
          f"({ne_hi/ne_mel:.1%} of melanoma) | HCC-low {ne_lo:,.0f} ({ne_lo/ne_mel:.1%})")

    m = mel.dropna(subset=["beta_outcome", "se_outcome", "beta_exposure"]).copy()
    m = m[(m.se_outcome > 0) & (m.beta_exposure != 0)]
    m["is_novel"] = (m.category == NOVEL_MEL)
    m[["chr", "pos"]] = m.SNP.str.split(":", expand=True)
    m["pos"] = m.pos.astype(int)
    # count independent loci, exactly as G1 does, so the two are comparable
    mel_snp = m.groupby(["chr", "pos"], as_index=False).size()
    mel_loci = assign_loci(mel_snp)
    m["locus"] = [mel_loci[(c, p)] for c, p in zip(m.chr, m.pos)]
    print(f"  melanoma strict records used: {len(m)} across {m.locus.nunique()} instrument loci")

    rows = []
    for target, name in ((ne_hi, "HCC_high"), (ne_lo, "HCC_low")):
        infl = np.sqrt(ne_mel / target)
        se_sim = m.se_outcome.values * infl
        extra = np.sqrt(np.maximum(se_sim ** 2 - m.se_outcome.values ** 2, 0))
        nloc, nknown, nnovel = [], [], []
        for _ in range(N_REP):
            b = m.beta_outcome.values + rng.normal(0, extra)
            q = bh(2 * norm.sf(np.abs(b / se_sim)))
            sig = m[q < 0.05]
            nloc.append(sig.locus.nunique())
            nknown.append(sig[~sig.is_novel].locus.nunique())
            nnovel.append(sig[sig.is_novel].locus.nunique())
        rows.append({"matched_to": name, "target_neff": target,
                     "median_sig_loci": float(np.median(nloc)),
                     "median_known_loci": float(np.median(nknown)),
                     "median_novel_loci": float(np.median(nnovel)),
                     "loci_p05": float(np.percentile(nloc, 5)),
                     "loci_p95": float(np.percentile(nloc, 95)),
                     "novel_p05": float(np.percentile(nnovel, 5)),
                     "novel_p95": float(np.percentile(nnovel, 95))})
        print(f"  melanoma @ {name} power (Neff {target:,.0f}): median "
              f"{np.median(nloc):.0f} significant LOCI "
              f"[5-95%: {np.percentile(nloc,5):.0f}-{np.percentile(nloc,95):.0f}] "
              f"= {np.median(nknown):.0f} known / {np.median(nnovel):.0f} novel "
              f"(novel 5-95% [{np.percentile(nnovel,5):.0f}, {np.percentile(nnovel,95):.0f}])")
    pd.DataFrame(rows).to_csv(f"{MR}/85d_matched_power.tsv", sep="\t", index=False)

    # ------------------------------------------------- distances, for transparency
    print("\n=== distance from each significant HCC locus to its nearest known lead SNP ===")
    for tag, df in (("HCC_high", hi), ("HCC_low", lo)):
        for L, s in df[df.fdr < 0.05].groupby("locus"):
            ch = str(s.chr.iloc[0])
            arr = kn_by_chr.get(ch, np.array([]))
            best, bestrs = None, ""
            for pos in s.pos.unique():
                for _, kr in known[known["chr"] == ch].iterrows():
                    dd = abs(int(kr.pos) - int(pos))
                    if best is None or dd < best:
                        best, bestrs = dd, kr.rsid
            print(f"  {tag:<9} {L:<8} genes={','.join(sorted(set(s.gene_id)))} "
                  f"nearest known = {bestrs} at {best/1000:.0f} kb"
                  if best is not None else f"  {tag} {L}: no known locus on chromosome")

    # -------------------------------- what happens to each locus at the other power
    print("\n=== cross-power behaviour of every significant HCC locus ===")
    for src, dst, sname, dname in ((lo, hi, "HCC_low", "HCC_high"), (hi, lo, "HCC_high", "HCC_low")):
        for g in sorted(set(src[src.fdr < 0.05].gene_id)):
            o = dst[dst.gene_id == g]
            if len(o):
                print(f"  {g} significant in {sname}: in {dname} min P = {o.p_mr.min():.3g}, "
                      f"min FDR = {o.fdr.min():.3g} -> "
                      f"{'still FDR<0.05' if o.fdr.min() < 0.05 else 'DROPS OUT'}")
            else:
                print(f"  {g} significant in {sname}: absent from {dname} strict set")
    print("\nwrote 85a-85d")


if __name__ == "__main__":
    main()
