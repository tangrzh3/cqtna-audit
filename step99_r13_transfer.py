#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 99 -- the FinnGen R13 transfer test.

Design, predictions and interpretation table: manuscript/PREREG_r13_transfer.md (S25)
Registered BEFORE any R13 association statistic was read.

This is NOT a sixth point on the R8->R12 power trajectory. R13 fails S9's step-0
gate: the endpoint changed from C3_MELANOMA_SKIN_EXALLC to C3_MELANOMA_SKIN_WIDE
and the control-exclusion rule with it. What is asked instead is the question a
reader would ask: re-run this audit on the release FinnGen ships today and does
the answer change?

Four cells (prereg §4), each paired with its R12 comparator computed by the SAME
code path so the comparison is not confounded by the machinery:

  C1  Soskic CD4    x  R13 melanoma      <- main cell
  C2  eQTLGen blood x  R13 melanoma
  C3  Soskic CD4    x  R13 HCC
  C4  eQTLGen blood x  R13 HCC

Everything else is held fixed: instruments, allele matching, Wald ratio, the
BH-FDR testing family, the known-locus reference lists, the 1 Mb single-linkage
locus rule. Only the outcome file moves.

Outputs: 99a_r13_trajectory.tsv, 99b_r13_lists.tsv, 99c_r13_attribution.tsv,
         99d_r13_records.tsv.gz -- per-record, so a later partition
         change can be recomputed without rescanning the sumstats
"""
import csv
import gzip
import os
import sys
import time
from math import lgamma, exp
from statistics import NormalDist

import numpy as np
import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
ND = NormalDist()
FDR_MAIN = 0.05
NOVEL = "潜在新位点"
LOCUS_KB = KNOWN_KB = 1000

R13_MEL = f"{MR}/r13/finngen_R13_C3_MELANOMA_SKIN_WIDE.gz"
R13_HCC = f"{MR}/r13/finngen_R13_C3_HEPATOCELLU_CARC_WIDE.gz"
R12_MEL = f"{MR}/finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz"
R12_HCC = f"{MR}/hcc/finngen_R12_C3_HEPATOCELLU_CARC_EXALLC.gz"

# ⚠ 92c_locus_annotated.tsv was built against meta_melanoma_final.tsv.gz (the
# FinnGen+Rashkin meta, N_eff 49,337), NOT against FinnGen R12 (N_eff 22,668).
# Its beta_out/se_out columns are therefore NOT the R12 comparator for the
# eQTLGen cells, and using them as one compares two releases AND two datasets at
# the same time -- exactly the confound this whole paper is about. The R12
# comparator for C2 is computed here from the R12 file; the meta cell is carried
# alongside as a separately labelled reference row (this is S22's 4.44x).

# metadata only (manifest), registered in prereg §1
META = {"R12_mel": (5753, 378749), "R13_mel": (6226, 372159),
        "R12_hcc": (947, 378749), "R13_hcc": (1070, 372159)}


# ----------------------------------------------------------------- statistics
def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, min(r1, c1) + 1)), 1.0)


def two_sided(z):
    return np.array([2 * ND.cdf(-abs(v)) if abs(v) < 37 else 0.0 for v in z])


def n_eff(nc, nk):
    return 4.0 / (1.0 / nc + 1.0 / nk)


def assign_loci(chrs, poss):
    """1 Mb single-linkage, identical rule to step59b / step94c."""
    order = sorted(range(len(chrs)), key=lambda i: (str(chrs[i]), int(poss[i])))
    out = [None] * len(chrs)
    lc, lp, lid = None, None, 0
    for i in order:
        c, p = str(chrs[i]), int(poss[i])
        if c != lc or p - lp > LOCUS_KB * 1000:
            lid += 1
        out[i] = lid
        lc, lp = c, p
    return out


def known_flagger(csv_path):
    """position-based known-locus membership, 1 Mb window, GRCh38 throughout"""
    k = pd.read_csv(csv_path)
    k["chr"] = k["chr"].astype(str)
    KN = {c: np.sort(s.pos.values) for c, s in k.groupby("chr")}

    def is_known(ch, pos):
        arr = KN.get(str(ch))
        if arr is None or not len(arr):
            return False
        i = np.searchsorted(arr, pos)
        return any(0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= KNOWN_KB * 1000
                   for j in (i - 1, i))
    return is_known, len(k)


def attribution(chrs, poss, fdr, known_flags, label):
    """share of FDR<0.05 independent loci carrying a known lead SNP for the disease"""
    loci = assign_loci(chrs, poss)
    d = pd.DataFrame(dict(locus=loci, known=known_flags, fdr=fdr))
    bg = d.groupby("locus").agg(known=("known", "any"))
    BT, BK = len(bg), int(bg.known.sum())
    sg = d[d.fdr < FDR_MAIN].groupby("locus").agg(known=("known", "any"))
    ST, SK = len(sg), int(sg.known.sum())
    fold = ((SK / ST) / (BK / BT)) if ST and BK else np.nan
    p = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)) if ST else np.nan
    print(f"    [{label}] background known {BK}/{BT} = {100*BK/BT:.1f}%  |  "
          f"FDR<0.05 {SK}/{ST} known"
          f"{f' = {100*SK/ST:.1f}%' if ST else ''}  ->  "
          f"{fold if not np.isnan(fold) else float('nan'):.2f}x, P = {p:.3g}")
    return dict(cell=label, bg_loci=BT, bg_known=BK, sig_loci=ST, sig_known=SK,
                pct_known=round(100 * SK / ST, 1) if ST else np.nan,
                fold=round(fold, 3) if ST and BK else np.nan, fisher_p=p)


# ----------------------------------------------------------------- outcome scan
def scan(path, want, label):
    """One pass over a FinnGen sumstats file, keeping every row at a wanted
    (chr,pos). Multi-allelic positions carry several rows and the allele choice
    must be made later against the instrument, not here -- this project already
    has one indel-representation ambiguity on record (step59b).

    Extracts are cached: each of these files is ~0.8 GB and a re-run that only
    changes the analysis must not have to pay for the scan again."""
    os.makedirs(f"{MR}/r13_extracts", exist_ok=True)
    cache = f"{MR}/r13_extracts/{label}.tsv"
    if os.path.exists(cache):
        got = {}
        with open(cache, encoding="utf-8") as fh:
            rd = csv.reader(fh, delimiter="\t")
            next(rd)
            for c, p, ref, alt, b, s in rd:
                got.setdefault((c, int(p)), []).append((ref, alt, float(b), float(s)))
        if set(got) >= want:
            print(f"  {label}: {len(got):,} positions from cache")
            return got
        print(f"  {label}: cache covers {len(got):,} positions but "
              f"{len(want - set(got)):,} are missing; rescanning")

    got, n, t0 = {}, 0, time.time()
    with gzip.open(path, "rt") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        ix = {c: i for i, c in enumerate(hdr)}
        ic, ip = ix["#chrom"], ix["pos"]
        ir, ia = ix["ref"], ix["alt"]
        ib, ise = ix["beta"], ix["sebeta"]
        for row in rd:
            n += 1
            try:
                key = (row[ic], int(row[ip]))
            except (ValueError, IndexError):
                continue
            if key in want:
                try:
                    got.setdefault(key, []).append(
                        (row[ir].upper(), row[ia].upper(),
                         float(row[ib]), float(row[ise])))
                except ValueError:
                    pass
            if n % 5_000_000 == 0:
                print(f"    {label}: {n/1e6:.0f}M rows, {len(got):,} positions, "
                      f"{time.time()-t0:.0f}s", flush=True)
    multi = sum(1 for v in got.values() if len(v) > 1)
    print(f"  {label}: {n:,} rows scanned, {len(got):,} positions matched "
          f"({multi} multi-allelic), {time.time()-t0:.0f}s", flush=True)
    with open(cache, "w", newline="", encoding="utf-8") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["chr", "pos", "ref", "alt", "beta", "sebeta"])
        for (c, p), v in got.items():
            for ref, alt, b, s in v:
                w.writerow([c, p, ref, alt, b, s])
    return got


def pick(cands, ea, oa):
    """allele-aware selection; returns (beta, se) oriented to the exposure
    effect allele, or None when no row at this position carries the pair"""
    for ref, alt, b, s in cands:
        if ea == alt and oa == ref:
            return b, s
        if ea == ref and oa == alt:
            return -b, s
    return None


# ----------------------------------------------------------------- exposures
def load_soskic_melanoma():
    """the FinnGen R12 melanoma round, exactly the table step59b reads"""
    d = pd.read_csv(f"{MR}/06_locus_annotation.tsv", sep="\t")
    d = d[["exposure", "SYMBOL", "category", "SNP", "beta.exposure",
           "beta.outcome", "se.outcome"]].dropna()
    d = d[(d["beta.exposure"] != 0) & (d["se.outcome"] > 0)].reset_index(drop=True)
    d["novel"] = d["category"].eq(NOVEL)
    d[["chr", "pos"]] = d.SNP.str.split(":", expand=True).iloc[:, :2]
    d["pos"] = d["pos"].astype(int)

    alle = {}
    with open(f"{MR}/01_harmonised_all.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            alle[(r["exposure"], r["SNP"])] = (r["effect_allele.outcome"].upper(),
                                               r["other_allele.outcome"].upper())
    d["ea"] = [alle.get((e, s), ("", ""))[0] for e, s in zip(d.exposure, d.SNP)]
    d["oa"] = [alle.get((e, s), ("", ""))[1] for e, s in zip(d.exposure, d.SNP)]
    d = d[(d.ea != "") & (d.oa != "")].reset_index(drop=True)
    print(f"Soskic x melanoma instruments: {len(d):,} records / "
          f"{d.SYMBOL.nunique():,} genes")
    return d


def load_soskic_hcc():
    d = pd.read_csv(f"{MR}/85a_HCC_low_annotated.tsv", sep="\t")
    d = d[["exposure", "gene_id", "chr", "pos", "ea", "oa", "beta",
           "beta_out", "se_out", "known"]].dropna()
    d = d[(d["beta"] != 0) & (d["se_out"] > 0)].reset_index(drop=True)
    d["chr"] = d["chr"].astype(str)
    d["pos"] = d["pos"].astype(int)
    d["ea"] = d["ea"].str.upper()
    d["oa"] = d["oa"].str.upper()
    print(f"Soskic x HCC instruments: {len(d):,} records")
    return d


def load_eqtlgen():
    d = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
    d = d[["gene_id", "symbol", "chr", "pos", "beta_exp", "beta_out",
           "se_out", "known"]].dropna()
    d["chr"] = d["chr"].astype(str)
    d["pos"] = d["pos"].astype(int)
    print(f"eQTLGen instruments: {len(d):,} records / {d.gene_id.nunique():,} genes")
    return d


# ------------------------------------------------------- per-record persistence
# step99 used to keep only the summary row for each cell. When the locus
# partition changed (PREREG_locus_partition.md) there was then no way to
# recompute these cells short of rescanning 1.5 GB of sumstats, and no way for a
# third party to check them at all. The standardised per-record table now goes to
# 99d_r13_records.tsv.gz with the columns that record fixes in its section 5:
#   cell, record_id, gene, chr, pos, p
RECORDS = []


def genes_of(sub):
    for col in ("SYMBOL", "symbol", "gene_id", "gene"):
        if col in sub.columns:
            return [str(x) for x in sub[col]]
    return ["unnamed_%d" % (i + 1) for i in range(len(sub))]


def stash(label, sub, pvals):
    for i, (g, c, po, pv) in enumerate(
            zip(genes_of(sub), sub.chr, sub.pos, pvals), 1):
        RECORDS.append((label, i, g, str(c), int(po), float(pv)))


# ----------------------------------------------------------------- cells
def soskic_melanoma_cells(d, got13, hcc_known_fn):
    """C1 and its R12 comparator: the trajectory-style list analysis + attribution.

    Main analysis is restricted to instruments present in BOTH releases (prereg
    §3); the secondary analysis uses each release's own coverage.
    """
    keep13, b13, s13 = [], [], []
    for i, r in enumerate(d.itertuples()):
        hit = pick(got13.get((r.chr, r.pos), []), r.ea, r.oa)
        if hit is None or hit[1] <= 0:
            continue
        keep13.append(i)
        b13.append(hit[0])
        s13.append(hit[1])
    matched13 = len(keep13)
    print(f"  C1 coverage: R12 {len(d):,} -> R13 {matched13:,} "
          f"({100*matched13/len(d):.1f}% of R12)")

    rows, lists, attrib = [], [], []
    for mode in ("common", "own"):
        idx13 = keep13
        idx12 = keep13 if mode == "common" else list(range(len(d)))
        for rel, idx, bo, so in (
                ("R12", idx12, None, None),
                ("R13", idx13, b13, s13)):
            sub = d.loc[idx].reset_index(drop=True)
            be = sub["beta.exposure"].values
            if rel == "R12":
                bout = sub["beta.outcome"].values
                sout = sub["se.outcome"].values
            else:
                bout, sout = np.array(bo), np.array(so)
            z = (bout / be) / (sout / np.abs(be))
            p = two_sided(z)
            f = bh(p)
            hit = f < FDR_MAIN
            loci = np.array(assign_loci(sub.chr.tolist(), sub.pos.tolist()))
            nc, nk = META[f"{rel}_mel"]
            rows.append(dict(mode=mode, cell="C1_Soskic_melanoma", release=rel,
                             cases=nc, controls=nk, n_eff=round(n_eff(nc, nk)),
                             n_tested=len(sub), n_nominal=int((p < .05).sum()),
                             n_hits=int(hit.sum()),
                             n_genes=sub.SYMBOL[hit].nunique(),
                             n_loci=len(set(loci[hit])),
                             n_novel_genes=sub.SYMBOL[hit & sub.novel.values].nunique()))
            if mode == "common":
                lists.append(dict(cell="C1_Soskic_melanoma", release=rel,
                                  genes_FDR05=",".join(sorted(set(sub.SYMBOL[hit]))),
                                  novel_genes_FDR05=",".join(
                                      sorted(set(sub.SYMBOL[hit & sub.novel.values]))) or "-"))
                stash(f"C1 Soskic x melanoma {rel}", sub, p)
                attrib.append(attribution(sub.chr.tolist(), sub.pos.tolist(), f,
                                          (~sub.novel).tolist(),
                                          f"C1 Soskic x melanoma {rel}"))
                if rel == "R13":
                    # NC: score the same cell against the MISMATCHED (HCC) list
                    nc_flags = [hcc_known_fn(c, p_) for c, p_ in zip(sub.chr, sub.pos)]
                    attrib.append(attribution(sub.chr.tolist(), sub.pos.tolist(), f,
                                              nc_flags,
                                              "NC C1 R13 scored with HCC list"))
    return rows, lists, attrib, matched13


def scanned_cell(inst, got, known_flags, label, allele_aware):
    """C2/C3/C4 against a scanned outcome file. Attribution is orientation-free,
    so only |z| matters and the allele pass is a filter, not a sign fix."""
    keep, z = [], []
    for i, r in enumerate(inst.itertuples()):
        cands = got.get((r.chr, r.pos), [])
        if not cands:
            continue
        if allele_aware:
            hit = pick(cands, r.ea, r.oa)
        else:
            # 92c carries no alleles; step94d's convention (first row) is kept
            hit = (cands[0][2], cands[0][3])
        if hit is None or hit[1] <= 0:
            continue
        keep.append(i)
        z.append(hit[0] / hit[1])
    sub = inst.iloc[keep].reset_index(drop=True)
    pv = two_sided(np.array(z))
    f = bh(pv)
    stash(label, sub, pv)
    res = attribution(sub.chr.tolist(), sub.pos.tolist(), f,
                      [known_flags[i] for i in keep], label)
    res["n_records"] = len(sub)
    res["coverage_pct"] = round(100 * len(sub) / len(inst), 1)
    return res


def precomputed_cell(inst, bcol, scol, known_flags, label):
    """a cell whose outcome columns are already in the table (the S22 meta row)"""
    pv = two_sided(inst[bcol].values / inst[scol].values)
    f = bh(pv)
    stash(label, inst, pv)
    res = attribution(inst.chr.tolist(), inst.pos.tolist(), f, known_flags, label)
    res["n_records"] = len(inst)
    res["coverage_pct"] = 100.0
    return res


def main():
    mel_known_fn, n_mel_known = known_flagger(f"{MR}/landi2020_known_loci_grch38.csv")
    hcc_known_fn, n_hcc_known = known_flagger(f"{MR}/84a_hcc_known_loci_grch38.csv")
    print(f"known-locus reference lists: melanoma {n_mel_known}, HCC {n_hcc_known}\n")

    d_mel = load_soskic_melanoma()
    d_hcc = load_soskic_hcc()
    d_eq = load_eqtlgen()

    want_mel = set(zip(d_mel.chr, d_mel.pos)) | set(zip(d_eq.chr, d_eq.pos))
    want_hcc = set(zip(d_hcc.chr, d_hcc.pos)) | set(zip(d_eq.chr, d_eq.pos))
    print(f"\npositions to look up: melanoma {len(want_mel):,}, HCC {len(want_hcc):,}\n")

    for p in (R13_MEL, R13_HCC, R12_MEL, R12_HCC):
        if not os.path.exists(p):
            sys.exit(f"missing outcome file: {p}")

    print("=== scanning outcomes ===", flush=True)
    got_mel13 = scan(R13_MEL, want_mel, "R13_melanoma")
    got_hcc13 = scan(R13_HCC, want_hcc, "R13_HCC")
    got_mel12 = scan(R12_MEL, set(zip(d_eq.chr, d_eq.pos)), "R12_melanoma")
    got_hcc12 = scan(R12_HCC, set(zip(d_eq.chr, d_eq.pos)), "R12_HCC")

    print("\n=== C1  Soskic CD4 x melanoma ===")
    rows, lists, attrib, matched13 = soskic_melanoma_cells(d_mel, got_mel13, hcc_known_fn)

    mel_flags = d_eq.known.astype(bool).tolist()
    hcc_flags = [hcc_known_fn(c, p) for c, p in zip(d_eq.chr, d_eq.pos)]

    print("\n=== C2  eQTLGen x melanoma ===")
    # the S22 reference row: same instruments, but the FinnGen+Rashkin META
    # outcome. Kept because it is the published 4.44x, and labelled so that it is
    # never read as R12.
    attrib.append(precomputed_cell(d_eq, "beta_out", "se_out", mel_flags,
                                   "REF eQTLGen x melanoma META (S22)"))
    attrib.append(scanned_cell(d_eq, got_mel12, mel_flags,
                               "C2 eQTLGen x melanoma R12", allele_aware=False))
    attrib.append(scanned_cell(d_eq, got_mel13, mel_flags,
                               "C2 eQTLGen x melanoma R13", allele_aware=False))

    print("\n=== C3  Soskic CD4 x HCC ===")
    hcc_inst_flags = d_hcc.known.astype(bool).tolist()
    attrib.append(precomputed_cell(d_hcc, "beta_out", "se_out", hcc_inst_flags,
                                   "C3 Soskic x HCC R12"))
    attrib.append(scanned_cell(d_hcc, got_hcc13, hcc_inst_flags,
                               "C3 Soskic x HCC R13", allele_aware=True))

    print("\n=== C4  eQTLGen x HCC ===")
    attrib.append(scanned_cell(d_eq, got_hcc12, hcc_flags,
                               "C4 eQTLGen x HCC R12", allele_aware=False))
    attrib.append(scanned_cell(d_eq, got_hcc13, hcc_flags,
                               "C4 eQTLGen x HCC R13", allele_aware=False))

    # ------------------------------------------------------------------ output
    traj = pd.DataFrame(rows)
    traj.to_csv(f"{MR}/99a_r13_trajectory.tsv", sep="\t", index=False)
    pd.DataFrame(lists).to_csv(f"{MR}/99b_r13_lists.tsv", sep="\t", index=False)
    pd.DataFrame(attrib).to_csv(f"{MR}/99c_r13_attribution.tsv", sep="\t", index=False)
    pd.DataFrame(RECORDS, columns=["cell", "record_id", "gene", "chr", "pos", "p"]) \
      .to_csv(f"{MR}/99d_r13_records.tsv.gz", sep="\t", index=False,
              compression="gzip")

    pd.set_option("display.width", 220)
    print("\n" + "=" * 92)
    print("C1 trajectory (main analysis = instruments shared by R12 and R13)")
    print("=" * 92)
    print(traj[traj["mode"] == "common"].to_string(index=False))
    print("\nsecondary (each release's own coverage)")
    print(traj[traj["mode"] == "own"].to_string(index=False))

    print("\n" + "=" * 92)
    print("FDR<0.05 gene lists")
    print("=" * 92)
    for r in lists:
        print(f"  {r['release']}: {r['genes_FDR05']}")
        print(f"       novel-locus genes: {r['novel_genes_FDR05']}")

    # ------------------------------------------------------------------ prereg
    c = traj[(traj["mode"] == "common")].set_index("release")
    r12, r13 = c.loc["R12"], c.loc["R13"]
    g12 = set(lists[0]["genes_FDR05"].split(",")) if lists[0]["genes_FDR05"] else set()
    g13 = set(lists[1]["genes_FDR05"].split(",")) if lists[1]["genes_FDR05"] else set()
    REF6 = {"CDK10", "CHMP1A", "CTU2", "PARP1", "SPATA33", "VPS9D1-AS1"}
    MC1R = {"VPS9D1-AS1", "CDK10", "SPATA33", "CHMP1A"}
    a = {x["cell"]: x for x in attrib}
    fold_c1 = a["C1 Soskic x melanoma R13"]["fold"]

    print("\n" + "=" * 92)
    print("PRE-REGISTERED CHECKS  (S25 §5-§7)")
    print("=" * 92)
    pc1 = (int(r12.n_hits) == 10 and int(r12.n_genes) == 6 and int(r12.n_loci) == 2
           and g12 == REF6)
    print(f"  PC1 R12 reproduction  : {int(r12.n_hits)} records / {int(r12.n_genes)} genes"
          f" / {int(r12.n_loci)} loci, list=={sorted(REF6)==sorted(g12)}"
          f"   -> {'PASS' if pc1 else 'FAIL'}")
    pc2 = MC1R <= g13
    print(f"  PC2 MC1R cluster in R13: {sorted(MC1R & g13)}   -> {'PASS' if pc2 else 'FAIL'}")
    cov = 100 * matched13 / len(d_mel)
    print(f"  PC3 coverage           : {cov:.1f}% of R12 instruments"
          f"   -> {'PASS' if cov >= 95 else 'FAIL'}")
    nc = a["NC C1 R13 scored with HCC list"]
    ncp = nc["fisher_p"]
    print(f"  NC  mismatched list    : fold {nc['fold']}, P = {ncp:.3g}"
          f"   -> {'PASS (null)' if not (ncp < .05 and (nc['fold'] or 0) > 1) else 'FAIL'}")

    checks = [
        ("P1 n records", int(r13.n_hits), 11, 9, 14),
        ("P2 n genes", int(r13.n_genes), 6, 5, 8),
        ("P3 recovery of R12's 6", len(REF6 & g13), 6, 5, 6),
        ("P4 novel-locus genes", int(r13.n_novel_genes), 0, 0, 1),
        ("P6 n independent loci", int(r13.n_loci), 2, 2, 4),
    ]
    print()
    for name, obs, pt, lo, hi in checks:
        ok = lo <= obs <= hi
        print(f"  {name:<24} observed {obs:>3}   registered {pt} [{lo},{hi}]"
              f"   -> {'MET' if ok else 'DEVIATION'}")
    ok5 = (fold_c1 is not None) and (not np.isnan(fold_c1)) and fold_c1 > 1
    print(f"  {'P5 attribution fold':<24} observed {fold_c1}   registered >1"
          f"   -> {'MET' if ok5 else 'DEVIATION'}")

    print("\n" + "=" * 92)
    print("attribution across all cells (R12 comparator vs R13)")
    print("=" * 92)
    print(pd.DataFrame(attrib)[["cell", "bg_known", "bg_loci", "sig_known",
                                "sig_loci", "pct_known", "fold", "fisher_p",
                                "coverage_pct"]].to_string(index=False))
    print("\nwrote 99a / 99b / 99c / 99d")


if __name__ == "__main__":
    main()
