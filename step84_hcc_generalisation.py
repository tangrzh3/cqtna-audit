"""Step 84 -- Part I generalisation to a second tumour type (hepatocellular carcinoma).

Pre-registered in manuscript/PREREG_hcc_generalisation.md (written before any HCC
outcome row was read). This script performs stage 1: build the GRCh38 known-locus
reference for HCC, extract the unchanged CD4 dynamic-eQTL instruments from both HCC
outcomes, harmonise, and run Wald-ratio MR.

Exposure side is byte-identical to the melanoma analysis: CD4_dynamic_top_eqtl_
instruments_all.csv, 8,064 instruments, GRCh38, one instrument per exposure.

Outcomes:
  HCC-high  GCST90809296 (harmonised, GRCh38)  3,748 EUR cases / 1,861,536 controls
  HCC-low   FinnGen R12 C3_HEPATOCELLU_CARC_EXALLC  947 cases / 378,749 controls

Outputs 84a-84d.
"""
import csv, gzip, os, sys, math

MR   = os.path.dirname(os.path.abspath(__file__))
HCC  = os.path.join(MR, "hcc")
EXP  = os.path.join(MR, "CD4_dynamic_top_eqtl_instruments_all.csv")
HIGH = os.path.join(HCC, "GCST90809296.h.tsv.gz")
LOW  = os.path.join(HCC, "finngen_R12_C3_HEPATOCELLU_CARC_EXALLC.gz")
KNOWN_RS = os.path.join(HCC, "hcc_known_rsids.csv")

# European lead variants from the outcome GWAS's own paper (Table 1, JHEP Rep 2025,
# PMID 40823170).  Added explicitly so the reference set cannot miss the loci that
# the primary outcome itself is powered to find.
PAPER_EUR_LEADS = {
    "rs2642442": "MTARC1", "rs7628416": "KLF15",  "rs4089": "HSD17B13",
    "rs10069690": "TERT",  "rs144861591": "HFE",  "rs58542926": "TM6SF2",
    "rs429358": "APOE",    "rs738408": "PNPLA3",
    # East Asian arm, kept for completeness of the reference set
    "rs9277534": "HLA-DPA1", "rs12971396": "IFNL4",
}


def log(*a):
    print(*a, flush=True)


def load_exposure():
    """One row per exposure (gene x cell x timepoint); GRCh38 chr:pos."""
    rows = []
    with open(EXP, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                b = float(r["beta.exposure"]); se = float(r["se.exposure"])
                p = float(r["pval.exposure"]); F = float(r["F_stat"])
                pos = int(r["pos"])
            except (ValueError, KeyError, TypeError):
                continue
            if b == 0:
                continue
            rows.append({
                "exposure": r["exposure"], "gene_id": r["gene_id"],
                "cell_type": r["cell_type"], "timepoint": r["timepoint"],
                "chr": r["chr"].replace("chr", ""), "pos": pos,
                "ea": r["effect_allele.exposure"].upper(),
                "oa": r["other_allele.exposure"].upper(),
                "beta": b, "se": se, "pval": p, "F": F,
            })
    log(f"exposure instruments: {len(rows)}  unique SNPs: {len({(x['chr'],x['pos']) for x in rows})}")
    return rows


def wanted_keys(exp):
    return {(x["chr"], x["pos"]) for x in exp}


def scan_finngen(keys, known_rs):
    """One pass: pull instrument rows by chr:pos AND map known rsIDs to GRCh38."""
    got, rsmap, n = {}, {}, 0
    with gzip.open(LOW, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        for line in f:
            n += 1
            p = line.rstrip("\n").split("\t")
            c, pos = p[ix["#chrom"]], p[ix["pos"]]
            k = (c, int(pos))
            if k in keys:
                try:
                    got[k] = {
                        "ea": p[ix["alt"]].upper(), "oa": p[ix["ref"]].upper(),
                        "beta": float(p[ix["beta"]]), "se": float(p[ix["sebeta"]]),
                        "pval": float(p[ix["pval"]]), "eaf": float(p[ix["af_alt"]]),
                        "rsid": p[ix["rsids"]],
                    }
                except ValueError:
                    pass
            rs = p[ix["rsids"]]
            if rs in known_rs and rs not in rsmap:
                rsmap[rs] = (c, int(pos))
    log(f"FinnGen HCC: {n} rows scanned, {len(got)} instruments matched, "
        f"{len(rsmap)}/{len(known_rs)} known rsIDs mapped")
    return got, rsmap


def scan_gcst(keys, known_rs, rsmap):
    got, n = {}, 0
    with gzip.open(HIGH, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        for line in f:
            n += 1
            p = line.rstrip("\n").split("\t")
            try:
                c, pos = p[ix["chromosome"]], int(p[ix["base_pair_location"]])
            except (ValueError, IndexError):
                continue
            k = (c, pos)
            if k in keys and k not in got:
                try:
                    eaf = p[ix["effect_allele_frequency"]]
                    got[k] = {
                        "ea": p[ix["effect_allele"]].upper(),
                        "oa": p[ix["other_allele"]].upper(),
                        "beta": float(p[ix["beta"]]), "se": float(p[ix["standard_error"]]),
                        "pval": float(p[ix["p_value"]]),
                        "eaf": float(eaf) if eaf not in ("", "NA") else float("nan"),
                        "rsid": p[ix["rsid"]],
                    }
                except (ValueError, IndexError):
                    pass
            rs = p[ix["rsid"]] if ix.get("rsid") is not None else ""
            if rs in known_rs and rs not in rsmap:
                rsmap[rs] = (c, pos)
    log(f"GCST90809296: {n} rows scanned, {len(got)} instruments matched, "
        f"known rsIDs mapped now {len(rsmap)}/{len(known_rs)}")
    return got


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def harmonise(e, o):
    """Align outcome to exposure effect allele. Returns (beta_out, se_out, action) or None."""
    ee, eo = e["ea"], e["oa"]
    oe, oo = o["ea"], o["oa"]
    if {ee, eo} == {oe, oo}:
        return (o["beta"] if oe == ee else -o["beta"]), o["se"], "direct"
    # try strand flip
    try:
        foe, foo = COMP[oe], COMP[oo]
    except KeyError:
        return None
    if {ee, eo} == {foe, foo}:
        return (o["beta"] if foe == ee else -o["beta"]), o["se"], "flipped"
    return None


def is_palindromic(a, b):
    return {a, b} in ({"A", "T"}, {"C", "G"})


def bh_fdr(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    out = [1.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        q = min(prev, pvals[i] * n / k)
        out[i] = q
        prev = q
    return out


def run_mr(exp, out, label):
    res, dropped = [], {"nomatch": 0, "harm": 0, "palin": 0}
    for e in exp:
        o = out.get((e["chr"], e["pos"]))
        if o is None:
            dropped["nomatch"] += 1
            continue
        h = harmonise(e, o)
        if h is None:
            dropped["harm"] += 1
            continue
        bo, so, action = h
        if action == "flipped" and is_palindromic(e["ea"], e["oa"]):
            dropped["palin"] += 1
            continue
        if so <= 0:
            continue
        b = bo / e["beta"]
        se = so / abs(e["beta"])          # first-order delta, as in the melanoma analysis
        z = bo / so                       # p depends only on the outcome, by construction
        p = math.erfc(abs(z) / math.sqrt(2))
        res.append({**e, "beta_out": bo, "se_out": so, "eaf_out": o.get("eaf"),
                    "rsid": o.get("rsid", ""), "b_mr": b, "se_mr": se,
                    "z": z, "p_mr": p, "action": action})
    log(f"[{label}] MR records: {len(res)}  dropped: {dropped}")
    return res


def main():
    exp = load_exposure()
    keys = wanted_keys(exp)

    known = {}
    with open(KNOWN_RS, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            known[r["rsid"]] = r.get("genes", "")
    for rs, g in PAPER_EUR_LEADS.items():
        known.setdefault(rs, g)
    log(f"known HCC lead rsIDs (GWAS Catalog p<5e-8 + outcome paper Table 1): {len(known)}")

    low_out, rsmap = scan_finngen(keys, set(known))
    high_out = scan_gcst(keys, set(known), rsmap)

    with open(os.path.join(MR, "84a_hcc_known_loci_grch38.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["rsid", "chr", "pos", "genes", "source"])
        for rs in sorted(known):
            if rs in rsmap:
                c, p = rsmap[rs]
                src = "paper_table1" if rs in PAPER_EUR_LEADS else "gwas_catalog"
                w.writerow([rs, c, p, known[rs], src])
    log(f"84a written: {len(rsmap)} known loci with GRCh38 coordinates")

    for tag, out in (("HCC_high_GCST90809296", high_out), ("HCC_low_FinnGenR12", low_out)):
        res = run_mr(exp, out, tag)
        strict = [r for r in res if r["pval"] < 5e-8 and r["F"] > 10]
        log(f"[{tag}] strict set: {len(strict)} records, "
            f"{len({(r['chr'], r['pos']) for r in strict})} unique SNPs")
        if strict:
            qs = bh_fdr([r["p_mr"] for r in strict])
            for r, q in zip(strict, qs):
                r["fdr"] = q
            log(f"[{tag}] FDR<0.05: {sum(1 for r in strict if r['fdr'] < 0.05)}  "
                f"nominal p<0.05: {sum(1 for r in strict if r['p_mr'] < 0.05)}")
        fn = os.path.join(MR, f"84b_mr_{tag}.tsv")
        cols = ["exposure", "gene_id", "cell_type", "timepoint", "chr", "pos", "rsid",
                "ea", "oa", "beta", "se", "pval", "F", "beta_out", "se_out", "eaf_out",
                "b_mr", "se_mr", "z", "p_mr", "fdr", "action"]
        with open(fn, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in sorted(strict, key=lambda x: x["p_mr"]):
                w.writerow(r)
        log(f"  wrote {fn}")


if __name__ == "__main__":
    main()
