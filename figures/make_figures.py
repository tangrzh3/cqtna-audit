"""Main figures — methodological findings + candidate evidence"""
import csv, gzip, os, math, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pyarrow.parquet as pq

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "figures")
PARQ = r"D:/Downloads/CD4_eqtl_step1_clean"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": .8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150,
})
CAT = {
    "潜在新位点":      ("#3B7DD8", "Novel locus"),
    "痣数目通路":      ("#E8A33D", "Nevus-count locus"),
    "色素/发色通路":   ("#C4453C", "Pigmentation locus"),
    "黑色素瘤已知位点": ("#7A5AA8", "Known melanoma locus"),
}
CHRLEN = {str(i): l for i, l in enumerate(
    [248956422,242193529,198295559,190214555,181538259,170805979,159345973,
     145138636,138394717,133797422,135086622,133275309,114364328,107043718,
     101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468], 1)}

def rd(fp):
    with open(os.path.join(MR, fp), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def f(x):
    try: return float(x)
    except Exception: return None

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300)
    plt.close(fig); print(name, "ok", flush=True)

# ---------------------------------------------------------------- Fig A
def fig_A():
    rows = rd("13_meta_locus_annotation.tsv")
    off, cum = {}, 0
    for c in map(str, range(1, 23)):
        off[c] = cum; cum += CHRLEN[c]
    pts = collections.defaultdict(list)
    for r in rows:
        p = f(r["pval"])
        if p is None: continue
        ch, pos = r["SNP"].split(":")
        if ch not in off: continue
        pts[r["category"]].append((off[ch] + int(pos), -math.log10(max(p, 1e-320)),
                                   r["SYMBOL"], f(r["FDR"])))
    CAP = 30.0                      # 截断，避免 p 下溢点（-log10=320）压扁全图
    fig, ax = plt.subplots(figsize=(11, 4.4))
    for c in map(str, range(1, 23)):
        if int(c) % 2 == 0:
            ax.axvspan(off[c], off[c] + CHRLEN[c], color="#F3F3F3", zorder=0)
    for cat, (col, lab) in CAT.items():
        d = pts.get(cat, [])
        if not d: continue
        low = [(q[0], q[1]) for q in d if q[1] <= CAP]
        hig = [(q[0], CAP) for q in d if q[1] > CAP]
        ax.scatter([q[0] for q in low], [q[1] for q in low], s=10, c=col,
                   alpha=.75, edgecolors="none", label=f"{lab} (n={len(d)})", zorder=3)
        if hig:
            ax.scatter([q[0] for q in hig], [q[1] for q in hig], s=34, c=col,
                       marker="^", edgecolors="black", linewidths=.4, zorder=4)
    sig = [q[1] for cat in pts for q in pts[cat] if q[3] is not None and q[3] < .05]
    if sig:
        ax.axhline(min(sig), ls="--", lw=.8, c="#666", zorder=2)
        ax.text(cum * .997, min(sig) + .5, "FDR = 0.05", ha="right", fontsize=7.5, color="#555")
    ax.axhline(CAP, ls=":", lw=.7, c="#AAA", zorder=2)
    ax.text(cum * .003, CAP + .5, r"$\blacktriangle$  $P$ below float precision",
            fontsize=7, color="#555")

    # 标注：先按 1 Mb 将 FDR 显著基因聚成位点，同一位点合并标注
    best = {}
    for cat in pts:
        for x, y, s, fd in pts[cat]:
            if fd is None or fd >= .05 or not s: continue
            if s not in best or y > best[s][1]:
                best[s] = (x, min(y, CAP), cat)
    clusters = []
    for s, (x, y, cat) in sorted(best.items(), key=lambda kv: kv[1][0]):
        if clusters and abs(x - clusters[-1]["x"]) < 1_500_000:
            c = clusters[-1]
            c["genes"].append(s); c["x"] = max(c["x"], x)
            c["y"] = max(c["y"], y); c["cat"] = cat if y >= c["y"] else c["cat"]
        else:
            clusters.append({"x": x, "y": y, "genes": [s], "cat": cat})
    for k, c in enumerate(clusters):
        g = c["genes"]
        txt = g[0] if len(g) == 1 else f"{len(g)} genes\n({g[0]} …)"
        ax.annotate(txt, (c["x"], c["y"]), textcoords="offset points",
                    xytext=(0, 12 + 11 * (k % 2)), ha="center",
                    fontsize=6.8, style="italic", color=CAT[c["cat"]][0],
                    arrowprops=dict(arrowstyle="-", lw=.5, color="#AAA",
                                    shrinkA=0, shrinkB=1))
    ax.set_xticks([off[c] + CHRLEN[c] / 2 for c in map(str, range(1, 23))])
    ax.set_xticklabels(range(1, 23), fontsize=7)
    ax.set_xlim(0, cum); ax.set_ylim(-1, CAP + 12)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(P_{\rm MR})$")
    ax.set_title("Locus attribution of CD4$^+$ T cell eQTL–melanoma MR signals",
                 fontsize=10, pad=8)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", markerscale=1.5,
              ncol=2, columnspacing=1.0)
    fig.tight_layout(); save(fig, "FigA_locus_attribution")

# ---------------------------------------------------------------- Fig B
def fig_B():
    col = {r["exposure"]: r for r in rd("16_coloc_meta_results.tsv")}
    smr = {r["gene"] + "|" + r["profile"]: r for r in rd("15_SMR_meta_results.tsv")}
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    hi = 0
    dat = []
    for e, c in col.items():
        s = smr.get(e)
        if not s: continue
        h4, ph = f(c["PP.H4"]), f(s["p_HEIDI"])
        if h4 is None or ph is None or ph <= 0: continue
        y = -math.log10(ph); hi = max(hi, y)
        dat.append((h4, y, CAT.get(c["category"], ("#999", ""))[0], c["SYMBOL"], c["category"]))
    ax.axvspan(0, .5, color="#FCEEED", zorder=0)
    ax.axhline(-math.log10(.05), ls="--", lw=.8, c="#666", zorder=2)
    ax.axvline(.7, ls="--", lw=.8, c="#666", zorder=2)
    ax.scatter([d[0] for d in dat], [d[1] for d in dat], s=18,
               c=[d[2] for d in dat], alpha=.8, edgecolors="none", zorder=3)
    n_lo = [d for d in dat if d[0] < .2]
    n_pass = [d for d in n_lo if d[1] < -math.log10(.05)]
    ax.text(.03, hi * .97,
            "coloc: distinct causal variants (PP.H4 < 0.2)\n"
            "HEIDI: not rejected\n"
            f"{len(n_pass)}/{len(n_lo)} = {len(n_pass)/len(n_lo):.0%} pass HEIDI",
            fontsize=8.5, va="top", color="#8B2E28")
    tag = {"VPS9D1-AS1", "CDK10", "SPATA33", "CHMP1A", "CTU2"}
    done = set()
    for h4, y, c, s, cat in dat:
        if h4 < .15 and y < -math.log10(.05) and s in tag and s not in done:
            done.add(s)
            ax.annotate(s, (h4, y), textcoords="offset points", xytext=(5, 2),
                        fontsize=6.5, style="italic", color="#8B2E28")
    ax.set_xlabel("coloc PP.H4  (posterior of a shared causal variant)")
    ax.set_ylabel(r"$-\log_{10}(P_{\rm HEIDI})$   (lower = HEIDI passes)")
    ax.set_title("HEIDI fails to reject LD-confounded signals identified by coloc",
                 fontsize=10, pad=8)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=5, color=c, label=l)
                       for c, l in CAT.values()], frameon=False, fontsize=7.5, loc="upper right")
    fig.tight_layout(); save(fig, "FigB_coloc_vs_HEIDI")

# ---------------------------------------------------------------- Fig C
def region(path, ch, lo, hi, ipv="pval"):
    xs, ys = [], []
    with gzip.open(os.path.join(MR, path), "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ic, ip, ipvi = h.index("#chrom"), h.index("pos"), h.index(ipv)
        for line in fh:
            g = line.rstrip("\n").split("\t")
            if g[ic] != ch: continue
            p = int(g[ip])
            if not (lo <= p <= hi): continue
            v = f(g[ipvi])
            if v is None: continue
            xs.append(p); ys.append(-math.log10(max(v, 1e-320)))
    return np.array(xs), np.array(ys)

def eqtl_region(prof, gid):
    t = pq.read_table(os.path.join(PARQ, f"{prof}_step1_clean.parquet"),
                      columns=["gene_id", "pos", "pval"],
                      filters=[("gene_id", "==", gid)]).to_pydict()
    x = np.array(t["pos"]); y = -np.log10(np.maximum(np.array(t["pval"]), 1e-320))
    return x, y

def fig_C():
    panels = [
        ("PARP1",   "ENSG00000143799", "CD4_Naive_stim_40h", "1",  0.92, 0.05),
        ("ZFYVE19", "ENSG00000166140", "CD4_Naive_stim_5d",  "15", None, 0.99),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 5.6), sharex="row")
    for i, (sym, gid, prof, ch, h4_old, h4_new) in enumerate(panels):
        ex, ey = eqtl_region(prof, gid)
        lo, hi = int(ex.min()), int(ex.max())
        peak_e = ex[np.argmax(ey)]
        fx, fy = region("finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz", ch, lo, hi)
        mx, my = region("meta_melanoma_final.tsv.gz", ch, lo, hi)
        for j, (X, Y, ttl, col) in enumerate([
                (ex, ey, f"{sym} cis-eQTL  ({prof.replace('CD4_','').replace('_',' ')})", "#2E7D5B"),
                (fx, fy, "Melanoma GWAS — FinnGen (5,753 cases)", "#8C8C8C"),
                (mx, my, "Melanoma GWAS — meta (12,530 cases)", "#C4453C")]):
            ax = axes[i, j]
            ax.scatter(X / 1e6, Y, s=5, c=col, alpha=.55, edgecolors="none")
            ax.axvline(peak_e / 1e6, ls=":", lw=1.1, c="#2E7D5B", zorder=5)
            if len(Y):
                pk = X[np.argmax(Y)]
                ax.scatter([pk / 1e6], [Y.max()], s=42, facecolors="none",
                           edgecolors="black", lw=1.1, zorder=6)
                if j:
                    d = abs(pk - peak_e) / 1000
                    ax.annotate(f"{d:.0f} kb from eQTL peak", (pk / 1e6, Y.max()),
                                textcoords="offset points", xytext=(0, -14),
                                ha="center", fontsize=7)
            ax.set_title(ttl, fontsize=8.5)
            if j == 0: ax.set_ylabel(r"$-\log_{10}P$")
            ax.set_xlabel(f"Chromosome {ch} (Mb)", fontsize=8)
        note = (f"PP.H4  {h4_old:.2f} → {h4_new:.2f}" if h4_old is not None
                else f"PP.H4 → {h4_new:.2f}")
        axes[i, 2].text(.97, .95, note, transform=axes[i, 2].transAxes,
                        ha="right", va="top", fontsize=9, weight="bold",
                        color="#8B2E28" if h4_old else "#1B5E20")
    # wording fixed to match manuscript 2.4: the two rounds differ in resolution AND
    # cohort, and the multiple-signal analysis adjudicates only PARP1, so "power
    # determines resolution" (the earlier title) is not a claim this figure supports
    fig.suptitle("Posterior support moves with the outcome dataset — in both directions",
                 fontsize=10.5, y=.99)
    fig.tight_layout(rect=[0, 0, 1, .96]); save(fig, "FigC_power_dependence")

# ---------------------------------------------------------------- Fig D
def fig_D():
    genes = ["TPI1", "SMC2", "ZFYVE19", "SPSB2", "HLA-C", "KIAA0040"]
    layers = ["MR\n(4-fold)", "coloc\nPP.H4", "Not affecting\nnevi", "No\npleiotropy",
              "scRNA ICB\nresponse", "Bulk cohort\nreplication", "Druggability\n& safety"]
    M = {"TPI1":[1,.51,1,1,1,.8,.3], "SMC2":[1,.86,1,.6,0,.2,.1],
         "ZFYVE19":[1,.99,0,1,0,0,1], "SPSB2":[1,.53,1,1,0,0,.6],
         "HLA-C":[1,.54,1,0,0,0,0], "KIAA0040":[.5,.85,0,0,.3,0,.6]}
    A = np.array([M[g] for g in genes])
    fig, ax = plt.subplots(figsize=(7.6, 3.5))
    im = ax.imshow(A, cmap="RdYlBu_r", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(layers))); ax.set_xticklabels(layers, fontsize=7.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=9.5, style="italic")
    for i in range(len(genes)):
        for j in range(len(layers)):
            v = A[i, j]
            ax.text(j, i, "✓" if v >= .8 else ("~" if v >= .4 else "✗"),
                    ha="center", va="center", fontsize=11,
                    color="white" if (v > .68 or v < .18) else "black")
    ax.set_xticks(np.arange(-.5, len(layers), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.5); ax.tick_params(which="minor", size=0)
    ax.set_title("Multi-layer evidence matrix for candidate genes", fontsize=10, pad=8)
    cb = fig.colorbar(im, ax=ax, fraction=.02, pad=.02)
    cb.set_label("Evidence strength", fontsize=8); cb.ax.tick_params(labelsize=7)
    fig.tight_layout(); save(fig, "FigD_evidence_matrix")

if __name__ == "__main__":
    fig_A(); fig_B(); fig_C(); fig_D()
