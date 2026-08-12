"""Fig S1 — 分析流程示意（含每一步的实际数量与判废条件）

不是通用流程图：每个框里写真实数字，每个判废条件写在它拦掉东西的那一步旁边，
右侧一列标出该步骤对应的方法学发现。这样流程图本身就是论证的索引。
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = r"D:/R_ex/MR/figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150})

BLUE, RED, GREEN, GREY, DARK = "#3B7DD8", "#C4453C", "#2E7D5B", "#8B9096", "#333333"

fig, ax = plt.subplots(figsize=(11.5, 9.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")


def box(x, y, w, h, title, body, ec=DARK, fc="white", tcol=DARK, fs=8.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6",
                                fc=fc, ec=ec, lw=1.2, zorder=3))
    ax.text(x + w / 2, y + h - 1.6, title, ha="center", va="top",
            fontsize=fs + .7, weight="bold", color=tcol, zorder=4)
    ax.text(x + w / 2, y + h - 5.0, body, ha="center", va="top",
            fontsize=fs, color=DARK, zorder=4, linespacing=1.45)


def arrow(x1, y1, x2, y2, col=DARK, lw=1.3, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=13, lw=lw, color=col, zorder=2))


def drop(x, y, txt):
    """判废/剔除标注：向右伸出的灰色分支"""
    ax.add_patch(FancyArrowPatch((x, y), (x + 5.5, y), arrowstyle="-|>",
                                 mutation_scale=10, lw=1, color=GREY,
                                 linestyle="--", zorder=2))
    ax.text(x + 6.3, y, txt, ha="left", va="center", fontsize=7.3, color=GREY,
            linespacing=1.4)


def finding(y, txt):
    ax.text(99, y, txt, ha="right", va="center", fontsize=7.6, color=RED,
            style="italic", zorder=4)


ax.text(1, 98.5, "Fig S1  |  Analysis flow, with the quantities and the stopping rules "
                 "that produced them", fontsize=11, weight="bold", va="top")
ax.text(1, 95.2, "Grey branches are what each step removed. Red text on the right marks "
                 "where a methodological finding comes from.",
        fontsize=8.2, color="#666", va="top")

# ---------------- 暴露侧
box(2, 84, 34, 8.5, "Exposure — dynamic CD4⁺ T cell cis-eQTL",
    "8 activation profiles (naive / memory × 0 h, 16 h, 40 h, 5 d)\n"
    "85–100 donors per profile · GRCh38, no liftover", ec=BLUE, tcol=BLUE)

box(2, 72.5, 34, 8.5, "Instruments (strict)",
    "top cis-eQTL per gene × profile, P < 5×10⁻⁸\n"
    "3,556 exposures / 2,141 variants · min F = 22.2", ec=BLUE, tcol=BLUE)
arrow(19, 84, 19, 81)
drop(36, 76.5, "F > 10 filter never binds\n(min F = 22.2)")

# ---------------- 结局侧
box(62, 84, 36, 8.5, "Outcome — melanoma GWAS",
    "FinnGen R12 5,753 cases  +  Rashkin 6,777 cases\n"
    "meta: 12,530 cases / 789,099 controls", ec=RED, tcol=RED)
box(62, 72.5, 36, 8.5, "Also analysed as outcome",
    "FinnGen R8–R12 (five nested power levels)\n"
    "five other cancers · naevi (negative control)", ec=RED, tcol=RED)
arrow(80, 84, 80, 81)

# ---------------- MR
box(24, 60, 52, 8.5, "Mendelian randomization — Wald ratio",
    "one instrument per exposure  ⇒  z = β_out ⁄ se_out\n"
    "the P value is supplied entirely by the outcome GWAS", ec=DARK)
arrow(19, 72.5, 35, 68.5)
arrow(80, 72.5, 65, 68.5)
finding(64.2, "← the constraint that\norganises the paper")

# ---------------- 筛选层
box(24, 47, 52, 8.5, "Colocalisation → SMR/HEIDI → multi-instrument sensitivity",
    "coloc.abf, and PP.H3+PP.H4 > 0.5 added to the usual ratio criterion\n"
    "IVW · IVW-MRE · weighted median", ec=DARK)
arrow(50, 60, 50, 55.5)
drop(76, 51.5, "weighted mode dropped:\n1 significant of 226")
finding(46.5, "② HEIDI passes 253/291\nrecords coloc assigns to\ndistinct causal variants")

# ---------------- 名单
box(24, 34.5, 52, 8.5, "Candidate list",
    "FinnGen round: 4 novel-locus genes    meta round: 6\n"
    "shared between the two rounds: none", ec=RED, tcol=RED, fc="#FDF6F5")
arrow(50, 47, 50, 43)
finding(38.7, "①④ where the signal is,\nand how little it survives")

# ---------------- 功能层
box(2, 18.5, 30, 12, "Function — CD4-resolved",
    "single-cell ICB cohorts\nlineage-purity control at the CELL level\n"
    "locked 16-gene signature, and without TPI1", ec=GREEN, tcol=GREEN)
box(35, 18.5, 30, 12, "Compartment attribution",
    "cell-type-resolved atlases (two cohorts)\nspatial · TCGA bulk survival\n"
    "positive controls: MLANA, PTPRC, HLA-C", ec=GREEN, tcol=GREEN)
box(68, 18.5, 30, 12, "Mechanism — in vitro",
    "purified CD4 multiome, 40,495 nuclei\nactivation-invariant peaks only\n"
    "motif enrichment by family", ec=GREEN, tcol=GREEN)
arrow(40, 34.5, 20, 30.5)
arrow(50, 34.5, 50, 30.5)
arrow(60, 34.5, 80, 30.5)

# ---------------- 判废
box(14, 4.2, 72, 11.5, "Tests discarded by their own controls — reported, not interpreted",
    "within-lymphoid spatial subset: positive control also failed  ·  "
    "pre-specified immune modules: genes not measurable\n"
    "peripheral-blood biomarker analysis: precondition test failed, analysis cancelled  ·  "
    "probe-based Visium: target genes absent from the panel\n"
    "chromatin test of the activation window: mismatched to a genotype × timepoint claim",
    ec=GREY, tcol=GREY, fc="#F7F7F8", fs=7.6)
for x in (20, 50, 80):
    arrow(x, 18.5, x, 15.7, col=GREY, lw=1)

ax.text(50, 1.4, "Every quantity above is a count from the analysis, not an illustration. "
                 "Stopping rules were fixed before the step they govern.",
        ha="center", fontsize=7.8, color="#666")

fig.savefig(os.path.join(OUT, "FigS1_flow.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "FigS1_flow.png"), dpi=300, bbox_inches="tight")
print("FigS1_flow ok")
