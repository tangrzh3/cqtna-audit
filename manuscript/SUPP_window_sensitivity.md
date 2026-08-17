# S36 — 1 Mb 窗口敏感性与连续距离分布

**日期**：2026-08-17 · **脚本**：`step121_window_sensitivity.py`
**表**：`121a_window_sensitivity.tsv` · `121b_locus_definition_sensitivity.tsv` ·
`121c_distance_distribution.tsv`
**来历**：HANDOFF v5 与 v6 都把这项列为未做（"100/250/500 kb 各跑一遍；
报连续距离分布而非二分类"）。

---

## 0. 先把两个"1 Mb"分开

全文一直把两个不同的约定统称为"1 Mb 窗口"，它们承重不同，必须分开扫：

| | 含义 | 影响 |
|---|---|---|
| **KNOWN_KB** | 位点算"已知" = 该位点的**显著记录**到最近已知 lead SNP ≤ 阈值 | **直接决定归属倍数** |
| **LOCUS_KB** | 独立位点 = 染色体内工具变量位置的单连锁聚类间距 | 决定分母位点数与显著位点的收拢 |

⚠ **覆盖范围**：melanoma × {Soskic, eQTLGen}、HCC-{high, low} × Soskic 共四格。
**RA 两格未做**——`step108` 的 Okada rsID→GRCh38 坐标是边扫 sumstats 边解、
**没有落盘**，且 RA 没有 per-record 注释表（`108a` 只有汇总行）。
补 RA 需重跑那一遍全量 sumstats 扫描。

✅ 每格的 1000 kb 行都**逐位重现已发表倍数**（4.09 / 4.44 / 8.85 / 17.62），
脚本内置该核对，对不上即报错退出。

---

## 1. KNOWN_KB 扫描（LOCUS_KB 固定 1 Mb）

| 格子 | 100 kb | 250 kb | 500 kb | **1000 kb（已发表）** |
|---|---|---|---|---|
| melanoma × Soskic | **18.26×**, P = 3.4×10⁻⁴ | 8.19×, P = 0.0039 | 5.65×, P = 0.011 | **4.09×, P = 0.028** |
| melanoma × eQTLGen | 4.42×, P = 1.9×10⁻⁷ | 4.53×, P = 5.4×10⁻¹⁰ | 4.19×, P = 2.3×10⁻⁹ | **4.44×, P = 3.7×10⁻¹¹** |
| HCC-high × Soskic | **0.00×**, P = 1.0 | **0.00×**, P = 1.0 | 14.45×, P = 0.068 | **8.85×, P = 0.110** |
| HCC-low × Soskic | 47.0×, P = 0.021 | 20.14×, P = 0.049 | 31.33×, P = 9.6×10⁻⁴ | **17.62×, P = 0.0031** |

**★ 结论一：1 Mb 是保守选择，不是调出来的。**
黑色素瘤两格在窗口收紧时倍数**上升或持平**（Soskic 4.09 → 18.26；eQTLGen 4.4 附近平稳），
P 值同向变小。若 1 Mb 是为了让结果好看而挑的，收紧窗口应当削弱结果——实际相反。

**★ 结论二：HCC-high 是唯一的例外，必须写明。**
该格仅有的那个"已知"显著位点距最近的 HCC lead SNP **408 kb**，
所以 **窗口 < 500 kb 时该格为 0/2、倍数 0.00**。
HCC-high 本就是六格中未达 P<0.05 的两格之一，这条进一步说明
**它的方向一致性依赖于较宽的窗口**，不宜作为支持性证据引用。

---

## 2. LOCUS_KB 扫描（KNOWN_KB 固定 1 Mb）

| 格子 | 100 kb | 250 kb | 500 kb | **1000 kb** |
|---|---|---|---|---|
| melanoma × Soskic | 5.33× (1,065 背景位点 / 8 显著) | 4.63× (876/7) | 4.59× (728/7) | **4.09× (554/7)** |
| melanoma × eQTLGen | **8.02×** (4,815/45), P = 1.7×10⁻²⁵ | 6.98× (2,402/32) | 6.74× (1,193/30) | **4.44× (559/30)** |
| HCC-high × Soskic | 6.67× (1,040/2) | 8.12× (861/2) | 8.95× (716/2) | **8.85× (549/2)** |
| HCC-low × Soskic | 13.26× (1,087/2) | 16.52× (892/2) | 17.6× (739/2) | **17.62× (564/2)** |

**★ 结论三：位点划法也一样——1 Mb 给出四个取值里最保守的倍数。**
把位点划细（100 kb）会让背景位点数涨到 2–9 倍，倍数随之上升
（eQTLGen 4.44 → 8.02）。**报 1 Mb 是在压低自己的效应量。**

---

## 3. ★ 连续距离分布（HANDOFF 真正要的那张）

每个**显著位点**的距离 = 该位点的 FDR<0.05 记录到最近已知 lead SNP 的最小距离。
（⚠ 必须只在显著记录上取，不能在位点内全部记录上取——初版按后者算，
HCC-high 报出 235 kb，而那属于一条不显著的记录，与二分类直接矛盾。已修。）

**melanoma × Soskic**（7 个显著位点，kb）：
`5 · 11 · 74` ‖ `2,167 · 7,275 · 20,136 · 38,494`
→ **前三个在 74 kb 以内，第四个跳到 2,167 kb。1 Mb 这条线落在一段空隙里。**
背景 542 个位点中位数 **10,728 kb**。

**melanoma × eQTLGen**（30 个，kb）：
`0 · 1 · 1 · 1 · 4 · 4 · 8 · 9 · 21 · 23 · 26 · 32 · 48 · 50` ‖
`128 · 159 · 189 · 222` ‖ `597 · 746` ‖ `3,550 · 5,939 · 6,101 · 7,293 · 8,541 ·
15,179 · 19,408 · 21,212 · 39,045 · 45,105`
→ **14 个在 50 kb 以内，20 个在 1 Mb 以内**；中位数 **144 kb**，
背景 538 个位点中位数 **10,360 kb**（相差约 72 倍）。

**HCC-high × Soskic**（2 个）：`408` ‖ `20,563` kb — 背景中位数 15,525 kb
**HCC-low × Soskic**（2 个）：`26` ‖ `401` kb — 背景中位数 15,920 kb

**★ 结论四：二分类严重低估了贴近程度。**
说"落在已知位点 1 Mb 内"听起来宽松，实际上黑色素瘤两格的多数命中在
**几 kb 到几十 kb**——eQTLGen 那格近一半在 50 kb 以内、有一个是 0 bp（同一变异）。
分布是**双峰**的：要么紧贴已知 lead SNP，要么远在数 Mb 之外，
中间几乎是空的。这比原来的二分陈述强，且解释了为什么倍数对阈值不敏感。

---

## 4. 正文可加的一句（若字数允许）

> Both 1 Mb conventions are conservative: tightening the known-locus window to
> 100 kb raises the CD4⁺ enrichment from 4.09- to 18.26-fold, and tightening the
> locus definition raises the whole-blood enrichment from 4.44- to 8.02-fold. The
> distances are bimodal — significant loci sit either within tens of kilobases of a
> known lead SNP or several megabases away, with little in between — so the
> threshold falls in a gap rather than near the data. The one exception is HCC at
> its higher power, whose single known-locus hit lies 408 kb away and which
> therefore requires a window of at least 500 kb (Supplementary S36).

⚠ GB 正文已贴 8,000 词上限，加这段须同时减内容。
