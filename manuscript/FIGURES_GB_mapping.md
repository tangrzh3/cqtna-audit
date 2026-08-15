# GB 九图 ↔ 图件文件 对照与缺口审计

日期：2026-08-15（HANDOFF_v6 §三-A1 的落实）
适用：`manuscript/MANUSCRIPT_GB.md` §Figures 所列的 **Fig 1–9**

> `FIGURES_plan.md` 记的是全文源 `MANUSCRIPT_v2_dual_thread.md` 的**八图**方案，
> 与 GB 的九图**不是同一套编号，也不是同一套分组**。
> GB 重排后，若干面板从原来的 `Fig9_part2` / `FigS7` 里被拆出来独立成图。
> **两份方案都有效，不要互相覆盖**；本文件只管 GB 这一套。

---

## 一、对照表

| GB | 图注要求 | 现有资产 | 数据源 | 状态 |
|---|---|---|---|---|
| **1** | 两个结局下的位点归属，**按独立位点** | `FigA_locus_attribution`<br>`make_figures.py::fig_A` | `13_meta_locus_annotation.tsv` | ⚠ 只画了 **meta 一个结局**，且是**按记录**的 Manhattan。GB 要"两个结局 + 按独立位点"。缺 FinnGen 轮面板（`06_locus_annotation.tsv`）与位点级计数（`36e_locus_level_attribution.tsv`，4.09×） |
| **2** | 泛化：五个嵌套 release · 迁移检验 · 第二疾病 · 非癌结局 · 第二暴露资源 · 双轴交叉网格 + 错配位点对照 | `FigG_release_trajectory` 只覆盖第一项 | 见 §二 | ❌ **须新建合成图**。⚠ 建前须先跑 `step119`（§三） |
| **3** | 按位点类别的功效-稳定性，**以及同一分层对全功效 \|z\|** | `FigF_power_stability`（a/b/c） | `53a` `53b` `54a` `54b` `55a` | ⚠ 缺第四面板：对全功效 \|z\| 的分层（`101c_stratified.tsv`） |
| **4** | 共定位 vs SMR/HEIDI，**含 in-sample 精细定位** | `FigB_coloc_vs_HEIDI` | `16_coloc_meta_results.tsv` `15_SMR_meta_results.tsv` | ⚠ 缺 in-sample 精细定位面板（`60a_susie_outcome_signals.tsv`：MC1R 三个 high-purity credible set，log₁₀BF 45.6/36.6/19.5） |
| **5** | 通路上三个层级的工具变量可得性 | — | `75a` `34a` `34b` | ✅ **已建**：`Fig5_instrument_ladder`（`make_gb_fig5_instrument_ladder.py`）。28 → 3 → 2 → 1，与正文逐字吻合 |
| **6** | TPI1 可工具化的活化窗口：**八个 profile 上效应量对精度** | `Fig9_part2` 面板 a、b | `27a_TPI1_eqtl_timecourse.tsv` | ⚠ 须从 `Fig9_part2` 拆出独立成图。注意 GB 要的是 **β vs se**，不是现 `Fig6_glycolysis_TPI1` 的"脉冲 + 非应答者通路表达" |
| **7** | CD4⁺ 代谢轴及其染色质签名 | `FigS7_axis_motif`（补充）+ `Fig9_part2` 面板 e | `44a`–`44c` `50a`–`50b` `48b` | ⚠ 须从补充**提升为正文主图**并与面板 e 合并 |
| **8** | 区室归属 | `Fig5_compartment`（a–e，五面板） | `22c` `23e` `26c` `69b` `70a` | ✅ 内容吻合，**只需改名**。（`FigE_spatial_compartment` 是更窄的纯空间版，不是这张） |
| **9** | 三个队列的患者 | `Fig9_part2` 面板 c | `74a`（发现）· `87b`（跨病种 HCC）· 同病种复制队列 | ⚠ 须拆出独立成图。数字已核：发现 Δ=+0.867 P=0.00055（`74a` Post/16/all）；跨病种 Δ=+0.8744 P=0.04286（`87b` H1/T/post/**nonTreg**/sig16）——正文的 +0.874 取的是 nonTreg 行，**不是 all 行**（all/sig16 = 0.8686） |

**小结**：真正"没有文件"的只有 Fig 2 与 Fig 5；但 Fig 1/3/4 各缺一个面板，
Fig 6/7/9 需要从现有合成图里拆分重组，Fig 8 只需改名。
HANDOFF_v6 的判断（"不是编号问题，是图件本身没按九图口径导出"）成立。

---

## 二、Fig 2 的六个组成部分与各自的权威数据源

| 组成 | 数字 | 表 |
|---|---|---|
| 五个嵌套 release | R8–R12，2,705→5,753 例；命中 7→10；已知位点恢复早且不动 | `59a_release_trajectory.tsv`（`own` 与 `common` 两口径）· `58a`（登记预测区间） |
| 迁移检验 R13 | 6,226 例，同六基因同两位点；归属 9.58×，P=0.0107 | `99c_r13_attribution.tsv` |
| 第二疾病 HCC | high 8.85×/P=0.110 · low 17.62×/P=0.0031 | `85c_G1_enrichment.tsv` |
| 非癌结局 RA | Soskic 3.459×/P=1.65e-5 · eQTLGen 3.566×/P=7.70e-14；去 MHC 后 3.026/3.319 | `108a_ra_attribution.tsv` |
| 第二暴露资源 eQTLGen | melanoma 4.44×/P=3.72e-11（vs Soskic 4.09×） | `92c` / `92d` / `92e`（匹配背景 4.14–4.48×，emp P=1e-4） |
| 交叉网格 + 错配对照 | 见下 | ⚠ **无单一权威表**，见 §三 |

### 主网格（三疾病 × 二资源 = 六格，四格 P<0.05）

| | Soskic CD4 | eQTLGen 全血 |
|---|---|---|
| melanoma | 4.09×，P = 0.028 ✔ | 4.44×，P = 3.7e-11 ✔ |
| HCC-high | 8.85×，P = 0.110 | 5.11×，P = 0.052 |
| RA | 3.46×，P = 1.7e-5 ✔ | 3.57×，P = 7.7e-14 ✔ |

功效敏感性（**不是第七、八格**）：HCC-low × Soskic 17.62×/P=0.0031 · HCC-low × eQTLGen 8.67×/P=0.017。

### 错配位点对照

| 组合 | 倍数 | P | 出处 |
|---|---|---|---|
| eQTLGen × melanoma 用 **HCC 名单** | 1.30× | 0.411 | ⚠ 仅见 `PREREG_generality_grid.md` §8 表，**无 TSV** |
| Soskic × HCC-low 用 **melanoma 名单** | 0.00× | 1.0 | ⚠ 同上，**无 TSV** |
| Soskic × melanoma R13 用 **HCC 名单** | 0.00× | 1.0 | `99c_r13_attribution.tsv` ✔ |
| Soskic × RA 用 **melanoma 名单** | 1.596× | 0.161 | `108a_ra_attribution.tsv` ✔ |

---

## 三、⚠ 画 Fig 2 之前必须先修的两处（新发现，不是审稿人提的）

### 1. `94d_grid.tsv` 里混着两行**已作废、且无任何作废标记**的结果

`PREREG_generality_grid.md` §8.2「偏离 1」写明：五个 FinnGen 癌种的已知位点坐标映射有 bug
（`step94c` 把已知位点 rsID 拿到**已被过滤成工具变量位置的**结局提取里查坐标，
致背景已知位点只有 2/559 与 5/559，而黑色素瘤为 58/554），
故 **lung（0/8，fold 0.00）与 colorectal（3/22，fold 15.25）"不得引用"**。

但这两行**原样躺在 `94d_grid.tsv` 里**，与有效行并列，没有 `void` 列、没有注释。
任何脚本（包括图件脚本）读这张表都会把它们当成有效格子。

**处置**：`step119` 重建时给每行加 `status` 列（`main` / `sensitivity` / `void_bad_coords`），
并在表头注明作废理由与出处。原表保留不动。

### 2. `94e_grid_summary.tsv` 是**旧网格**的汇总，与现行口径不符

现内容 `n_testable=6, n_fold_gt1=5, pct=83.3`。这 6 格是
Soskic×{melanoma, HCC-high, HCC-low} + eQTLGen×{melanoma, **lung**, **colorectal**}
——**含两行作废数据，且不含 RA**。而现行主网格是 melanoma/HCC-high/RA × 2 资源，
预注册 §8 记的 U1 是 **6/6 = 100%**，不是 5/6 = 83.3%。

**处置**：`step119` 重算并另存 `119b_grid_summary.tsv`；`94e` 保留 + 加日期批注。

### 3. 两个错配对照数字没有生成脚本

`1.30×/P=0.411` 与 `0.00×/P=1.0` 只存在于预注册 §8 的记录里，
仓库中**没有任何 TSV 承载**，也没有能重跑出它们的脚本
（`step94c` 里编码的 NC 是"melanoma 名单套胰腺癌"，因胰腺癌名单为空而被放弃，
改成两名单互错配的那次运行**没有落盘**）。

这与上一窗口发现的 `96a` 无生成脚本属同一类问题，且**可以廉价补上**：
两份名单都在仓库里（`landi2020_known_loci_grch38.csv` 157 个位点 ·
`84a_hcc_known_loci_grch38.csv` 73 个位点），逐格重打分即可。

**处置**：`step119` 重算两个错配对照，落 `119c_mismatch_controls.tsv`，
并与预注册 §8 记录的数字逐位比对；**若对不上，以重算为准并披露差异**。

---

## 四、执行顺序

1. `step119`：重建网格主表 + 汇总 + 错配对照（→ `119a`/`119b`/`119c`）
2. **Fig 2** 合成图（只读 `119a`–`119c`、`59a`/`58a`、`99c`、`108a`、`92d`）
3. Fig 1 补 FinnGen 轮面板 + 改为位点级
4. Fig 3 补 \|z\| 分层面板 · Fig 4 补 in-sample 精细定位面板
5. Fig 6 / 7 / 9 从 `Fig9_part2`、`FigS7` 拆分重组
6. 统一导出 `Fig1`–`Fig9`（PDF + PNG）

⚠ 文件名一律在 `save()` 传参处改，**不要手动重命名 PDF**（脚本之间互相引用）。
