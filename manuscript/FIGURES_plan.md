# 图件方案（定稿编号）

日期：2026-08-11
原则：图序 = 论证序（①→②→③④→收口→归属方法→应用示例→证据矩阵），不按数据类型分组。

## 主图

| 定稿 | 现文件 | 生成脚本 | 承担的论点 | 状态 |
|---|---|---|---|---|
| **Fig 1** | `FigA_locus_attribution` | `figures/make_figures.py::fig_A` | ① 信号从哪来：Manhattan，按位点类别着色 | 已有 |
| **Fig 2** | `FigB_coloc_vs_HEIDI` | `figures/make_figures.py::fig_B` | ② HEIDI 放行 coloc 判为 LD 混淆的信号（253/291） | 已有 |
| **Fig 3** | `FigC_power_dependence` | `figures/make_figures.py::fig_C` | ④ coloc 后验随结局数据集移动，误差双向（PARP1 假阳性 / ZFYVE19 假阴性）。⚠ **不得写成"由功效决定"**——两轮同时换了功效与队列，多信号分析只判得了 PARP1（正文 §2.4）| 已有 |
| **Fig 4** | `FigF_power_stability` | `figures/make_fig_power_stability.py` | ★★ 收口：功效-稳定性曲线 + 按位点类别分层的恢复率 | **本次新增** |
| **Fig 5** | `FigE_spatial_compartment` | `figures/make_fig_spatial.py` | ⑤ 区室归属（含 HLA-C 镜像阳性对照） | 已有 |
| **Fig 6** | `Fig6_glycolysis_TPI1` | `figures/make_fig6_glycolysis.py` | 应用示例：16h 窗口 + 通路级应答分层 | 已有 |
| **Fig 7** | `Fig7_crosscancer` | `step30e_locus_level_and_fig.py` | 跨癌种：按独立位点的富集（a）、TPI1 特异 vs SMC2 泛癌（b）、命中数随结局功效（c） | 已有 |
| **Fig 8** | `FigD_evidence_matrix` | `figures/make_figures.py::fig_D` | 候选证据矩阵：无一候选三维度全优 | 已有 |

⚠ 文件名保持不变（脚本互相引用），仅在正文与图注中使用 Fig 1–8 的编号。
投稿前统一改名时，改 `save()` 的传参即可，勿手动重命名 PDF。

## Fig 4 三个面板（新增图）

- **a** 黑色素瘤降功效曲线：命中数（左轴）+ 与全功效名单的 Jaccard（右轴，含区间）。
  五角星 = 独立测得的 FinnGen 轮真实值，落在模拟区间内 → **外部校准，非自证**。
- **b** 五癌种同一曲线（横轴为各自观测病例数的比例），50% 恢复线。
- **c** ★★ 按位点类别分层：已知色素/痣位点 vs 新位点。灰带 = 同类研究常见的 3,000–8,000 例区间。

数据源：`53a`、`53b`、`54a`、`54b`、`55a`。曲线**只画到观测功效为止**（两次外推尝试均被 sanity check 判废，见 FINDINGS Step 53），图注须写明不可外推。

## 两处待统一的数字（写正文前须定）

1. **Fig 4a 的校准点**：`53b_calibration.tsv` 给出 FinnGen 功效下模拟命中 **10.2** vs 真实 10（比值 1.02）。
   HANDOFF v2 第七节记的 "9.6（比值 0.96）" 实为 `53a` 中 **5,000 例**那一行，不是校准点本身。
   两个数都对，但**校准应引 10.2**；HANDOFF 该行建议改。
2. **发现①报几倍**：按记录计数为 12/21 = 57.1%、5.26×、p=2.7e-7；
   按独立位点计数为 3/7 = 42.9%、**4.09×、p=0.028**（`36e`）。
   Step 30 的自设纪律是"应改为按独立位点计数"，故**正文与 Fig 1/Fig 7a 一律报 4.09×**，
   记录级数字仅在补充材料中给出。Discussion 草稿 §3① 目前用的是记录级，需改。

## 补充图（Supplementary，暂列）

| 编号 | 内容 | 来源 |
|---|---|---|
| S1 | 工具变量与 harmonise 流程图 + F 统计量分布 | `01`/`04` |
| S2 | Steiger 方向性与有界 R²（`rsq.exposure.bounded`） | `09` |
| S3 | 敏感性分析四法对比（含 weighted mode 不适用的说明） | `18`/`18b` |
| S4 | 痣（阴性对照表型）与共病 MR | `20` |
| S5 | 纯度对照：谱系匹配前后效应量塌缩 | `29a`–`29h` |
| S6 | 二次归一化错误的诊断与重算前后对照 | `31`–`32h` |
| S7 | GSE282266 轴的无偏刻画 + motif（全部 peak vs 活化不变 peak） | `44a`–`44c`、`50a`–`50b` |
| S8 | GSE199994 ICC 与随机模块地板（方向关闭的依据） | `52a`–`52d` |
