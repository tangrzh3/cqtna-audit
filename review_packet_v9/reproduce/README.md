# `reproduce/` — 从中间表重算全部位点级结果

**加入原因**：v9 第一版复核包声称"随包提供 per-record 中间表，可重建全部位点结果"，
但**实际没有附上**，复核者只能核对最终 TSV。这是包装疏漏，已补齐。

---

## 怎么跑

```r
install.packages("../cqtna_0.3.0.tar.gz", repos = NULL, type = "source")
setwd("<解压目录>/reproduce")
```

脚本里的 `MR <- "D:/R_ex/MR"` 与 `setwd(MR)` 需改成你的解压路径，
或直接把本目录内容放到一个工作目录里再跑。

| 跑什么 | 得到 | 大约耗时 |
|---|---|---|
| `Rscript step124_full_grid.R` | `123d` 主网格（六格 + 两功效敏感性 × 三种口径）| < 1 分钟 |
| `Rscript step125_mismatch_loci.R` | `125a/b` 逐位点错配诊断 | < 1 分钟 |
| `Rscript step128_window_sensitivity_fixed_anchor.R` | `128a/b/c` 窗口扫描与距离分布 | 约 2 分钟 |
| `Rscript step130_multilist_control.R` | `130a/b/c` 多名单对照 | 约 3 分钟 |
| `Rscript step126_recompute_on_fixed_anchor.R` | `126a/b/c` 非网格格子 + 置换（10,000 次 × 8 格 × 5 容差）| **约 10 分钟** |
| `python step127_audit_manuscript_numbers.py` | 正文 ↔ 表格对账（需 `../MANUSCRIPT_GB.md`）| 秒级 |

`step129_selfcheck_on_fixed_anchor.R` 需要 `96a`/`92c` 之外的历史表，**未附**；
它只回答"自检表是否随分区变"，不产生正文数字。

---

## 附了什么

**per-record 中间表**（列口径 `cell, record_id, gene, chr, pos, p`，可直接喂 `as_cqtna_mr()`）

| 文件 | 内容 |
|---|---|
| `123b_ra_records.tsv.gz` | RA 五个变体的逐记录表（35,795 条）|
| `123c_eqtlgen_hcc_records.tsv.gz` | eQTLGen × HCC 两格（25,059 条）|
| `99d_r13_records.tsv.gz` | R13 迁移检验九张表 |

**已知位点参照系**：`landi2020_…`（melanoma 157）· `84a_…`（HCC 73）·
`known_loci_colorectal_grch38.csv`（642）· `known_loci_prostate_grch38.csv`（35）·
`123b_ra_known_positions.tsv`（Okada 定位出的 83 个 GRCh38 坐标）

**逐记录注释表**：`13_meta_locus_annotation.tsv`（CD4 × melanoma）·
`92c_locus_annotated.tsv`（eQTLGen × melanoma）· `85a_HCC_{high,low}_annotated.tsv`

**脚本**：`step124_cells.R`（共享格子定义，其余脚本 `source()` 它）+ 上表六个 step。

---

## 仍然没附的

- **FinnGen / eQTLGen 原始 sumstats（约 1.5 GB）**。
  它们只被 `step99` / `step108` / `step94d` 用来**生成**上面那三张 per-record 表；
  这三个脚本因此也没附。
  ⚠ 也就是说：**从 per-record 表往后的每一步都可复算，
  但"per-record 表本身是否从 sumstats 正确抽出"无法在本包内验证。**
  若要核这一层，需要原始文件与那三个脚本，可另行索取。
- `step85e_matched_background.py`（强度/频率匹配背景），依赖 `85a` 之外的列。
- 单细胞、coloc/HEIDI、患者分层部分——不在本次复核范围。

---

## 一致性自检

跑完后，下面几个数应当与包内 TSV 完全一致（浮点逐位）：

| 来源 | 应得 |
|---|---|
| `123d` melanoma × CD4，`analysis_status = prereg_primary` | 4/8，4.96×，P = 0.00484 |
| `123d` RA × eQTLGen 主行 | 错配 2.05×，P = 2.9×10⁻⁴，`control = FAILED` |
| `123d` legacy 行 RA × Soskic_CD4 | 3.46×，P = 1.65×10⁻⁵（复现 S33 §8 登记值）|
| `130c` RA × Soskic_CD4 | `F_own` 3.91，`F_max_all` 3.05（HCC），`margin_all` 1.28 |
| `128a` RA × CD4，`known_kb = 100` | 错配 4.17×，`control = FAILED` |

若有任何一项对不上，**先看 R 版本与 `cqtna` 是否为随包的 0.3.0**，
再看工作目录是否包含全部输入表。
