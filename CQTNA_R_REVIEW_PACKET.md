# cqtna R 包 —— 第三方评审材料

**日期**：2026-08-17 · **版本**：0.2.0 · **状态**：`R CMD check` OK，106 项断言全过
**仓库根**：`D:/R_ex/MR/`（本文所有路径相对于此）

本文件给外部评审者：**它是什么、文件在哪、为什么这样设计、
以及应该往哪里使劲挑毛病**。最后一节列了我自己认为最薄弱的地方，
请优先打那里，不要浪费时间在我已经知道的问题上。

---

## 一、这是什么，从哪来

一个审计工具，输入是**别人已经做完的** cis-eQTL 工具化 MR 候选名单，
输出是"每个基因究竟拿到了哪一类证据"。

它不是从零设计的：

1. 上游是一篇方法学审计论文（CD4⁺T 动态 eQTL → 黑色素瘤靶点提名），
   其 Discussion 列出**八条检查**，说"每一条被推翻的主张，都是被一条
   同类研究做得起、却通常不做的检查推翻的"。
2. 该论文已附带一个 **Python 实现 `cqtna/cqtna.py`**（595 行，仅 numpy+pandas），
   把八条里**能自动化的五条**做成六个模块，
   **不能自动化的三条点名列出、每次报告都印**。
3. 本 R 包是把它移植到 R，因为目标用户（TwoSampleMR / coloc 生态）在 R 里，
   并新增一个模块 G（窗口敏感性）。

**Python 版仍然保留，作为 R 版的参照实现（oracle）。** 见第五节。

---

## 二、文件路径与规模

### R 包本体：`cqtna_r/`

| 路径 | 行 | 内容 |
|---|---|---|
| `DESCRIPTION` | 26 | 元数据。依赖只有 stats / utils / graphics |
| `NAMESPACE` | 27 | 手写，13 个导出函数 + 4 个 S3 方法 |
| `R/utils.R` | 112 | 统计与区间原语：BH、Simes、单侧 Fisher、位点聚类、已知位点索引 |
| `R/io.R` | 123 | 输入强制转换 + **基因组构建版本守卫** |
| `R/modules.R` | 264 | **模块 A–G**，各自可独立调用 |
| `R/audit.R` | 128 | 总入口 cqtna_audit、证据分层、三条不可自动化诊断 |
| `R/report.R` | 244 | print / summary / plot / cqtna_report 渲染 |
| `R/cqtna-package.R` | 11 | 包级文档 |
| **R 代码合计** | **882** | |
| `tests/testthat/test-paper-numbers.R` | 104 | 论文自身数字的回归测试 |
| `tests/testthat/test-python-oracle.R` | 97 | 与 Python 参照实现比对 |
| `tests/testthat/test-guards.R` | 55 | R 版新增守卫 |
| **测试合计** | **260** | 运行时 **106 项断言** |
| `man/*.Rd` | 410 | roxygen 生成，15 个 |
| `README.md` | 104 | |
| `inst/extdata/*.tsv` | — | demo 夹具，见下 |
| `.github/workflows/R-CMD-check.yaml` | 32 | 4 平台 CI 矩阵 |
| `CITATION.cff` / `.zenodo.json` | 21 / 8 | 发布元数据（OWNER 是占位符） |

### 夹具：`cqtna_r/inst/extdata/`

| 文件 | 大小 | 是什么 |
|---|---|---|
| `mr_meta.tsv` | 274 KB | 3,556 条 MR 记录（meta 结局，12,530 例黑色素瘤） |
| `mr_finngen.tsv` | 273 KB | 同一暴露数据、换 FinnGen 结局（模块 C 用） |
| `known_melanoma.tsv` | 3.5 KB | 157 个已知黑色素瘤 lead SNP |
| `known_hcc_mismatched.tsv` | 1.7 KB | 73 个 HCC lead SNP，作**错配阴性对照** |
| `instruments.tsv` | 0.4 KB | 28 个糖酵解基因的三级可用性 |
| `expression.tsv` | 0.1 KB | TPI1 的分细胞类型表达 |
| `reference_python.tsv` | 2.2 KB | **Python 参照实现的输出**，84 个键值对 |

### Python 参照实现：`cqtna/`

| 路径 | 行 | |
|---|---|---|
| `cqtna.py` | 595 | 参照实现，`python cqtna.py --demo` 可跑 |
| `README.md` | 96 | |
| `cqtna_manual.md` | 84 | **三条不可自动化诊断的失败模式**，写得比 R 包好，值得单读 |
| `demo/` `demo_out/` | — | demo 输入与已提交的输出（确定性，逐字节可复现） |

### 设计文档

| 路径 | 内容 |
|---|---|
| `CQTNA_R_PACKAGE_PLAN.md` | 动工前的方案构思 + 文末第十一节建成记录 |
| `manuscript/SUPP_window_sensitivity.md` | 模块 G 的方法学依据（S36） |

---

## 三、数据流与 API

用户的 MR 结果（record_id, gene, chr, pos, p）经 `as_cqtna_mr()`：
重算 FDR、按 1 Mb 单连锁聚类成独立位点、记录基因组构建版本。
已知位点名单（chr, pos[, source]）经 `as_cqtna_known()`：缺 source 列警告。
两者交给 `cqtna_audit()` 跑模块 A–G，再经 print / summary / plot / cqtna_report 输出。

可选输入：错配名单（→ A、G 的阴性对照）、第二个结局的 MR 结果（→ C）、
通路基因表（→ D）、表达表（→ E）、峰位置表（→ F）。

**核心口径：凡是计数要承担论证的地方，一律按独立位点算，不按基因记录算**
——因为一个显著位点并不指名一个基因。

---

## 四、七个模块各算什么

| | 函数 | 统计定义 |
|---|---|---|
| **A** | `cqtna_attribution()` | 位点按染色体内 1 Mb 单连锁聚类；某位点"已知" = 其上任一变异距最近已知 lead SNP ≤ known_kb。倍数 = (显著位点中已知比例)/(背景位点中已知比例)，P = **单侧 Fisher 上尾**（phyper 上尾）。换一份别的疾病名单再跑一次 = 错配阴性对照 |
| **B** | `cqtna_unit_sweep()` | 同一份结果在 record / variant / gene / locus 四种推断单位下重算：用 **Simes** 合并同组 p，再 BH |
| **C** | `cqtna_stability()` | 换第二个结局 GWAS 后，基因层与位点层的 Jaccard，以及丢失/新增名单 |
| **D** | `cqtna_ladder()` | 三级计数：可工具化 → 与本结局 harmonise 后仍可分析 → 名义显著 |
| **E** | `cqtna_compartment()` | 提名基因"其他细胞类型的最高表达 / 目标细胞类型表达"之比 |
| **F** | `cqtna_peak_distance()` | eQTL 峰与 GWAS 峰的碱基距离 |
| **G** | `cqtna_window_sweep()` | **两个都叫"1 Mb"的约定分开扫**（100/250/500/1000 kb）：known_kb 直接动倍数，locus_kb 动分母。外加**连续距离分布**，并判断所选阈值落在数据空隙里还是踩在数据上 |

**模块 G 的论证**：为迎合结果而挑的阈值，收紧时会变弱。
若倍数随窗口收紧而**上升**，则你报的是保守值。
demo 数据上：4.09 倍 → 18.26 倍（收紧到 100 kb），且 1 Mb 落在 74 kb 与 2,167 kb 之间的空隙里。

---

## 五、验证策略（最该被评审的部分）

三套测试，**106 项断言，0 失败，2 条预期警告**：

### 1. test-paper-numbers.R —— 回归测试就是论文自己的结果

直接断言已发表数字：4.09 倍 / 单侧 P = 0.0281（12 位精度）、错配 0.00 倍 / P = 1、
record 层 10 基因 vs locus 层 28 基因、Jaccard 0.455、28 → 3 → 2 → 1、
两个窗口扫描的全部 8 个倍数、7 个连续距离。

### 2. test-python-oracle.R —— 跨实现比对

`inst/extdata/reference_python.tsv` 由 `python cqtna.py --demo` 导出并随包提交，
覆盖模块 A/B/C/D/G 的 84 个值。**测试环境不需要装 Python。**
任一侧实现漂移，测试立刻失败。

> 这是当初选"纯 R 移植 + Python 留作 oracle"而非"Python 退休"的全部理由。
> 它确实抓到了一处真实分歧，见第六节。

### 3. test-guards.R —— R 版新增行为

构建版本必须声明、不一致报错；缺 source 列警告；缺错配名单警告；
TwoSampleMR 风格列名可接受；**FDR 一律重算**（有一条测试往输入里塞
fdr = 0 投毒，验证不被采信）。

### 独立复现方式

```
cd cqtna && python cqtna.py --demo
R CMD build cqtna_r && R CMD check cqtna_0.2.0.tar.gz --no-manual
Rscript -e 'testthat::test_local("cqtna_r")'
```

实测：3,556 条记录跑完整 audit 用 **0.31 秒**；3.5 万条强制转换 0.12 秒。

---

## 六、oracle 抓到的一处真实分歧（已处置，请复核处置是否正确）

位点扫描在 100 kb 处的真值是 **5.325000000000000177**。

- R 的 `round(v, 2)` 给 **5.32**（R 会向十进制字面量方向"修正"）
- `sprintf("%.2f", v)` 与 Python 的 `round()` 给 **5.33**（按实际二进制值取整）

真值高于中点，我判断 **5.33 正确**，也是论文与 Python 参照所报。

**处置**：模块返回值一律**全精度**，四舍五入只在渲染层用 sprintf 做。
另修 `cq_fmt_kb`：formatC(format="d") 是截断，会把 74 kb 显示成 73。

⚠ **请评审这个判断**：也可以主张"两种取整都合法，包不该在 0.01 上较真"。
我选择与论文和参照实现一致，但这是判断，不是定理。

---

## 七、刻意保留的三条"拒绝"（不是遗漏）

1. **没有任何路径能走到 target-supported。** 报告只发
   screened / unresolved / state-informative。
   能授予那个词的检查，恰好是它跑不了的那三条。
2. **三条不可自动化的诊断每次都印**，且同时挂在返回对象的
   not_automated 字段上——因为 R 用户会 `as.data.frame()` 然后把警告丢掉。
3. **缺错配阴性对照时警告而非静默跳过**：没有它就分不清
   "结局特异的归属"与"在所有疾病里都密集的位点"。

---

## 八、★ 请优先打这些地方（我自己认为最薄弱的，已排序）

1. **模块 G 的 in_gap 判据是拍脑袋的。**
   `R/modules.R` 里是 `nearest_above >= 4 * nearest_below`，
   这个 **4 倍没有任何理论依据**，纯粹是在 demo 数据上看着合理，
   却直接决定报告里印"落在空隙里"还是"踩在数据上"。
   → 该不该有这个自动判语？还是只报距离、让人自己看？

2. **证据分层里 `ratio > 2` 判 state-informative 同样是拍脑袋的**
   （`R/audit.R` 的 cq_tiers）。

3. **模块 A 的 Fisher 检验把位点当独立检验，而位点是聚类出来的。**
   聚类阈值一变位点数就变（100 kb 下 1,065 个 vs 1 Mb 下 554 个），P 值随之变。
   模块 G 报了这个敏感性，但**没有解决**独立性假设本身。
   → 有没有更该用的检验？置换？

4. **模块 B 用 Simes 合并同位点的 p，再对合并后的 p 做 BH。**
   Simes 在正相依下有效，同位点记录确实正相依，
   但"Simes 之后再 BH"这一层的 FDR 控制**没有证明**。
   包里印了警告说"最短的名单不等于 FDR 被控制住了"，那是免责不是解决。

5. **错配名单是弱阴性对照，论文自己也这么说。**
   它能排除"在所有疾病里都密集的位点"，但排除不了"本身就是疾病特异的密度"。
   后者需要按工具变量密度与基因密度做置换，本工具**不做**。

6. **夹具全部来自同一篇论文的同一个分析。** 没有外部数据集验证。
   一个只在自己数据上验证过的审计工具，泛化性是未知的。
   → 建议评审者拿自己的一份 MR 结果跑一遍。

7. **未测的边界**：X/Y 与非常染色体的处理、单一位点、全部显著、
   百万级记录的性能（已测 3.5 万条，再往上未测）。
   零显著记录已测：返回 NA 而非报错。

8. **`cqtna_unit_sweep()` 用 `interaction(..., sep="\r")` 拼分组键**，
   若基因名含 `\r` 会串组。实际不会发生，但是隐患。

9. **模块 E / F 在 demo 里是空的或未跑**（demo 没配 peak 表）。
   新用户跑 demo 会看到两个模块 "configured but empty"，容易以为坏了。
   发布前该配小样例。

---

## 九、还没做的（需要账号，我没动）

1. 建 GitHub 仓库并 push；DESCRIPTION / README.md / CITATION.cff /
   .zenodo.json 里的 **OWNER 是占位符**，建库后统一替换。
2. Zenodo 打开仓库开关 → GitHub 打 release（建议 v0.2.0）→ 自动铸 DOI。
3. 正文 code availability 需提一句 R 包（目前只写了 CQTNA）。
4. 是否投 CRAN 未定：R CMD check 已 OK，技术上够格，
   但 CRAN 带来长期维护义务，而 GitHub + Zenodo 已满足可引用性。

---

## 十、给评审者的最短上手路径

```
# 1. 看设计立场（10 分钟）
cat cqtna_r/README.md
cat cqtna/cqtna_manual.md

# 2. 看核心算法（30 分钟）
cat cqtna_r/R/utils.R
cat cqtna_r/R/modules.R

# 3. 看它拿什么证明自己（20 分钟）
cat cqtna_r/tests/testthat/test-paper-numbers.R
cat cqtna_r/tests/testthat/test-python-oracle.R

# 4. 自己跑一遍
R CMD check cqtna_0.2.0.tar.gz --no-manual
Rscript -e 'testthat::test_local("cqtna_r")'
```

---

最后一句照抄 Python 版 README 的收尾，它同样适用于本包：
**一个用来检测选择性强调的工具，应该先拿它自己去审提出它的那项研究。**
