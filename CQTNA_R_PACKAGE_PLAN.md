# 把 CQTNA 做成 R 包 —— 方案构思

**日期**：2026-08-17 · **状态**：构思，未动工
**前提**：`cqtna/` 已有一个可跑的 Python 实现（469 行，仅依赖 numpy + pandas），
正文 code availability 已写明"CQTNA, a runnable implementation of the diagnostics
in the Discussion, is included in the deposit"。

---

## 一、先把问题问对

**不是"能不能做一个包"，而是"要不要给已有的 CQTNA 换一个 R 前门"。**

CQTNA 已经把最难的三件事做完了，这三件事恰恰是包做不好时最容易垮的地方：

1. **划清了能自动化与不能自动化的边界**：六个模块 A–F 自动跑，
   三条诊断（共定位+多信号+匹配 LD、细胞层面 lineage 匹配、release 间 endpoint 逐码核对）
   **点名列出、每次报告都印**，而不是悄悄省略。
2. **输入契约已固定且极窄**：必需的只有两张表（`mr_results` 与 `known_loci`），
   其余全可选。没有 Seurat、没有 GWAS 对象、没有数据库。
3. **demo 锁在已发表数字上**，且**确定性**——我刚重跑一遍，
   输出与仓库里已提交的 `demo_out/` **逐字节相同**。

所以剩下的工作量比"从零做一个包"小一个量级。

---

## 二、结论与推荐

> **推荐：移植成纯 R 包，保留 Python 版本作为 CI 里的参照实现（oracle），
> 包名沿用 `cqtna` 不改。**

三条路的取舍：

| | 方案 | 评价 |
|---|---|---|
| (a) | 纯 R 重写，Python 版退休 | 会永远存在"两处实现悄悄分叉"的风险，但退休后无从校验 |
| (b) | R 包用 `reticulate` 包 Python | 只有一处实现，但**引入 Python 依赖 = 放弃采用率**，与做 R 包的初衷矛盾 |
| **(c)** | **纯 R 重写 + Python 留作 CI oracle** | **推荐。** 两处实现，但 CI 每次跑 demo 数据比对二者输出，分叉当场暴露 |

(c) 正是本项目一贯的纪律（step119/120/121 都是"重算一遍并自动核对"）。
让人不敢选 (a) 的那个风险，恰好被 (c) 消掉。

**包名不改**：正文已经写了 CQTNA，改名会让论文与工具对不上。

---

## 三、为什么值得做 R 包（采用率论证）

CQTNA 现在是"导出 TSV → 跑 Python CLI → 读回结果"。它的目标用户是 cis-MR 从业者，
而那个生态**整个在 R 里**：`TwoSampleMR`、`coloc`、`MendelianRandomization`、`gwasglue`。

没有人会为了跑一次审计，把 `TwoSampleMR::mr()` 的结果落盘、切到 Python、再读回来。
**摩擦不在算法，在换语言。** R 包能直接吃 `TwoSampleMR` 的 harmonised data frame，
这一步省掉，工具才有可能真被人用。

顺带的好处：论文的核心批评是"你们缺这些检查"，
而一个能直接接在 `TwoSampleMR` 后面的包，把批评变成了可执行的替代方案。
这是回应"那你说该怎么办"最有力的形式。

---

## 四、API 设计（R 原生入口）

```r
library(cqtna)

au <- cqtna_audit(
  mr       = dat,              # TwoSampleMR 输出，或含 chr/pos/gene/p 的 data.frame
  known    = known_loci,       # 你这个结局自己的已知 lead SNP
  mismatch = known_loci_hcc,   # 错配名单（缺失时**大声警告**，不静默跳过）
  build    = "GRCh38",         # ★ 见 §六
  locus_kb = 1000, known_kb = 1000, fdr = 0.05
)

au                          # print：分层报告 + 三条不能自动化的诊断
summary(au)                 # 计数表
plot(au, "attribution")     # Fig 1 式
plot(au, "distance")        # 连续距离分布（S36）
cqtna_report(au, "out.md")  # 与 Python 版同格式的 markdown 报告
```

各模块也单独可用（对应 Python 的 `module_a`–`module_f`）：

| R 函数 | Python | 做什么 |
|---|---|---|
| `cqtna_attribution()` | `module_a` | 按独立位点的归属 + 错配对照 |
| `cqtna_unit_sweep()` | `module_b` | 同一份结果在 record / variant / gene / locus 四种单位下重算 |
| `cqtna_stability()` | `module_c` | 换第二个结局 GWAS 后的名单稳定性 |
| `cqtna_ladder()` | `module_d` | 工具变量三级衰减 |
| `cqtna_compartment()` | `module_e` | 各细胞类型表达比 |
| `cqtna_peak_distance()` | `module_f` | eQTL 峰与 GWAS 峰的距离 |
| **`cqtna_window_sweep()`** | **（新增）** | **见 §五** |

---

## 五、R 版应当比 Python 版多出来的两件事

**1. `cqtna_window_sweep()` —— 把 step121 收进来。**
Python 版没有窗口敏感性。本窗口刚做出的 `step121` 证明这件事值得每个用户都跑一遍：
两个"1 Mb"约定分开扫，收紧窗口后倍数**上升**（4.09 → 18.26×），
并给出**连续距离分布**而非二分类。
对使用者的意义很直接——**它能告诉你"1 Mb"这个阈值是落在数据的空隙里，还是正踩在数据上**。
本文四格里三格落在空隙里，HCC-high 那格踩在数据上（唯一的已知命中在 408 kb 外），
这正是应该被报出来的差别。

**2. 构建版本必须显式声明并强制一致。**
Python 版的 README 写"any consistent build"，靠使用者自觉。
本项目自己就被这件事咬过（v6 坑 #35：eQTLGen 是 GRCh37 而 eQTL Catalogue 是 GRCh38）。
R 版应当**要求 `mr` 与 `known` 各带 build 属性，不一致直接报错**，
而不是算出一个静悄悄错到底的结果。这类错误的特征是数字看着完全正常。

---

## 六、测试套件 = 论文自己的数字

这是这个包不寻常的地方：**它的回归测试是一篇已发表论文的结果。**
`testthat` 直接断言下列各项（全部已在 `demo_out/` 里，且确定性）：

| 断言 | 值 |
|---|---|
| 归属：背景 | 58/554 = 10.47% |
| 归属：显著 | 3/7 = 42.9% |
| 归属：倍数与单侧 Fisher P | **4.09×，P = 0.0281** |
| 错配对照 | **0.00×，P = 1** |
| 单位扫描 | record 10 基因 / locus 28 基因（**同一份结果**） |
| 名单稳定性 | 10 vs 6，共享 5，Jaccard 0.455 |
| 工具变量阶梯 | **28 → 3 → 2 → 1** |
| 窗口扫描（新增） | 100 kb → 18.26×；1000 kb → 4.09× |
| eQTLGen 那格 | 20/30，4.44×，P = 3.7×10⁻¹¹ |

外加一条 **oracle 测试**：同一份 demo 输入分别跑 R 与 Python，
断言 `cqtna_results.json` 等价。分叉当场炸。

R 侧的等价函数都存在且是精确的，**移植风险低**：

| Python | R |
|---|---|
| `fisher_greater`（超几何上尾） | `phyper(a-1, ..., lower.tail = FALSE)` |
| `bh` | `p.adjust(method = "BH")` |
| `np.searchsorted` | `findInterval` |
| 1 Mb 单连锁聚类 | `cumsum(diff(pos) > w)` |
| `simes` | 三行 |

---

## 七、移植中**必须原样保留**的东西

这几条是 CQTNA 的设计立场，不是实现细节。换语言时最容易被"顺手优化"掉：

1. **没有任何路径能走到 `target-supported`。** 报告只发
   `screened` / `unresolved` / `state-informative` 三档。
   能授予那个词的检查，恰好是它跑不了的那三条。
2. **三条不能自动化的诊断每次都印**，且带各自的 failure mode
   （`cqtna_manual.md` 里写得很好，直接搬成 vignette）。
   ⚠ R 里尤其危险：用户会 `as.data.frame(au)` 然后把警告丢掉。
   所以对象里要带 `not_automated` 字段，`print` 与 `cqtna_report()` 都必须渲染它。
3. **错配名单缺失时大声警告**，不静默跳过——没有它就分不清
   "结局特异的归属"与"在所有疾病里都密集的位点"。
4. **已知位点名单的 provenance（`source` 列）与循环性提醒**：
   若名单里有来自结局 GWAS 自己的论文，那这个检验就循环到那个程度；
   且**剔除这些位点只会抬高富集**，过度修正不是更严格的检验。
5. **"A tool for detecting selective emphasis should be run first against the study
   proposing it."** README 最后那句留着。

---

## 八、工作量与风险

**工作量**：约 1–2 个专注窗口。
469 行 Python → 估计 600–700 行 R（含 roxygen）+ 测试。
不大，因为边界已经划好、契约已经固定、demo 已经是测试夹具。

**风险与对策**：

| 风险 | 对策 |
|---|---|
| 两处实现分叉 | CI oracle 测试（§六） |
| 构建版本静默错配 | 强制声明 + 报错（§五-2） |
| CRAN 要求：无网络测试、数据够小 | demo 夹具需裁剪；**建议先 GitHub + Zenodo DOI**（正文本来就写"deposited at ⟨repository DOI⟩"），CRAN 之后再说 |
| 依赖膨胀 | 硬性上限：**base + stats 足够**，画图用 `graphics`。不引 Seurat、不引 tidyverse |
| 模块 E/F 在 demo 里是空的 | 不是缺陷（demo 没配那两张表），但**发布前该配一份小样例**，否则用户以为坏了 |

---

## 九、需要你拍板的三件事

1. **要不要做**——若这一篇的目标是尽快投出去，R 包属于加分项而非阻塞项，
   可以放到投稿之后；但它作为"审稿人问'那该怎么办'"的答复很有力。
2. **范围**：只做 A–F + 窗口扫描（推荐），还是也把单细胞侧的
   (v) 区室归属 / (vii) lineage 匹配做成 R 函数（会把 Seurat 拖进依赖）。
   —— 推荐**不做**：CQTNA 现在用纯表格接口绕开了这个依赖，那是个好决定，别推翻。
3. **发布路径**：GitHub + Zenodo（推荐，与正文的 repository DOI 一致），还是直奔 CRAN。

---

## 十、若开工，第一步

不是写 R 代码，是**先给 Python 版补 `module_g`（窗口扫描）并让 demo 覆盖它**——
这样 R 移植过去时，oracle 里已经有这一项可比对。
顺序反过来会导致新功能只有一处实现、没有参照。
