# 稿件工作笔记（不进正文）

> 从 `MANUSCRIPT_v2_dual_thread.md` 抽出，2026-08-12（Step 90）。
> 抽出的理由：审稿 R3 指出正文文件里混着中文重构说明与待办，不是投稿稿的样子。
> **正文文件从此只含正文**；本文件是它的工作笔记。

---

## 一、v2 重构说明与备选标题

**Draft v2 — 2026-08-11 ｜ 双主线重构，取代 `MANUSCRIPT_full_v1.md`**

> **重构说明（中文，不进正文）**
> v1 是单主线（方法学审计 + TPI1 作为有界示例）。v2 拆成两条主线，接缝定义在
> **暴露侧／结局侧**：凡依赖结局 GWAS 的推断都不稳（主线一），凡只依赖暴露数据的观察
> 都稳（主线二）。这不是特批，是结构性的——单工具变量下 `z = β_out/se_out`，
> MR 的一切由结局决定；反过来，"工具变量在哪个时点出现"完全由暴露数据决定。
> TPI1 在两条主线各出现一次，证据等级不同，互不借力。
> 数字口径见文末；⟨…⟩ 为待补。
>
> **备选标题**
> 1. *Dynamic eQTLs are more reliable for locating when genetic regulation can be measured than for nominating what to target*
> 2. *Target nomination fails, timing survives: an audit of dynamic-eQTL Mendelian randomization in melanoma*


---

## 二、9. 重构后的待办（不进正文）

**结构层面**
1. **Fig 9 需新作**（Part II 的主图），面板见上
2. Part II 的三节目前引用了原 Results 的数据，需确认无一处经由黑色素瘤 MR
3. Abstract 已改为双主线；Title 待定（三选一或另拟）

**内容层面**
4. 三篇对照文献著录；16 条 ⟨verify⟩ 引用
5. Fig 5 图注与实际面板对齐
6. Methods §5b 已含 release 轨迹；多信号敏感性已写入 §5.5
7. **B1 已执行（§2.2、Methods §9b、S20）**：泛化到 HCC，登记为"方向一致但主检验未达显著"。
   待定：是否给 Fig 7（跨癌种）加一个 HCC 面板，或另作一图并入 —— 与 B2 的压缩要求冲突，
   建议**不加图**，留表在 S20
8. HCC 结局的引用（Gellert-Kristensen et al., *JHEP Rep* 2025, PMID 40823170）
   与 FinnGen R12、GWAS Catalog GCST90809296 三条著录需补入 `REFERENCES.md`

**口径（勿再改动）**
- 位点富集报**独立位点级 4.09×（P=0.028）**，记录级入补充
- HCC 泛化：**G1 主检验 P=0.110 未过**，只能写"方向一致的支持性证据"；
  HCC-high/low **不得**称两个独立队列；不得用 HCC 支持 TPI1 或糖酵解任何主张
- 功效曲线校准点 **10.2 vs 10**
- 单细胞糖酵解模块 **22 基因（不含 ENO3）**；GSE282266 的轴 **16 基因**
- 单细胞患者层面数字取 `32a`–`32h`
- release 轨迹：已知位点支有真实纵向证据；**新位点支只能说"不出现"，不得说"反复更替"**
- 多信号：**PARP1 排除，ZFYVE19/TPI1 无法判定**，不得写成"整体排除"
- **Part II 任何一句都不得出现黑色素瘤风险**

---

## 三、待核实后才可写入正文的一条（2026-08-13）

审稿人指出可用 **Chen KY, Kibayashi T, Giguelay A, et al. Genome-wide CRISPR screen
in human T cells reveals regulators of FOXP3. *Nature* 2025;642:191–200,
doi:10.1038/s41586-025-08795-5**（PMID 40140585，**论文身份已核实**）作为
TPI1 的真实实验约束：据称 TPI1 有 4 条 sgRNA，在该 FOXP3 诱导表型上
**不是显著命中**（gene-level FDR 约 0.999 与 0.480）。

⚠ **这两个 FDR 数字尚未核实**：该文的补充表在 Springer 静态路径上只取到一个 17 KB 的
附件，基因组级筛选结果不在其中（可能在 Source Data 或需订阅）。
**已从正文撤下，不得在核实之前写入。**

核实途径（按优先级）：① 期刊 Source Data 文件；② 作者 GitHub/Zenodo 沉积；
③ BioGRID ORCS 收录的该筛选。核实后应写成"**constrains but does not exclude**"，
且只能放 Discussion 的边界段，**不得**当作 TPI1 的功能验证。
