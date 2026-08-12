# 项目交接文档

**项目**：CD4+ T 细胞活化动态 eQTL → 黑色素瘤的孟德尔随机化研究
**目标期刊**：Research / JAR 一档（IF ~10 综合性）
**最后更新**：2026-08-08

> 新会话请先读本文件，再读 `FINDINGS_step5_pigmentation.md`（按步骤记录的完整分析日志）。

---

## 一、论文定位（重要，已多次调整后确定）

**主线是方法学，靶点是应用示例。** 不押单一靶点。

原因：候选名单在分析过程中被推翻过三次（见第四节），说明"找到 N 个免疫靶点"这类主张
在当前数据功效下不稳定。而四个方法学发现互相支撑、可迁移，且都指向三篇参考文献的共同盲区，
N=1 也成立。

参考文献（`D:/文章/孟德尔/pdf_text/`）：
- `s41467-024-53621-7.txt` — Nat Commun 2024，CD4+T 动态 eQTL → T2D/CAD（框架原型）
- `s41698-025-01236-6.txt` — npj Precis Oncol 2026，同框架 → 结直肠癌（**与本方案几乎同模子**）
- `ADVS-12-e07451.txt` — Adv Sci 2025，六癌种 + MR-DEG 方法
- `s41586-022-04713-1.txt` — Nature 2022 Nathan et al.，连续细胞状态 eQTL（非同一数据）

---

## 二、数据与工具位置

### 工作目录 `D:/R_ex/MR/`

| 类别 | 文件 |
|---|---|
| 暴露 eQTL（8 profile 全 cis 窗口） | `D:/Downloads/CD4_eqtl_step1_clean/*.parquet`（各约 840 万行）|
| 工具变量（每 exposure 1 SNP） | `D:/Downloads/CD4_eqtl_step2_instruments/`、`CD4_dynamic_top_eqtl_instruments_all.csv` |
| 结局：FinnGen R12 黑色素瘤 | `finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz`（5,753 例）|
| 结局：Rashkin 2020 | `rashkin2020_melanoma.h.tsv.gz`（6,777 例）|
| **结局：meta（主用）** | **`meta_melanoma_final.tsv.gz`**（12,530 例 / 789,099 对照，N=801,629，s=0.015631）|
| 已知位点参照系 | `landi2020_known_loci_grch38.csv`（Landi 2020 三个 GWAS 的 157 个 lead SNP，带表型标签）|
| 共病表型 GWAS | `pheno/*.gz`（7 个 FinnGen 表型）|
| bulk ICB 队列 | `icb_bulk/GSE78220_FPKM.xlsx`、`icb_bulk/GSE91061_fpkm.csv.gz` + series matrix |
| LD 参考面板 | `ref/1kg_eur_meta.{bed,bim,fam}`（1000G hg38，525 EUR 非亲缘 × 116,600 SNP）|
| 二进制工具 | `bin/plink2.exe`(v2.0.0-a.7.2)、`bin/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe` |

### 单细胞数据

- `D:/Downloads/GSE120575_*`（Sade-Feldman，**有应答标注**，在 `_patient_ID_single_cells.txt.gz` 第20行起）
- `D:/数据/黑色素瘤单细胞人/GSE115978`（Jerby-Arnon，自带 cell.annotations，**无应答标注**）
- 同目录还有 GSE72056、**GSE123139**（黑色素瘤增殖性 T 细胞专题，未用，与假说契合）、
  GSE139829（葡萄膜黑色素瘤，不建议用）等
- `D:/Downloads/ST-Melanoma-Datasets_1/`（Thrane 空转，4 例 ×2 重复，**已用于 Step 19**）

### 关键分析脚本（均在 `D:/R_ex/MR/`）

`meta_finngen_rashkin.py`、`fix_meta_pvalues.py`、`rerun_mr_meta.py`、
`build_besd_meta.py`、`build_relaxed_meta.py`、`comorbidity_mr.py`、
`step6_coloc.R`、`figures/make_figures.py`

---

## 三、四个方法学发现（论文主线）

1. **位点归因**：FDR<0.05 的信号 100% 落在已知色素/痣位点 1 Mb 内（MC1R rs1805007 R151C 红发变异 48 kb、PARP1 5–14 kb）。三篇参考文献均未做此检查。→ **Fig A**

2. **HEIDI 失效**：coloc 判为不同因果变异（PP.H4<0.2）的 291 条记录中，**253 条（86.9%）被 HEIDI 放行**。原因是 nsnp_HEIDI 仅 8–20、GWAS 信号弱。→ 证据分级中 Level 3（MR+SMR 无 coloc）不安全，应以 coloc 为主判据。→ **Fig B**

3. **MR 与 coloc 反向变化**：GWAS meta 后 MR 增强（FDR<0.05 从 10→21）但 coloc 减弱（同一批 exposure 的 PP.H3+H4 中位 0.249→0.207）。因 MR 只依赖单个 top SNP 的 z 值，coloc 依赖整个区域信号形状。

4. **★ coloc 分辨能力由 GWAS 功效决定，误差双向**（最重要）：
   - PARP1：GWAS 峰 p 1.5e-6→1.2e-11、峰位移 25 kb，**PP.H4 0.92→0.05**（假阳性被消除）
   - ZFYVE19：meta 峰精确落在 eQTL 峰（0 kb），**PP.H4→0.99**（假阴性被揭示）
   - **两轮候选名单交集为零** → **Fig C**

5. **归属问题：MR 之后如何用功能数据定位效应**（两个配对的子发现）
   - **共位点基因归属**：chr12p13 上 SPSB2 与 TPI1 相距 3–20 kb，MR 无法区分。
     三条来源不同的正交证据均指向 TPI1：单细胞应答分层（SPSB2 p=0.27/0.13）、
     GSE115978 检出率（5% vs 60.6%）、**eQTL 时间模式（SPSB2 组成型 vs TPI1 16h 脉冲）**
   - **区室归属**：候选基因若广泛表达，组织水平复现**在原理上无法支持细胞类型特异的因果主张**。
     空转裁决 TPI1 的组织信号由肿瘤糖酵解驱动（vs 糖酵解模块 ρ=+0.174 8/8 切片、
     vs 淋巴 ρ=−0.080），**裁决结果为否定** —— 阳性对照 HLA-C 呈完全镜像（淋巴 +0.187 8/8）
     证明该检验有分辨力。→ **Fig E**

6. **工具变量可得性、而非生物学重要性，决定 MR 点名哪个基因**（Step 28–29）
   糖酵解通路：功能上 22 个酶中 20 个在非应答者升高、8 个过 FDR（Post），
   **遗传上 28 个基因中只有 TPI1 与 ENO1 有工具变量，比例 1/22**。
   基线最强的 PGAM1（FDR=0.032）根本没有工具变量。
   通路级富集检验亦为阴性（mean |z| 0.890 vs 匹配背景 0.937，p=0.52，但 n=9 功效低）。
   → **通路级主张不成立，必须停在单位点表述**；三篇参考文献同样未声明此边界。

7. **在细胞群内部按评分切分会同时切出细胞类型纯度梯度**（Step 24）
   深度匹配**不能**修复；必须额外按谱系评分匹配。
   且 **cluster 层面对照会通过而细胞层面对照失败** —— 对照必须做在细胞层面。
   实例：GSE120575 中 Glyco-hi/lo 的 CD4 检出率 40.1% vs 19.5%、CD8A 47.6% vs 63.9%；
   按 CD8ness+深度双重配对后（610 细胞、28 样本，残余 SMD 0.0097/0.0000）
   全部模块归零，去谱系基因的 EffectorClean 方向 14/14 恰好随机。

---

## 四、候选基因现状

**四重验证（MR+coloc+SMR/HEIDI+敏感性）存活 11 条 / 新位点 6 基因**：
ZFYVE19、SMC2、KIAA0040、SPSB2、HLA-C、TPI1（5 个独立位点，SPSB2/TPI1 同位点）

### 最终排序（**TPI1 领先**）

**scRNA 应答一列已按 Step 26 修正版更新（`32c`），旧的 `19h` 数字作废。**

| 基因 | MR | 不影响痣 | 零多效 | scRNA 应答（FDR）| 组织水平 | 成药性 |
|---|---|---|---|---|---|---|
| **TPI1** | OR=1.305 | ✓ | ✓ | **✓ Pre 0.046 / Post 0.017**（唯一两时点均过）| ~ 方向一致但为肿瘤源 | ✗ 必需基因 |
| **SMC2** | OR=1.415 | ✓ | ~ | **✓ Post 0.017**（Pre ns） | ✗ | ✗ 必需基因 |
| **HLA-C** | OR=1.078 | ✓ | ✗ | **✓ Post 0.048**（Pre ns）★新增 | （空转阳性对照）| ✗ SJS/TEN |
| ZFYVE19 | OR=0.934 | ✗ 影响痣 | ✓ | ✗ 两时点均 ns | ✗ | ✓ 最优 |
| SPSB2 | OR=0.957 | ✓ | ✓ | ✗ 两时点均 ns | ✗ | ~ |
| KIAA0040 | — | ✗ | ✗ | ✗ Post 0.127（较原报告更弱）| 方向相反 | — |

另：糖酵解模块本身 Post FDR=**0.0152**，是该家族中最强信号，强于任何单基因。

### ⚠ SMC2 定位已下调（Step 30 跨癌种，2026-08-10）

**SMC2 是泛癌的，不能再作为"CD4 免疫监视"或"细胞分裂模块支撑免疫假说"的佐证。**

| | 黑色素瘤 | 肺 | 乳腺 | 前列腺 |
|---|---|---|---|---|
| SMC2 | OR=1.415 (1.5e-5) | **1.384 (5.7e-4)** | **1.304 (2.5e-5)** | **1.232 (0.0044)** |

四个癌种同向且显著 → 走的是**普适的增殖/基因组稳定性**机制。
与既有线索一致（DepMap 必需基因、pLI=0.99999、与增殖评分的相关近乎循环论证），
现在有了外部遗传学证据。

**连带影响**：FINDINGS Step 13 提出的「主打 ZFYVE19 + SMC2 细胞分裂模块、
两个独立位点收敛到同一生物学、支撑 CD4 增殖影响黑色素瘤易感性」这一定位**作废**——
ZFYVE19 在胰腺/乳腺/前列腺同样显示保护效应，两者都偏泛癌，
收敛的是**泛癌增殖生物学**而非黑色素瘤免疫监视。
SMC2/ZFYVE19 可作为「同一流程也会捞出泛癌增殖基因」的对照来写，不作免疫论据。

患者层面相关性（`32g`，修正版）：SMC2~增殖 ρ=0.693/0.812（FDR 0.008/9.8e-6）、
TPI1~增殖 0.549/0.746、TPI1~耗竭 0.726/0.617（FDR 0.004/0.004）。
**HLA-C~增殖的强负相关（旧版 −0.612/−0.720）已证实为二次归一化假象，现为 +0.14/−0.12 (p≈0.55)。**

**TPI1 证据链方向一致**：表达↑ → 黑色素瘤风险↑（MR）→ ICB 应答差。

⚠ **但支撑"CD4 特异"的证据只有两条**（Step 19 空转后收窄）：
1. MR 工具变量按构造来自 CD4+T 动态 eQTL —— 这是因果主张的真正依据
2. GSE120575 CD4+T 分辨的应答分层（治疗前 FDR=0.034 / 治疗后 0.015）—— **仅此一个队列**

bulk 复现（GSE78220 p=0.045、GSE91061-On p=0.030）**处理方式待定**：
空转（Step 19）显示组织水平 TPI1 由肿瘤糖酵解主导，据此应降级；
但降不降、怎么表述，等生物学线（Step 22/23）跑完一并决定。数据本身见 FINDINGS Step 19。

### 曾被推翻的候选（写进 Discussion 作为筛选有效性的例证）
- PADI4 —— 敏感性分析中效应塌缩（OR 1.087→1.004）
- IMPA1 / PRPSAP2 / GCC2 —— 换 meta 结局后全部掉出
- HLA-C 的"增殖机制"解释 —— GSE115978 中方向反转，**已撤回**

---

## 五、必须写入论文的局限

1. **GSE91061 治疗前未复现 TPI1**（p=0.814），而治疗前是临床最有价值的比较
2. **bulk 复现确实反映肿瘤糖酵解（Warburg）** —— 空转已证实（Step 19），不再是"可能"。
   TPI1 的 CD4 特异证据仅剩 MR 工具变量构造 + GSE120575 单队列。
   空转本身另有局限：100 µm spot 含 10–40 个混合细胞，CD4+T 为少数群体，
   故只能下"组织信号由肿瘤糖酵解主导"这一正面结论，**不能**推出"TPI1 在 CD4+T 中无作用"
3. TPI1 与 SMC2 均为 DepMap 必需基因，成药性差
4. 当前候选名单同样受功效限制，更大 GWAS 下可能再次更替
5. 单 SNP Wald ratio 下 p 值由结局 GWAS 唯一决定 → 细胞类型/时间特异性不能用 p 值差异主张
6. Steiger 全部通过属数据结构必然（R² 差 4 个数量级 + winner's curse），不应作为因果方向强证据
7. 工具变量 r²<0.1 相关 → IVW 显著性被高估（IVW FDR<0.05 有 157 个 vs 单 SNP 严格集 21 个）
8. Rashkin 无 SE，由 `se=|log(OR)|/|Φ⁻¹(p/2)|` 反推（Cochran Q 显著率 5.33%，校准良好）

---

## 六、技术坑（重新踩会很贵）

1. **坐标是 GRCh38，不要做 liftover** —— 曾误判为 hg19 并 liftover，导致匹配率从 90% 掉到 2%
2. **不要用 OpenGWAS 的 `finn-b-C3_MELANOMA_SKIN_EXALLC`** —— R5 快照仅 98 例
3. **Landi 2020 无完整 sumstats**，只有 76 个位点，不能用于 MR
4. **公开数据无可用白癜风 GWAS**（11 项研究全无 sumstats；FinnGen 仅 391 例）→
   已改用「黑色素细胞痣」(13,357 例) 作色素通路阴性对照，阳性对照 5/5 成立
5. R 中 `sum(x & NULL)` 静默返回 0，会伪装成"无显著结果"；`harmonise_data` 已带出额外列，
   不要再 merge（会产生 `.x`/`.y`）
6. `AnnotationDbi::select` 会遮蔽 `dplyr::select`，arrow 管道需写 `dplyr::select`
7. TwoSampleMR 的 `units="SD"` R² 公式 `2β²·EAF·(1−EAF)` **无上界**，会算出 R²>1；
   已加 `rsq.exposure.bounded` 列用 `F/(F+N−2)`
8. meta 后 1,740 个变异 p 下溢为 0，已加 `mlogp` 列并对 p 设下限 1e-300
9. Cochran Q 两研究 meta 的 df=1，上尾是 `erfc(√(Q/2))` 而非 `exp(-Q/2)`
10. Open Targets 旧 Genetics API 已停用；新 Platform API 字段是 `drugAndClinicalCandidates`
    （无 `size` 参数），且按 rsID 查询对部分变异会映射到错误染色体
11. matplotlib 默认字体无中文字形 → 图全部用英文标注
12. **`statsmodels 0.14.4` 与本机 scipy 不兼容**（`scipy._lib._util._lazywhere` 已移除，
    `import statsmodels.api` 直接报错）。Step 19 脚本已用 numpy 自实现 OLS 与 BH-FDR，
    后续沿用即可，**不要花时间修环境**
13. Thrane 空转数据**无独立坐标文件**，坐标编码在列名（`2x9` = x2,y9）；
    行名为 `SYMBOL ENSG...`，需取第一段
14. **`GSE120575_CD4_clusters_1_5_6_12_15.rds` 被二次归一化**（把 log2(TPM+1) 当 counts
    喂进了 `NormalizeData()`）。与正确对象的非零值 ρ 仅 0.874。
    正确对象是 `24_checkpoint_GSE120575_annotated_obj.rds`（反解每细胞和 9.94e5 ≈ 1e6）。
    该 rds 仅可用于取细胞 ID，**不可用于取表达值**。
    修正后的定论数字见 `32a`–`32h`（Step 26），旧的 `19*`/`28*` 患者层面数字作废。
    ⚠ **GSE115978 没有这个问题**：其原始文件每细胞总和 CV=1.15、跨 46 倍，本就未归一化，
    故 `GSE115978_TCD4.rds` 的 LogNormalize 是正确处理，Step 18 结论有效。
    **判据是原始文件是否已归一化，不是 rds 与原始文件是否一致**
15b. **10x Multiome 的 RNA 是 snRNA-seq（细胞核）**。GSE282266 中 MALAT1 以 142 万
    counts 居首（96.5% 核）、LINC00486/TALAM1 等核内 lncRNA 高，而 GAPDH 只在 18.7%、
    ACTB 52.1% 的核中检出——**这不是错误，是核提取的正常特征**。
    含义：胞质 mRNA（含糖酵解酶）在核里天然偏低，**per-cell 检出率不可与 scRNA 比较，
    RNA 必须按样本 pseudobulk 分析**；ATAC 侧不受影响且很深（peak 2.11 亿 vs 基因 3,380 万）。
15c. **GSE282266 每个样本的 peak 集不同**（rest_set_1 有 97,119 个，act_15_set_1 有 132,780 个），
    且 `pooled_peaks.bed.gz` 只有 54,506 个并含重复行。
    跨样本比较必须**按基因组窗口定义区间、用各样本自己的 peak 求和**，不能按行号对齐。
    features 文件里基因条目自带坐标，可直接取 TSS。
15d. **工具变量都不落在 peak 内**：TPI1 rs12302749 距最近 peak 92 bp，
    SPSB2 的三个距 1.1–4.1 kb；且离 TPI1 最近的那个 peak 同时是 SPSB2 某工具变量的近邻。
    **原设想的"TPI1 元件 vs SPSB2 元件"内部对照在此 peak 分辨率下做不干净**，
    已改为基于 TSS 窗口（±2 kb 主、±10 kb 敏感性）的基因水平可及性。
15e. **`19i_patient_tp_with_response.tsv` 的应答标注不可用**：GEO 中 Post_P1 是 Responder、
    Post_P1_2 是 Non-responder（P5 亦然，方向相反），19i 把它们并成一个样本一个标签。
    应答标注必须取自对象的 GEO `response` 字段，样本 ID 用 GEO 的 `patient_tp`。
    同一患者 Pre 与 Post 的标签也可能不同（序贯治疗），Methods 需说明

---

## 七、进度与下一步

### 已完成
数据准备 / harmonise / MR / 位点归因 / coloc / SMR-HEIDI（官方 v1.3.1）/ Steiger /
GWAS meta / 敏感性分析 / 共病与通路归属 / 多效性与成药性 / 单细胞（GSE120575+GSE115978）/
bulk ICB 复现（GSE78220+GSE91061）/ **空转区室归属（Step 19，结果为阴性，见下）** / 五张主图

### 待做（按优先级）
0. **TPI1 生物学线** —— Step 22（eQTL 16h 时间窗）、Step 23（糖酵解通路）、
   Step 20（外部完成的 Post-ICB 表型）、Step 24（纯度对照）**均已完成**，见 FINDINGS。
   ⚠ **Step 20 主结论已被纯度对照推翻**，但 TPI1 应答分层稳健（Step 24c）
1. **跨癌种 MR** —— 回答特异性问题，工具变量/BESD/LD 面板全现成，只需换结局 GWAS。
   线索：SMC2 的 rs2417487 在 GWAS Catalog 中关联胰腺癌 (p=6e-7)
2. **Fig F**（TPI1 三队列应答分层）、**Fig G**（痣 vs 黑色素瘤通路归属对比）—— 数据已在手
   （原计划的 Fig E/F 顺延，Fig E 已被空转图占用）
3. TCGA-SKCM 生存/免疫浸润（可选，标准动作）

### ⚠ 已撤回的结论（不要重新捡起）
**GSE120575 Post-ICB「糖酵解高 CD4 呈 helper/调节倾向、细胞毒性弱」** —— 该表型是
CD4/CD8 谱系构成假象。按 CD8ness + 深度双重配对后（610 细胞、28 样本，残余 SMD
0.0097/0.0000），CytotoxicScore p=0.206、HelperRegScore p=0.104、去谱系基因的
EffectorClean p=0.882（方向 14/14 恰好随机）；FOXP3 p=0.727。
**教训：在细胞群内部按评分切分会同时切出纯度梯度，深度匹配不能修复，必须按谱系匹配。
且 cluster 层面对照会通过而细胞层面对照失败——对照必须做在细胞层面。**
详见 FINDINGS Step 24。CellChat 的 CCL/CCR5 结果随之失效（所比较的两组构成本就不同）。

### 已做但结果为阴性（不要重做）
**空转分析**（Thrane ST，8 切片 2,317 spots）—— 目的是补 TPI1 的 CD4 特异性漏洞，
**结论是漏洞确实存在**：TPI1 跟肿瘤/糖酵解区室（8/8 切片），与淋巴区室负相关。
阳性对照 HLA-C 呈镜像（淋巴 +0.187，8/8）证明检验有分辨力。
净效果：方法学主线 +1 个发现（区室归属），靶点侧 bulk 证据降级。详见 FINDINGS Step 19。

### 明确不做
**药筛（CMap + 分子对接）** —— TPI1 是 DepMap 必需基因、pLI=0.87，自己的成药性分析已判定
毒性风险高，再找抑制剂等于自相矛盾；且与方法学主线无关，属稀释重点。

**46-profile 扩展** —— 有数据证明无效：瓶颈在结局 GWAS 分辨率而非暴露侧覆盖
（coloc PP.H1 中位数仍达 0.83，多数区域根本无 GWAS 信号）。且本地只有 8 个 profile，
其余 38 个需重新获取清洗。

---

## 八、主要输出文件索引

| 文件 | 内容 |
|---|---|
| `12_MR_meta_strict.tsv` | meta 结局下严格集 MR（3,556 检验，FDR<0.05 有 21） |
| `13_meta_locus_annotation.tsv` | 上表 + category/near_locus/SYMBOL |
| `15_SMR_meta_results.tsv` | SMR/HEIDI（369 条） |
| `16_coloc_meta_results.tsv` | coloc（371 条，均质化版） |
| `17_relaxed_instruments_meta.tsv` | 放宽工具变量（1,601 个，354 exposure ≥2 SNP） |
| `18_sensitivity_meta.tsv` | IVW / IVW-MRE / weighted median |
| `19h_responder_vs_nonresponder.tsv` | GSE120575 应答分层 |
| `20_comorbidity_MR.tsv` | 13 基因 × 7 表型共病 MR |
| `22b_GSE115978_correlation.tsv` | GSE115978 复现 |
| `23a`–`23h_ST_*.tsv` | 空转区室归属（`23g` 功效不足，勿引用）|
| `26c_ST_compartment_colocalization.tsv` | 区室间空间共定位（已跑，见第九节）|
| `figures/Fig[A-E]*.{png,pdf}` | 五张主图 |

---

## 九、进行中：TPI1⁺CD4 的功能故事（CellChat → 空转）

**思路**：不再用空转证 CD4 特异性（Step 19 已证伪），改为回答「TPI1 高的 CD4 在做什么」。
先在单细胞里做 CellChat 找 TPI1-hi CD4 的差异通讯，再用空转验证配体-受体是否空间共定位。

### 脚本（已写好，未跑）

| 脚本 | 用途 | 环境 |
|---|---|---|
| `step20_install.R` | 装 CellChat + NMF（跑一次）| R 4.4.1 |
| `step20_cellchat_GSE120575.R` | 免疫-免疫网络，**有应答标注** | R 4.4.1 |
| `step20b_cellchat_GSE115978.R` | CD4↔肿瘤轴，**无应答标注** | R 4.4.1 |
| `step21_spatial_LR_coloc.py` | 空转验证 LR 共定位（已冒烟测试通过）| python |

**分组方式（Step 23 后已改）**：主分组是**糖酵解模块**（检出率≥30% 的基因、剔除 TPI1、
`AddModuleScore` 表达分箱匹配），TPI1 单基因分组保留作**敏感性分析**。
理由：功能表型是通路级的，模块比单基因稳、功效更好、且分箱匹配本身削弱深度混杂。
脚本一次跑完两种分组，并输出 `24g_split_concordance_pathways.tsv` 比对通路重叠——
**两种分组结论不一致则两个都不可信**。

空转脚本接命令行参数：
```
python step21_spatial_LR_coloc.py 24f_Glyco_top_LR_for_spatial.tsv
python step21_spatial_LR_coloc.py 24f_TPI1_top_LR_for_spatial.tsv
```

⚠ **R 环境**：干活的是 `D:/R/R-4.4.1`（Seurat 5.5.1 在这里），
**不是** R 4.5.0（那个库是空的）。

### 三个必须守住的设计点

1. **GSE120575 是 CD45⁺ 分选** —— 无肿瘤/成纤维/内皮细胞，做不了 CD4↔肿瘤。
   该轴只能用 GSE115978（但它无应答标注）。**两个队列的结论不能混着讲。**
2. **按 TPI1 高低劈 CD4 会引入测序深度混杂** —— TPI1 是高表达管家基因，
   TPI1-hi 细胞多半只是文库更深，而 CellChat 通讯强度随检出率上升，
   不校正则"TPI1-hi 通讯更活跃"是必然的假结果。
   脚本内做**样本内按 nFeature 的 1:1 贪心匹配**，并强制打印匹配前后的 Wilcoxon p，
   **匹配后不显著才能往下解读**。
3. **空转不要去找"TPI1 高的 CD4 spot"** —— 100 µm spot 做不到，且 Step 19 已证明
   spot 层面的 TPI1 就是肿瘤读数。可验证的命题是**配体与受体是否占据相邻组织**，
   不是 TPI1 本身。

### 空转验证用了两个零假设（关键）

- 空间置换：只检验"是否共定位"。几乎所有丰度基因都会通过（组织结构使然），**不足以支撑主张**
- **表达量匹配的随机基因对零假设**：检验"是否比同等丰度的任意基因对更共定位"。
  **这个才是能写进论文的判据**，脚本按它出 FDR

另有两条硬约束，脚本会自动分开报告：
- CellChatDB 标为 `Cell-Cell Contact` 的互作在 200 µm spot 间距下**原理上无法验证**，单列不得主张
- 大量免疫配体/受体（趋化因子受体、IFNG 等）低于 5% 检出下限 → 输出 `26d_ST_untestable_pairs.tsv`。
  **"测不了"不等于"被否定"**，两者不可合并

### 已跑出的真实结果：区室空间共定位（`26c`）

| 区室对 | 双变量 Moran's I | 8 张切片中显著 |
|---|---|---|
| 淋巴 – B 细胞 | **+0.245** | 6/8 |
| 淋巴 – 髓系 | +0.125 | 7/8 |
| 髓系 – B 细胞 | +0.043 | 4/8 |
| 髓系 – 肿瘤 | −0.030 | 1/8 |
| **淋巴 – 肿瘤** | **−0.205** | **0/8** |
| B 细胞 – 肿瘤 | −0.224 | 1/8 |

**T 细胞与肿瘤在空间上互斥（0/8 切片），且淋巴与 B 细胞聚集共定位（TLS 样）。**
这是讲故事的空间底色：若 TPI1⁺CD4 有作用，其空间语境是淋巴聚集灶而非肿瘤接触面。
也意味着 CD4↔肿瘤的直接接触型互作在这批标本里**本就不该期待被空转验证**。
