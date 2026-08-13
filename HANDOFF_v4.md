# 项目交接文档 v4

**项目**：CD4⁺T 细胞活化动态 eQTL → 靶点提名的审计（黑色素瘤 + HCC，两套暴露资源）
**最后更新**：2026-08-13（**Step 99 后**）
**取代**：`HANDOFF_v3.md`（其第六节技术坑、第七节方法论纪律仍有效，本文件继承并扩充）

> **新会话阅读顺序**
> 1. 本文件
> 2. 两版投稿稿：`manuscript/MANUSCRIPT_GB.md`、`manuscript/MANUSCRIPT_NC.md`
> 3. 全文源 `manuscript/MANUSCRIPT_v2_dual_thread.md`（**唯一正文源**，两版都是派生物；
>    Step 99 的回填已把它恢复到这个地位，见 §三）
> 4. `FINDINGS_step5_pigmentation.md` 的 Step 84 起
> 5. **七份预注册**（S9/S18/S20/S21/S22/S24/**S25**）+ `SUPP_attempt_timeline.md`

---

## ★ Step 99（本轮新增，最重要的一件事）

**第七份预注册 S25（`manuscript/PREREG_r13_transfer.md`）已写、已跑、已填结果。**
判读落在 **A（转移成功）**：六项登记预测全部符合、三项过程对照全过、错配阴性对照为空。

| | R12 | R13 |
|---|---|---|
| 病例 | 5,753 | 6,226 |
| FDR<0.05 记录 / 基因 / 位点 | 10 / 6 / 2 | **10 / 6 / 2（逐字相同）** |
| 新位点基因 | 0 | **0** |
| 归属倍数（Soskic×melanoma）| 9.58×，P=0.011 | **9.58×，P=0.011** |

**为什么这不是"数据没变"**：3,434 条 z 值**无一相同**（corr|z|=0.937），
se 中位收缩 **3.85%**，与病例数 +8.2% 预期的 3.8% 吻合；
名义 p<0.05 的 275 条里**换掉了 51 条**（Jaccard 0.687）。
→ **阈上纹丝不动，阈下换掉近五分之一。**

**四条必须一起读的限定**（S25 §10.4–10.5，正文已按此写）：
1. P1 点预测 11 **没中**，观测 10 = 与 R12 完全相同；落在区间 [9,14] 但在"零增长"那一端。
   **"落在区间内" ≠ "预测准确"。**
2. C1 与 C3 的倍数/P 逐位相同是**算术**（同一张 2×2 表必然同一个 Fisher P），
   **不是第二次独立确认**。看起来像 bug，不是。
3. 275 = 275 是**巧合**，不是同一批记录。
4. **我自己的一个 bug，已修正并登记**：`92c_locus_annotated.tsv` 的 `beta_out/se_out`
   是对着 **meta**（FinnGen+Rashkin）算的，**不是 R12**。初版脚本拿它当 R12 基线，
   于是 C2 实际比的是 meta→R13，**同时换了 release 和数据集**——正是本文通篇在讲的混杂。
   症状是"功效上升而显著位点从 30 掉到 16"。修正后 C2 = 3.70×→3.74×（基本不动）。

**第八条诊断已立**，两版投稿稿都已加：
> **跨 release 比较候选名单之前，先逐条核对 endpoint 定义是否同一。**
> R13 manifest 把该 endpoint 写作 "including Hilmo"，照字面读会以为新纳入了住院病例；
> 逐条编码定义显示**病例定义相同、变的是 controls 排除规则**。
> 只读摘要行就比，会把纯粹的表型漂移记成候选名单的不稳定。

GB/NC 两版的 "seven checks" 已全部改为 "eight"（正文 5 处；
另有 3 处 "seven" 指的是**十一次尝试中七次推翻**，**不得改**）。

**ORCS 已从一节压成一段**，两版都是。措辞锁死在
"a knockout cannot isolate a CD4-specific role"，未改成"CRISPR 支持结论"。

---

## 一、当前定位

**一条方法学主线 + 一个被严格限定的 worked example。**（不再是"双主线"——
审稿人 R2 指出第二条线证据不足以对等，已接受。）

主线：**靶点提名是结局 GWAS 的性质，不是暴露侧的性质。** 支撑它的是一张
**疾病 × eQTL 资源的交叉网格**，六格全部同向，配一个干净的错配位点阴性对照；
Step 99 后再加**一次跨 release 的转移测试**（S25，判读 A），
把"换成读者今天能拿到的那份结局数据会怎样"这个最容易被当场证伪的问题也堵上了。

worked example：TPI1，定性为 **instrument-visibility 的示例**，不承担任何
动态遗传调控或靶点主张。

---

## 二、这一轮（Step 84–99）新增了什么

| Step | 内容 | 结果 |
|---|---|---|
| 84–86 | **HCC 泛化**（预注册 S20）| 4 个显著位点中 3 个是已知 HCC 位点（含 PNPLA3）；**主检验 P=0.110 未过** |
| 87 | **HCC 患者分层**（S21）| post 支 Δ=+0.874、p=0.043（与发现队列 +0.87 同量级）；pre 支事先不可能显著 |
| 89 | B2 结构压缩第一阶段 | Part I 10 节→7 节、Discussion 8→5、新增 Box 1 |
| 90 | 应第二轮审稿的"距离 8 分"七条 | 三层口径拆分、11 项四分类、统一 primary test、10 处矛盾清零 |
| 91 | **系统文献审计**（S23）| 152 篇全文：位点归属至多 **7.9%**，名单-结局依赖 **0.7%** |
| 92 | **换暴露资源**（S22）| eQTLGen 全血 n=31,684：**4.44× vs 4.09×**，P=3.7e-11 |
| 94 | **疾病×资源网格**（S24）| 六格全部 fold>1；错配阴性对照 1.30×/P=0.41 与 0.00×/P=1.0 |
| 95 | **更正一处科学事实错误** | genotype×pseudotime 数据其实公开；三个 TPI1 lead variant 均无交互 |
| 96 | **自我施加的 ground truth** | 公认因果基因位点上 **6/10 归对**；MC1R 点名 16 个基因却漏掉 MC1R |
| 97 | **ORCS 全库 CRISPR** | TPI1 命中率 42.7%，与 GAPDH/PGAM1 同带；**PGAM1 47.2% ≈ TPI1** |
| 98 | **FinnGen R13** | 数据已下载；endpoint 映射见第四节 |
| **99** | **R13 转移测试（预注册 S25）** | **判读 A**：名单逐字不变、四格倍数全部同向，而 3,434 条 z 无一相同；**第八条诊断由此确立** |

---

## 三、两版投稿稿的状态

| | GB | NC |
|---|---|---|
| 文件 | `manuscript/MANUSCRIPT_GB.md` | `manuscript/MANUSCRIPT_NC.md` |
| 摘要 | 382（结构化，目标 ~350）| **141**（限 150 ✅）|
| 正文 | 4,987（GB 惯例 6,000–8,000，**偏短，有余量**）| 2,830（限 5,000 ✅）|
| 结构 | 保留完整论证链，七条诊断为 Discussion 骨架 | 单一叙事线，无 Part I/II 之分 |
| Box 1 | 保留（紧凑分类版）| 无 |

**两版都已按 R3 的五类移出过程性内容**：失败全过程与修改历史、stopping-rule 实例、
处理错误的技术经过、候选选择时间线、重复的 claim boundaries —— 全部只在 Supplementary。

### ✅ "唯一正文源"已恢复（选了方案 a，回填完成）

Step 99 一度出现"三份文件各自为源"：ORCS（Step 97）与第八条诊断 + R13（Step 99）
只在 GB/NC，全文源没有。**已按方案 (a) 回填，现在全文源重新是唯一源，
GB/NC 都是纯派生物，没有任何自己独有的内容。**

回填的四处：

| 位置 | 内容 |
|---|---|
| §2.7「What cannot be claimed」之后 | **ORCS 段**（接在 DepMap-essential 那句后面，是它的延伸）。措辞锁死 "a knockout cannot isolate a CD4-specific role"，并明写 "not evidence for the nomination" |
| §2.2 五个 release 段之后 | **「Why the series stops at R12」+ R13 转移测试**两段，含 S25 §1 的方向不对称论证 |
| §4.5 | **第八条诊断 ⑧**；标题改为 "Eight diagnostics"；摘要加一句 R13 转移 |
| Supplementary 表 + Methods | **S25 行**（含四条自我登记项）；`METHODS_draft.md` 新增 **§5c**，`assemble.py` 会自动renumber 成 §5.5c |

⚠ **Step 96 的自我 ground truth（6/10）仍然三份都没写**，是独立的待办，不在本次回填范围。

### 字数现状

| | 词数 | 摘要 | 状态 |
|---|---|---|---|
| 全文源 | 15,880（回填前 14,664）| — | — |
| `MANUSCRIPT_assembled.md` | 23,099 | — | 由 `assemble.py` 生成，**永不手改** |
| GB | 6,005 | **370** | 正文进入 6,000–8,000 惯例区间 ✅；摘要仍超 ~350 目标 **20 词** |
| NC | 3,235 | **141** | 限 150 ✅；正文限 5,000 ✅ |

GB 摘要**回填前是 382**，加了 R13 那句之后压到 **370**——比原来还短，
但仍超目标。再砍就要动"In this analysis, and within the power range we could
observe"这类限定语，**那与本文的立场相悖，故停在 370，留给投稿时定夺。**

`assemble.py` 负责把 Methods 与 References 装配进全文版（生成 `MANUSCRIPT_assembled.md`，
**永不手改**）。`build_gb.py` 是早期的减法脚本，**已被手写 GB 版取代，可删**。

---

## 四、FinnGen R13：endpoint 映射已定，**分析已做**（结果见本文件开头的 Step 99 节）

### 映射结论（用户查 Risteys 逐项核对，已采纳）

| | R12（本文所用）| R13 |
|---|---|---|
| endpoint | `C3_MELANOMA_SKIN_EXALLC` | **`C3_MELANOMA_SKIN_WIDE`** |
| 病例 / 对照 | 5,753 / 378,749 | 6,226 / 372,159 |
| 疾病编码 | C43 · 172 · C44 + melanoma morphology | **相同** |
| controls 排除 | C3_CANCER + C3_CANCER_WIDE_EXALLC | C3_CANCER + **C3_CANCER_WIDE** |

→ **病例定义相同，差别在 controls 排除规则。**
⚠ 我在 Step 98 曾写"病例认定口径变了（纳入住院登记）"，那是从 manifest 一行描述
**推断**的，**已更正**：Risteys 的逐项核对显示 R12 的 `_EXALLC` 本来就含住院编码。

**不得改用裸 `C3_MELANOMA_SKIN`**（R13 把它列为更窄的相似 endpoint；
且已发布的 manifest 里根本没有这个 endpoint，只有 `_WIDE` 与 `C3_MELANOMA_WIDE`）。

本地已核实：`_EXALLC` 后缀在整个 R13 manifest 出现 **0 次**。

### 但它仍不是纯功效点，预注册必须事先处理

**对照池从 378,749 缩到 372,159，而队列总量是增长的** —— 排除规则确实变宽了。
所以 R13 = 功效增长（+473 例）**叠加**一个残余的对照定义差异。

**S25 §1 已给出解法，且这是那份预注册最有用的一节**：两处变化**无法分离**，
但**方向可以事先判定**——病例 +8.2%（↑）、对照 −1.7%（N_eff 里只占 1.6%，≈0）、
排除规则变宽使对照更干净（↑ 或 ≈0）。**三项没有一项指向"更弱"。**
由此得到不对称判读：**转移成功 = 部分被混杂**（不得说成功效所致）；
**转移失败 = 两处变化都解释不了**（不得拿定义差异当开脱，纪律第 9 条）。
实际结果是成功，故按前者写。

### 数据状态

`r13/finngen_R13_C3_MELANOMA_SKIN_WIDE.gz`（804 MB，gzip OK，**已扫描**）
`r13/finngen_R13_C3_HEPATOCELLU_CARC_WIDE.gz`（792 MB，gzip OK，**已扫描**）
`r13_manifest.tsv`（820 KB）
`r13_extracts/*.tsv`（四份提取缓存：R13/R12 × melanoma/HCC，各约 0.5 MB）
—— **step99 重跑会直接读缓存，不再扫 3.2 GB，秒级完成。**

R13 HCC 的对应数字（S25 的 C3/C4 格用到）：`C3_HEPATOCELLU_CARC_WIDE` =
**1,070 / 372,159**，对 R12 的 947 / 378,749 为病例 +13.0%。

---

## 五、下一轮的工作顺序

**v4 原列的 1–4 项已全部完成**，"三份文件各自为源"也已按方案 (a) 回填解决（见 §三）。
剩下的：

1. **Step 96 的自我 ground truth（6/10）三份稿子都没写。** 决定进不进正文。
   进的话注意它的限定：样本小、公认基因清单是我们自己定的（跑前定死但非预注册）。
   若要进，**回填顺序仍是先全文源、再 GB/NC**——这是本轮刚恢复的纪律，别再破。
2. **GB 摘要 370 词，超 ~350 目标 20 词。** 能砍的都砍过了；再砍会动到限定语。
   投稿时若期刊硬性卡 350，优先删"Down-sampling shows the loss falls unevenly…"
   那句的百分比细节，**不要删限定语**。
3. GB 正文 6,005 词，已进惯例区间但仍偏下沿，可加回：
   coloc 窗口/先验敏感性的完整表述、匹配功效模拟的校准细节、患者部分的三种重复样本处理。
4. §六待核实清单仍未清（Chen et al. 的两个 FDR、DepMap release 号、6 处 `[ref]`）。

---

## 六、⚠ 待核实清单（**未核实不得写入正文**）

| 项 | 状态 |
|---|---|
| Chen et al. *Nature* 2025;642:191–200 的 TPI1 gene-level FDR（约 0.999 / 0.480）| **论文身份已核实**（PMID 40140585），**两个 FDR 未核实**。Springer 静态路径只取到 17 KB 附件；ORCS 发布包无索引，无法定位 SCREEN_ID。详见 `MANUSCRIPT_worknotes.md` |
| DepMap 的具体 release 号 | 须与脚本核对后定稿（`REFERENCES.md` 条目 34）|
| 正文 6 处 `[ref]` 标记 | 对应 `REFERENCES.md` 条目 35–40，装配前替换为期刊编号格式 |

**BioGRID ORCS API**：用户提供的 access key 在**主库有效、ORCS 被拒**——
ORCS 需单独开通。**但不必再追**：批量下载（免 key、MIT 许可）已给出更强的全库统计。

---

## 七、技术坑（v3 第六节之外的新增）

20. **GEO 的 h5ad 里没有的临床标注，往往在 GSM 层的 `!Sample_title` / `characteristics` 里**。
    只看 `obs` 就断言"无标注"会出错（预注册 S20 §7 犯过一次）。
21. **NCBI FTP 大文件容易断在中途** —— 用 `curl -L -C - --retry 10 --retry-all-errors` 续传。
22. **BioGRID ORCS 的批量包免注册**（MIT 许可），但**第一次可能下到一半就截断且续传无效**
    （服务器不支持 range）——**必须 `gzip -t` 校验，删掉重下**。发布包**无索引文件**，
    SCREEN_ID → 论文/细胞类型的映射不在其中。
23. **FinnGen R13 的两个坑**：manifest 路径要带 `summary_stats/` 前缀（根目录 404）；
    该桶在默认 HTTP/2 下报 curl exit 35（TLS），**必须加 `--http1.1`**。
24. **GWAS Catalog 的 associations 端点只返回 rsID，不给坐标**；
    `efoTraits/{id}/associations` 已 404，可用的是
    `associations/search/findByEfoTrait?efoTrait=<名称>`（用疾病名，不是 EFO ID）。

---

## 八、方法论纪律（v3 第七节之外的新增）

13. **撤回一条主张时必须全文搜其同义表述，且范围包含 `figures/*.py` 的字符串字面量与图注。**
    Step 87–88 共抓到 5 处残留，**3 处在图里**。
14. **撤回一个统计量时，产出它的结果表列名也要标记。**
    `63c_patient_locked.tsv` 的列就叫 `binom_p`，只要有人再拿它画图，撤回的量就会从数据侧复活。
15. **任何"参照系"用之前先看它在背景上的覆盖率。**
    覆盖率异常低说明是映射错了，不是这个疾病的已知位点本来就少
    （Step 94 的 lung 0.36% / colorectal 0.89% vs melanoma 10.5%）。
16. **审稿人给的数字也要核实过才能进正文。**
    Step 95 差点把未核实的 CRISPR FDR 写进 Limitations；Soskic 的四个数字则逐位核对后确认无误。
17. **不要从 manifest 的一行描述推断 endpoint 定义**（Step 98）。
    要查 Risteys 的逐项编码定义。推断出来的结论会把"对照排除规则变了"说成"病例口径变了"。
    → 这一条已升格为**正文的第八条诊断**。
18. **一张结果表的 outcome 列是对着哪个结局算的，用之前必须查产出它的脚本**（Step 99）。
    `92c_locus_annotated.tsv` 的 `beta_out/se_out` 是对着 **meta** 算的，不是 R12；
    拿它当 R12 基线，就在无意中同时换了 release 和数据集。
    **症状是"功效上升而显著数下降"这类方向反常**——出现方向反常时先查基线身份，别急着解释。
19. **两个 release 的结果如果逐位相同，先判断它是不是算术的必然。**
    Fisher P 只依赖那张 2×2 表，位点计数没变就必然给同一个 P。
    不是 bug，但**也不是第二次独立确认**，两种误读都要防。
20. **计数相同不等于集合相同。** R12 与 R13 的名义显著记录都是 275 条，
    但只共享 224 条。**报告 Jaccard，不要只报 n。**

---

## 九、各主张的当前强度

### 最强

| 主张 | 证据 |
|---|---|
| **位点归属现象** | 独立位点级 4.09×（P=0.028）· 五个嵌套 release 无新位点基因 · **换暴露资源 4.44×（P=3.7e-11）** · **六格网格全部同向 + 错配阴性对照干净** · **R13 转移测试：名单逐字不变而 3,434 条 z 无一相同（S25，判读 A）** |
| **可工具化率的三层口径** | 3/28 instrumentable · 2/28 analysable · 1/28 nominal（**不得压成一个数**）|
| **区室归属** | 6.7× / 2.4× / 方向一致，三队列，深度匹配后仍在 |
| **功效-稳定性** | 已知 vs 新位点 85.8% vs 22.8%；2,506 vs 8,771 例；12/12 预注册命中 |
| **暴露/结局双轴** | 换 300 倍暴露资源：黑色素瘤 7→30，HCC 几乎不动 —— **暴露定工具变量数，结局定其中多少能显著** |

### 中等（须带限定）

| 主张 | 限定 |
|---|---|
| HCC 泛化 | **主检验 P=0.110 未过**；只能写"方向一致的支持性证据" |
| 网格普适性 | 只有 **2 疾病 × 2 资源**，六格中三格未达 P<0.05，HCC 两层不独立 |
| 自我 ground truth | 6/10 归对；样本小，公认基因清单是我们自己定的（跑前定死但非预注册）|

### 弱（勿再拔高）

| 主张 | 现状 |
|---|---|
| 患者糖酵解应答分层 | 三队列**无任何一支被确认两次**；发现队列两支都显著，两次后续互为镜像 |
| TPI1 作为候选 | FDR=0.119 · PP.H4=0.51 依赖先验（8/16）· 无 credible set · **ORCS 命中率 42.7% 属核心必需带** · Open Targets 里黑色素瘤不在其前列 |

**TPI1 的表述定死**：只能是 "instrument-visibility worked example"，
**不得**隐含 dynamic genetic regulation（Step 95 的交互检验已否定），
**不得**称最佳候选。证据阶梯见 GB 版 Discussion。

---

## 十、文件地图（本轮新增）

**预注册**：`PREREG_hcc_generalisation.md`(S20) · `PREREG_hcc_part2_generalisation.md`(S21) ·
`PREREG_exposure_resource.md`(S22) · `PREREG_generality_grid.md`(S24) ·
**`PREREG_r13_transfer.md`(S25)**
**补充材料**：`SUPP_literature_audit.md`(S23) · `MANUSCRIPT_worknotes.md`（工作笔记，不进正文）
**脚本**：`step84`–`step98`、**`step99_r13_transfer.py`**、`manuscript/assemble.py`
**结果表**：`84a`–`85e`（HCC）· `87a`–`87c`（患者）· `91a`–`91e`（文献）·
`92a`–`92e`（eQTLGen）· `94d`–`94f`（网格）· `96a`（自我检验）· `97a`（ORCS）·
**`99a`–`99c` + `99_console.log`（R13 转移测试）**
**数据**：`eqtlgen/` · `hcc/` · `r13/` · `orcs/`（752 MB）· `cancer_extracts_eqtlgen/` ·
**`r13_extracts/`（工具变量位置的提取缓存，重跑 step99 秒级完成，已 gitignore）**

### git 状态（本轮已解决）

`.git/index.lock` **已不存在**，无需再 `rm`。本轮三个 commit：

| commit | 内容 |
|---|---|
| `69af5ad` | **S25 预注册 + step99 脚本，在任何 R13 结果产出之前提交**（时间戳即证据）|
| `5be2afb` | Step 94–98 全部文件 |
| `3598647` | Step 99 结果 + 第八条诊断 + ORCS 压缩 |

⚠ **v4 原来建议的 `git add -A` 会提交约 300 MB 的 LD 矩阵**（`regions/*.vcor1` 四个文件
合计 292 MB，外加 plink 的 pgen/pvar/psam）。`.gitignore` 本轮已补上
SMR/plink 工作树、`soskic/`、`orcs/`、`icb_bulk/`、`cancer_extracts_eqtlgen/`、
`r13/`、`r13_extracts/` 等目录。补完后未跟踪总量 3.3 MB，`git add -A` 才是安全的。
仓库现 490 余文件，无 gz/rds/h5ad/bam/fastq/vcf。
