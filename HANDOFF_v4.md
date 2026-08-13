# 项目交接文档 v4

**项目**：CD4⁺T 细胞活化动态 eQTL → 靶点提名的审计（黑色素瘤 + HCC，两套暴露资源）
**最后更新**：2026-08-13（Step 98 后）
**取代**：`HANDOFF_v3.md`（其第六节技术坑、第七节方法论纪律仍有效，本文件继承并扩充）

> **新会话阅读顺序**
> 1. 本文件
> 2. 两版投稿稿：`manuscript/MANUSCRIPT_GB.md`、`manuscript/MANUSCRIPT_NC.md`
> 3. 全文源 `manuscript/MANUSCRIPT_v2_dual_thread.md`（**唯一正文源**，两版都是派生物）
> 4. `FINDINGS_step5_pigmentation.md` 的 Step 84 起
> 5. **六份预注册**（S9/S18/S20/S21/S22/S24）+ `SUPP_attempt_timeline.md`

---

## 一、当前定位

**一条方法学主线 + 一个被严格限定的 worked example。**（不再是"双主线"——
审稿人 R2 指出第二条线证据不足以对等，已接受。）

主线：**靶点提名是结局 GWAS 的性质，不是暴露侧的性质。** 支撑它的是一张
**疾病 × eQTL 资源的交叉网格**，六格全部同向，配一个干净的错配位点阴性对照。

worked example：TPI1，定性为 **instrument-visibility 的示例**，不承担任何
动态遗传调控或靶点主张。

---

## 二、这一轮（Step 84–98）新增了什么

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

⚠ **全文源 `MANUSCRIPT_v2_dual_thread.md` 一个字没动**，两版都是派生物。
`assemble.py` 负责把 Methods 与 References 装配进全文版（生成 `MANUSCRIPT_assembled.md`，
**永不手改**）。`build_gb.py` 是早期的减法脚本，**已被手写 GB 版取代，可删**。

---

## 四、★ FinnGen R13：endpoint 映射已定，分析未做

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
第七份预注册的判读表必须事先写死：两者混杂时怎么读，以及
"若结果与 R8–R12 轨迹外推不符，是否允许归因于定义差异"（按纪律第 9 条，**不允许**当开脱）。

### 数据状态

`r13/finngen_R13_C3_MELANOMA_SKIN_WIDE.gz`（804 MB，gzip OK）
`r13/finngen_R13_C3_HEPATOCELLU_CARC_WIDE.gz`（792 MB，gzip OK）
`r13_manifest.tsv`（820 KB）

---

## 五、下一轮的工作顺序（建议）

1. **写第七份预注册（R13 转移测试）并跑** —— 数据已就位，纯计算。
   判读表按"同资源、更新 outcome dataset 的转移测试"设计，**不是**第六个功效点。
2. **git 解锁并提交** Step 94–98 全部文件：
   `rm -f .git/index.lock` 然后 `git add -A && git commit`
3. **ORCS 那一节从一节压成一段**（现在给了它一整节，偏重；它是边界性观察 + 一条旁证，
   不是独立发现）。⚠ 压缩时**不得**改成"CRISPR 数据支持我们的结论"——
   正确措辞是 "a knockout cannot isolate a CD4-specific role"。
4. **R13 定义映射写进 §2.2 轨迹段的限定**，并考虑立为**第八条诊断**：
   *跨 release 比较之前，先核对 endpoint 定义是否同一。*
   （任何做"随 GWAS 更新看名单变化"的研究，不核对定义就比，测到的"不稳定"
   里会混进纯粹的表型漂移。这条比现有七条里的好几条更普适。）
5. GB 版还有 1,000–3,000 词余量，可加回：coloc 窗口/先验敏感性的完整表述、
   匹配功效模拟的校准细节、患者部分的三种重复样本处理。

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

---

## 九、各主张的当前强度

### 最强

| 主张 | 证据 |
|---|---|
| **位点归属现象** | 独立位点级 4.09×（P=0.028）· 五个嵌套 release 无新位点基因 · **换暴露资源 4.44×（P=3.7e-11）** · **六格网格全部同向 + 错配阴性对照干净** |
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
`PREREG_exposure_resource.md`(S22) · `PREREG_generality_grid.md`(S24)
**补充材料**：`SUPP_literature_audit.md`(S23) · `MANUSCRIPT_worknotes.md`（工作笔记，不进正文）
**脚本**：`step84`–`step98`、`manuscript/assemble.py`
**结果表**：`84a`–`85e`（HCC）· `87a`–`87c`（患者）· `91a`–`91e`（文献）·
`92a`–`92e`（eQTLGen）· `94d`–`94f`（网格）· `96a`（自我检验）· `97a`（ORCS）
**数据**：`eqtlgen/` · `hcc/` · `r13/` · `orcs/`（752 MB）· `cancer_extracts_eqtlgen/`

⚠ **git 有一个未释放的 `.git/index.lock`**，Step 94–98 的文件**尚未提交**。
仓库首个提交（402 文件）与第二个提交已在，安全检查通过（无 rds/h5ad/bam/fastq/vcf）。
