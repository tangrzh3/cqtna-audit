# S40 前提一：合格外部队列的筛查

**日期**：2026-08-21
**规则来源**：`EXTERNAL_VALIDATION_PROTOCOL_v9_3.md` §1.3、§2、§7（S40 全文采纳，冻结）
**产物**：`132b_exclusion_roster.tsv`（排除名单）·
`132a_power_gate_curve.tsv`（读法 A）· `133a_projection_reading_b.tsv`（读法 B）·
`134a_candidate_evaluation.tsv`（逐队列判定）· `S40_POWER_GATE_SPEC.md`（实现说明）

> **一句话**：排除名单已经查全了，**没有一个已发表的皮肤黑色素瘤 GWAS 资源在名单之外**；
> 名单外**符合 §2.2 表型要求**的最大队列，是 MVP 的 5,364 例组织学确诊侵袭性黑色素瘤，
> 其 N_eff 只有发现集的 **0.43 倍**（MVP 另有两个更大的 phecode 口径，但表型不合格）。§2 的两条取数路径里，第二条（consortium
> leave-cohort-out）**在结构上已经作废**。

---

## 1. 排除名单已经查全（26 行，逐条有出处）

§1.3 要求与三个来源的**全部贡献队列**无 participant 重叠。三个来源的队列谱系：

| 来源 | 贡献队列 | 出处 |
|---|---|---|
| **FinnGen R12/R13** | 全部芬兰生物库 | 主结局本身 |
| **Rashkin GCST90011809** | **UK Biobank + GERA**（Kaiser Permanente 北加州）| GWAS Catalog study record 的 `cohort` 字段：`UKB\|GERA` |
| **Landi 2020 参照名单** | **21 个 confirmed 集 + 2 个 self-report 集** | Landi 2020 Supplementary Table 1（已取得原表） |

⚠ **这份名单查全到"数据集"这一层，不是"招募中心"这一层。**
"全部芬兰生物库"是一个类别而不是一份枚举；GenoMEL Phase 1/2 本身是个联盟，
它的成员中心没有在 Landi 的 Supplementary Table 1 里逐一列出。
**§1.3 要的保管方签署谱系表，正是用来补这一层的**——本项目补不了。

Landi 的 23 个数据集逐一列在 `132b_exclusion_roster.tsv`。其中三个值得单独指出，
因为它们是"看上去像独立队列、其实已在名单里"的那类：

- **Michigan**（1,198 例 / 26,211 对照）＝ **Michigan Genomics Initiative**。
  MGI 是常被当作独立验证集的美国医院生物库，**在这里不能用**。
- **NCI CPS-II + PLCO + Rose**（171 例）—— 三个美国队列研究已被并入。
- **UK Biobank 在名单里出现两次**（Landi confirmed 3,499 例 + Landi self-report 1,802 例），
  并且**再一次**以 Rashkin 的 `cohort` 标签出现。

### 1.1 一个结构性后果：§2 的第二条取数路径作废

§2 允许两种取数方式：

> (i) a de novo GWAS in an independent biobank, or
> (ii) a consortium-generated **leave-cohort-out** meta-analysis excluding all
> overlapping discovery/reference cohorts.

**(ii) 已经没有东西可留**：参照名单来自 Landi，而 Landi 的 21 个 confirmed 集
**就是**全世界已确诊黑色素瘤 GWAS 的主体（30,134 例中的全部）。
把重叠队列全部剔除之后，consortium 手里剩下的是空集。

**因此 S40 只剩 (i) 一条路：在一个从未参与过上述任何一项的生物库里新做一个 GWAS。**

---

## 2. 名单外还剩什么（按可查证程度排序）

GWAS Catalog 的两个 trait（`MONDO_0005105` melanoma、`MONDO_0005012` cutaneous
melanoma）合计 **72 项唯一研究**，逐条查过 `cohort` 字段。

⚠ **`cohort` 字段只对较新的条目填写**：2020 年以前的条目大多是 `None`，
那是"未标注"而不是"无排除队列"，必须回到出版物去判。所以分两段说：

**2020 年之后**：带排除标签（UKB / FinnGen / QSKIN / Q-MEGA / MIA / GERA）的之外，
**只剩 MVP 的 7 项**（其中欧洲祖源 2 项）与 `GCST90132200`（Pflugfelder 2022，
标注 `NR`，仅 556 人）。**MVP 是唯一有实质体量的例外。**

**2020 年以前**：有分量的条目基本都已被 Landi 吸收——
例如 `GCST004142`（Ransohoff 2017，4,842 例 / **286,565 对照**）
的对照数与 Landi 的 23andMe 行（4,824 / **286,565**）相同；
`GCST003061`（Law 2015）是 Landi 的前身 meta。
未被吸收的只有极小的几项（Song 2014 共 494 例、Teerlink 2011 共 156 例高危家系）。

**结论：除 MVP 外，没有已发表的皮肤黑色素瘤 GWAS 资源可以直接拿来当外部验证集。**

| 候选 | 口径 | 病例 | 对照 | N_eff | R = N_eff/49,337 | 证据 |
|---|---|---:|---:|---:|---:|---|
| **MVP** | PheCode 172.1 "dx **or hx**" | 45,274 | 389,597 | 162,242 | **3.29** | ✅ GCST90475577，**全 P 值公开** |
| **MVP** | PheCode 172.11 | 10,468 | 436,283 | 40,891 | **0.83** | ✅ GCST90475578，**全 P 值公开** |
| **MVP** | 组织学确诊**侵袭性** | 5,364 | ~436,283 | 21,195 | **0.43** | ✅ PMID 39853431（须新做 GWAS）|
| All of Us | — | 未知 | — | — | — | ⚠ 须向保管方索取 |
| Estonian Biobank | — | 未知 | — | — | — | ⚠ 须向保管方索取 |
| deCODE / 冰岛 | — | 未知 | — | — | — | ⚠ 商业授权 |
| Copenhagen Hospital Biobank / DBDS | — | 未知 | — | — | — | ⚠ 须向保管方索取 |
| HUNT（挪威） | — | 未知 | — | — | — | ⚠ 须向保管方索取 |
| BioVU · MGB · Penn · Geisinger · BioMe | — | 各自很小 | — | — | — | ⚠ 须逐个索取并相互查重 |
| Our Future Health（英国） | — | 未知 | — | — | — | ⚠ **与 UKB 的参与者重叠须先排查** |
| BioBank Japan · 台湾 · CKB | — | — | — | — | — | ❌ §2.3 把非欧洲祖源排除出主分析 |

### 2.1 MVP 是唯一现实的候选，但卡在**表型**上，不是卡在准入上

MVP 的两个黑色素瘤 GWAS **全 P 值已经公开**（Verma 等 2024 *Science*，PMID 39024449），
今天就能跑，不需要任何数据传输协议。问题在 §2.2：

> The primary phenotype is invasive cutaneous melanoma or a definition
> **demonstrably compatible** with it.

两个公开口径**都是 phecode**：

- `GCST90475577` = PheCode **172.1**"Melanomas of skin, dx **or hx**"，45,274 例
- `GCST90475578` = PheCode **172.11**"Melanomas of skin"，10,468 例

而 Wheless 等 2025（*Arch Dermatol Res*, PMID 39853431）在 **MVP 自己身上**做过这件事：

| MVP 口径 | 病例 | 全基因组显著变异 |
|---|---:|---:|
| phecode ≥ 2 次 | 45,665 | 20,457 |
| 组织学确诊**侵袭性** | **5,364** | 2,582 |
| 组织学确诊原位 | 4,792 | 1,989 |

> "Most of the variants identified in the phecode cohort **did not replicate**
> in the histologically-confirmed cohorts."

**所以 45,274 那个数字不能拿来充功效**——否定它的证据是 MVP 研究者自己发表的。
§2.2 意义上可用的是 **5,364 例**那一档，而那一档**还没有做成 GWAS**，
要做就是 §2 的路径 (i)：在一个独立生物库里新做一个。

⚠ **两条必须写进谱系表、但都不构成排除理由的事**：

1. **MVP 与 GERA 的重叠不为零**（退伍军人同时是 Kaiser 会员是可能的）。
   §1.3 要的是保管方签署的重叠表，不是我们的估计。
2. **Landi MT 与 Brown K（NCI）是 MVP 黑色素瘤那篇的合作者**。
   这是**作者重叠，不是参与者重叠**，§1.3 管的是后者。
   写下来是为了防止日后有人把它误当成独立性问题——也防止有人以为我们没注意到。

---

## 3. 功效闸门要求什么（`step132`、`step133`）

§7 的判据（已冻结）：在 **f = 0.50** 的 effect-retention 下，
≥80% 的模拟要产生 **≥8 个显著 bounded locus 且通过主 Fisher 终点**。

发现集：**N_eff = 49,337**（12,530 例 / 789,099 对照），实得 8 个显著位点、
其中 4 个 known、fold 4.96、单侧 Fisher P = 0.0048。

### 3.1 ★★ §7 没有规定外部效应估计怎么抽，两种读法差一倍

§7 规定了 se 怎么投影、规定了三个 retention 因子、规定了跑一万次模拟，
**但没有一句话说 beta_ext 是怎么来的**。两种自然读法：

| | 模型 | 随机性来源 |
|---|---|---|
| **读法 A** | `beta_ext ~ Normal(f · beta_disc, se_ext²)` | 抽样噪声 + 覆盖率 |
| **读法 B** | `beta_ext = f · beta_disc`（确定性）| **只有**覆盖率 |

**读法 A 不能用，因为它通不过自己的自校准。**
在 f = 1.00、R = 1.00 下——外部队列与发现集**效应量相同、样本量相同**——
读法 A 预测**中位 102 个显著位点**。本研究实得 **8 个**。

原因是它把 `beta_disc` 当成真值又在上面加一层噪声，于是
Var(z_ext) = Var(z_disc) + 1。发现集里绝大多数记录本来就是噪声，
这一步把噪声的方差翻倍，尾部被人为撑开。

**读法 B 逐项复现发现集**（g = f·√R = 1 时 ST = 8、SK = 4、fold 4.962、
P = 0.0048407，见 `step133` 的自校准断言）。**因此下面用读法 B。**

⚠ 但读法 B 也不是一个合格的功效计算：它没有抽样方差，
"≥80% 的模拟通过"只能靠覆盖率的波动来兑现。
**统计上正确的做法是第三种**——先从发现集的 z 估出真实 non-centrality 的分布，
再从那个分布抽——**而那是一个设计选择，S40 §2.1 不允许本项目单方面做。**

**➜ 这一条要带回给方案作者：§7 需要补一句话，规定 beta_ext 怎么抽。**

### 3.2 ★ 判定不是 N_eff 的单调函数：队列可以**大到不合格**

读法 B 是确定性的，所以整条轨迹可以逐格看清楚。判定只取决于
**g = f·√R** 这一个量（全覆盖）：

| g | 显著位点 | **fold** | Fisher P | 判定 |
|---:|---:|---:|---:|---|
| 0.354 – 0.919 | 2 – 6 | 9.92 – 4.96 | — | ❌ 不足 8 个位点 |
| **1.000** | **8** | **4.962** | **0.00484** | ✅ ← **恰好是发现集本身** |
| 1.225 | 44 | 3.158 | 3.0×10⁻⁵ | ✅ |
| 1.414 | 120 | 2.233 | 4.7×10⁻⁶ | ✅ ← **P 的最低点** |
| 1.732 | 257 | 1.352 | 0.0117 | ✅ |
| 2.000 | 326 | 1.248 | 0.0232 | ✅ |
| 2.828 | 463 | 1.136 | 0.0444 | ✅ 勉强 |
| **3.162** | 494 | 1.105 | **0.0736** | ❌ |
| 4.000 | 538 | 1.051 | 0.223 | ❌ |

**fold 从 4.96 单调衰减到 1.05；Fisher P 是 U 形的**，在 g ≈ 1.41 触底
（4.7×10⁻⁶），之后回升，在 g ≈ 2.9 与 3.2 之间越过 0.05。

原因是结构性的，与读法无关：主终点是"显著位点集"相对"全部可检验位点"的富集。
功效越大，显著集越接近背景（g = 4 时 538 个显著位点 / 655 个背景位点），
**fold 必然趋向 1**，Fisher 随之失去可判读的对比。

**➜ 合格窗口是 g ∈ [1.00, ≈2.9]，即 R ∈ [1/f², ≈8.4/f²]。**
f = 0.50 时是 **R ∈ [4, 34]**。上界在实践中够宽，
但**"下界之上一律通过"这个直觉是错的**，而 §6 的写法看不出这件事。

**➜ 这也要带回给作者：§6 第 3 条（fold > 1 且 P < 0.05）在高功效下会自我否定。**

### 3.3 门槛（读法 B，全覆盖）

| f | 通过所需 | N_eff | 若对照池极大，约等于病例数 |
|---:|---:|---:|---:|
| **0.50**（判定用） | **R ≥ 4.0**（g ≥ 1.000）| **197,348** | **≈ 49,300** |
| 0.75 | R ≥ 2.0（g = 1.061）| 98,674 | ≈ 24,700 |
| 1.00 | R ≥ 1.0（g = 1.000）| 49,337 | ≈ 12,300 |

判据恰好是 **g = f·√R ≥ 1**：外部队列必须把 retention 的损失**完全用样本量补回来**。
f = 0.50 意味着 **R ≥ 1/f² = 4**。

---

## 4. 结论

**没有合格的外部队列，而且缺口不是边际的。**

把闸门直接套在每个候选队列自己的 N_eff 上（`step134`，f = 0.50，全覆盖）：

| 队列 | R | g | 显著位点 | fold | Fisher P | 判定 |
|---|---:|---:|---:|---:|---:|---|
| 门槛 | 4.00 | 1.000 | 8 | 4.962 | 0.00484 | — |
| MVP PheCode 172.1（45,274 例）| 3.29 | 0.907 | **6** | 4.962 | 0.0157 | ❌ **且表型不合 §2.2** |
| MVP PheCode 172.11（10,468 例）| 0.83 | 0.455 | 2 | 9.924 | 0.010 | ❌ |
| MVP 组织学确诊侵袭性（5,364 例）| 0.43 | 0.328 | 2 | 9.924 | 0.010 | ❌ |

⚠ **值得说清楚的是它是怎么失败的**：MVP 最宽那个口径**不是在富集上失败**——
fold 仍然是 4.96、P = 0.0157——**它是在信息闸门上失败**：
6 个显著位点，而 §4 要求至少 8 个。**证据方向没错，量不够判读。**

（自洽性核对：发现集自己在 f = 1.00、R = 1.00 下返回 8 个位点、fold 4.962、
P = 0.0048407，与实测逐位相同。）

用 MVP 的对照池（436,283）反解：§2.2 合规的**组织学确诊侵袭性**队列
需要 **55,627 例**才能达到 R = 4；对照池无限时也要 **49,339 例**。
MVP 现有 5,364 例，**差 9.2–10.4 倍**。
作为参照：**史上最大的黑色素瘤 meta（Landi 2020）确诊病例总共 30,134 例，而且全部在排除名单里。**

### 4.1 三条路，作者选一条

1. **接受 §7 自己的降级条款**——"Otherwise the cohort may be analysed only as a
   **labelled pilot**"。用 MVP 的公开全 P 值今天就能跑，
   但**只能写成 pilot，不能写成 external validation**，也不能用来升级任何主张。
2. **把 §7 退回给方案作者**，请其补上 §3.1 的抽样规定、并重新审视 §3.2 的高功效自我否定。
   在那之前，**任何 N_eff 门槛都不是确定的数**。
3. **新建一个多生物库 meta**（§2 的路径 (i)）：MVP 组织学确诊 5,364 例
   ＋ All of Us ＋ EstBB ＋ 北欧登记 ＋ 美国医院生物库。
   ⚠ 即使乐观地凑到 15,000–25,000 例，R 也只有 1.2–2.0，**仍然通不过 f = 0.50**。

**➜ 本项目的建议：走第 2 条，同时把第 1 条作为可交付的下限。**
把"S40 已采纳但无合格队列"改写成一句更具体、更可核验的话：
**门槛是 R ≥ 4；现存最好的独立资源是 R = 0.43；缺口约 10 倍。**

### 4.2 正文措辞可以更新（仍不得称 pre-registered）

S40 §4 现在的措辞说"no qualifying cohort has been secured"。
现在可以说得更硬、也更有信息量：

> An external validation protocol has been adopted and its design frozen
> (Supplementary S40). Screening every contributing cohort of FinnGen R12/R13,
> Rashkin GCST90011809 and the 23 datasets of the Landi reference meta-analysis
> leaves **no published cutaneous melanoma GWAS resource eligible**, and the
> protocol's own power gate requires an effective sample size four times that of
> the discovery meta-analysis. The largest eligible resource we could identify —
> 5,364 histologically confirmed invasive cases in the Million Veteran Program —
> reaches 0.43 times it. The protocol is therefore reported as a commitment,
> with the size of the gap stated, rather than as evidence.

⚠ **仍然不得写 pre-registered**（S40 §3：还差公开注册取得 DOI）。

---

## 5. 交给作者去做的事（本项目做不了的）

1. 向 **MVP / VA** 申请：以 Wheless 2025 的组织学确诊口径新做一个 GWAS
   （联系点：Wheless L, Hartman RI；注意 Landi MT 与 Brown K 是该文合作者——
   **作者重叠不是参与者重叠**，但要写进谱系表）；
2. 向 **All of Us、Estonian Biobank、deCODE、Copenhagen Hospital Biobank/DBDS、HUNT**
   索取组织学确诊侵袭性黑色素瘤的病例数与队列谱系表；
3. **Our Future Health 必须先查与 UKB 的参与者重叠**，不能默认它是新队列；
4. 决定走 §4.1 的哪一条路——这是科学判断，不是本项目能替作者做的。
