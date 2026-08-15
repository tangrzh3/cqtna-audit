# Findings — MANUSCRIPT_NC recent MR benchmark

Key manuscript observations and verified comparator findings will be recorded here.

## Manuscript audit
- `MANUSCRIPT_NC.md` is a compact Nature Communications-formatted version: 140-word abstract and ~2,444-word main text; 7 proposed display items.
- The revised draft now correctly distinguishes dynamic expression, time-specific eQTL significance, and genotype-by-pseudotime interaction for TPI1.
- Central empirical contribution: five nested FinnGen releases; two diseases; two exposure resources; crossed grid with mismatched disease-locus control; outcome-locus enrichment and candidate-list instability.
- Strongest remaining methodological limitation stated in the draft: coloc.abf single-causal assumption; external-LD SuSiE attempt failed by over-splitting; key novel loci cannot be fine-mapped.
- Main biological arm is deliberately negative/bounded: TPI1 is nominal, prior-sensitive, unresolved genetically, tumour-dominant in tissue, and patient programme not consistently replicated.
- Current title says nomination is set by outcome GWAS "not by exposure", but the manuscript itself shows exposure power sets instrument availability; title is therefore sharper than the evidence permits.

## Claim-calibration issues
- The abstract should label TPI1 as a nominal association (FDR 0.119), not simply “an association”.
- “This does not depend on outcome dataset, disease, or exposure resource” is too absolute when only 3/6 crossed-grid enrichment directions reach P<0.05. The defensible statement is that the direction recurs across tested settings.
- “The two diseases behave alike once power is equalised” shows compatibility inside a prediction interval, not statistical equivalence.
- A defensible headline is: exposure power governs instrument availability; outcome architecture and power strongly govern which instruments cross the nomination threshold.

## Initial comparator pattern
- Recent high-impact papers rarely publish summary-statistic MR as a self-contained endpoint. MR is typically embedded in a new molecular resource, independent replication, prospective phenotyping, multi-omic convergence, or functional experimentation.
- A 2025 Nature Medicine HIV multi-omics study enrolled 1,342 participants, measured five omics layers, mapped 5,962 molecular QTLs, and used MR to connect molecular features with ex-vivo cytokine phenotypes. It included instrument validation filtering, bidirectional/sensitivity analyses, public summary statistics, a web app, and code. Its evidential advantage is a newly generated cohort/resource plus triangulation.
- A 2025 Nature Communications pan-disease proteogenomic study combined 7,626 healthy participants with 28,064 patients, 2,901 proteins across 42 disease states, 25,987 pQTLs, and 110 reported high-confidence causal proteins across 21 diseases. Its advantage is scale, breadth, and connection of causal inference to measured disease-state proteomes.

## Verified recent benchmark designs
- European Heart Journal (14 Jan 2026), Meena et al.: proteome-wide cis-MR (2,923 proteins; pQTL GWAS up to ~54,219), FDR, weak-instrument filtering, reverse-direction/Steiger checks, Bayesian colocalisation, network MR mediation, UK Biobank observational analyses, and PheWAS. Only 4 proteins survived the MR-plus-colocalisation chain for cardiovascular mediation. This is a useful example of progressively narrowing claims rather than treating every MR hit as a target.
- Circulation (2025), AAA study: MR across proteins, metabolites, and lipoproteins was followed by three hypertriglyceridaemic mouse models, RNA-seq/cell mechanistic work, local mediator rescue, and an antisense-oligonucleotide intervention. It supports a mechanistic/therapeutic claim because perturbation and rescue sit downstream of MR.
- Circulation (16 Feb 2026), VTE proteomics: prospective discovery, independent cohort replication on a different proteomics platform, MR with an explicit tiered significance rule, colocalisation, and sensitivity MR restricted to replicated cis-pQTL instruments. The key lesson is that both the phenotype association and the instrument were replicated before causal prioritisation.
- Nature Genetics (15 Apr 2026) “Challenges and future directions for Mendelian randomization” is a current methodological benchmark, not an original comparator. It supports emphasizing empirical calibration, robustness to instrument/LD assumptions, population matching, and triangulation rather than method-count accumulation.

## Manuscript-specific implications from comparators
- The manuscript should adopt a formal evidence ladder, for example: screened record → independent locus → outcome-known/novel annotation → FDR MR → locus-resolved colocalisation → release/resource replication → compartment-compatible biology. “Target” should be reserved for the final tier; TPI1 stops much earlier.
- The paper's genuine differentiator is empirical calibration of nomination instability, not nomination itself. The audit should therefore be the primary thread; TPI1 should be a bounded worked example, not an equal second discovery claim.
- Add a flow diagram/table reporting attrition at every gate and list why each candidate fails. This mirrors the narrowing pipelines of recent high-impact studies and makes the negative result constructive.
- Do not add virtual perturbation as if it were functional validation. It can be a hypothesis-generating consistency analysis, but cannot raise TPI1 to causal target status because perturbation models inherit training-data correlations and do not resolve the genetic attribution problem.
- Public release should include executable code, frozen exposure/outcome accession versions, machine-readable known-locus lists, preregistration timestamps, complete candidate lists at each FinnGen release, and all null/failed analyses. Recent high-impact resource papers commonly ship portals or reusable summary outputs.

## Time-sensitive dataset opportunity
- FinnGen DF13 was publicly released 2 June 2026 with 500,186 participants, 2,755 endpoints and greater longitudinal follow-up. Because the manuscript's core claim is about outcome-power-dependent list instability, DF13 is an unusually valuable prospective update.
- However, DF13 cancer endpoints changed: `C3_XXX_EXALLC` was renamed `C3_XXX_WIDE`, and cancer case/control definitions were standardized. Therefore DF13 can only be called a further nested power level after exact endpoint-definition identity is verified; otherwise it is a combined power-plus-phenotype sensitivity analysis.

## Final review priorities
- P0 claim repair: replace the binary outcome-versus-exposure title; label TPI1 as nominal; replace universal/equivalence language with recurrence/compatibility language; define the estimand as candidate-list composition conditional on exposure instrument selection.
- P0 statistical repair: because first-order single-SNP Wald z equals the outcome z, collapse repeated gene×time records to unique SNP/locus for primary multiplicity control or use hierarchical FDR; retain profile mapping as descriptive. Quantify how duplicate profiles change BH thresholds and list membership.
- P0 locus-resolution repair: do not let unresolved single-causal coloc support target claims. Use matched in-sample LD/multiple-signal results where available; otherwise explicitly classify a locus as unresolved and run prior/window specification curves.
- P1 generalisation: eQTLGen versus dynamic CD4 changes sample size, cell composition and platform simultaneously, so it is exposure-resource sensitivity, not a pure exposure-power experiment. Either simulate/downsample exposure precision with assumptions clearly labelled, or moderate the causal language about exposure power.
- P1 enrichment null: extend matched-background/permutation testing to all disease-resource cells, matching cis-region/gene density, MAF, LD score, instrument strength and outcome coverage. Use pre-outcome or leave-one-source-out known-locus lists to reduce circularity.
- P1 literature audit: add a reproducible PRISMA-style flow, exact searches, dual manual adjudication or a validated sample with inter-rater agreement, and machine-readable coding. Otherwise call it a text-mining survey of reporting, not a systematic audit of practice.
- P1 prospective validation: add FinnGen DF13 after endpoint identity audit and freeze predictions before analysis.
- P2 narrative: TPI1 no longer demonstrates activation-dependent genetic regulation because genotype×pseudotime interactions are null. It demonstrates the distinction among dynamic expression, time-specific significance and dynamic genetic effects; virtual perturbation cannot substitute for functional validation.
- P2 reproducibility: replace the repository DOI placeholder before submission and deposit every release-specific candidate list, null analysis, preregistration timestamp and executable environment.

## Fresh severe review — current 345-line version
- Title/abstract make an absolute claim: nomination is “set by” / “a property of” the outcome GWAS. The Results later concede a two-axis model: exposure power sets instrument availability and changes the number of significant loci, while outcome power controls outcome z and threshold crossing. The headline erases this conditional estimand.
- The claim “This does not depend on the outcome dataset, disease, or exposure resource” is contradicted by the manuscript’s own nuance: only 3/6 crossed combinations have P<0.05; HCC-high enrichment is P=0.110; eQTLGen differs in cell mixture, platform, and sample size; the disease comparison is compatibility with a prediction interval, not equivalence.
- The exact single-instrument identity `z_Wald = beta_out / se_out` is mathematically true for first-order Wald MR, but it does not show that nomination is biologically a property of outcome alone. Exposure selection, allele harmonisation, winner’s curse, instrument strength, profile duplication, locus coverage and the denominator of tested hypotheses all define the candidate universe and FDR threshold.
- Primary unit ambiguity is severe: 3,556 “exposure records” from eight activation profiles, then FDR-significant records, genes and independent loci are interchanged across the story. The manuscript must explicitly define which unit receives multiplicity correction and how repeated SNP/gene/profile records are collapsed.
- The strongest, genuinely interesting result is not a target discovery but an empirical audit: known-outcome loci dominate, novel nominations turn over with outcome power, and MR/SMR can disagree with colocalisation in multi-signal loci. This could be a valuable cautionary/resource paper if framed as calibration rather than an over-totalized theorem.
- The Methods confirm the primary instruments are one top cis-eQTL per gene × profile: 3,556 records but only 2,141 unique variants. The text says BH-FDR was controlled “within each analysis family” without defining those families here. This is insufficient because repeated profiles can map the same outcome z/P to multiple records and materially alter the BH universe; a locus-level enrichment after record-level selection does not repair selection multiplicity.
- Exposure-resource comparison is not an exposure-power intervention: dynamic CD4 eQTLs versus eQTLGen whole blood changes N, cell composition, assay/study design, instrument universe and likely gene coverage together. It supports robustness across two exposure resources, not the pure statement that exposure power has been isolated.
- Disease “power equalisation” uses 200 summary-space simulations and checks whether HCC falls inside melanoma-derived prediction intervals. This shows non-contradiction under a melanoma effect architecture; it does not establish that diseases behave alike or that effective sample size alone explains differences.
- The colocalisation section is admirably candid but structurally limits the strongest claim: exposure-side fine-mapping was not possible; `coloc.abf` is single-signal; proxy-LD SuSiE visibly over-splits. Therefore the manuscript can show that MR/HEIDI agreement is not sufficient at calibration loci, but cannot generally adjudicate shared-versus-distinct causal variants for novel nominations.
- Literature audit says “at most” 12/152 and explicitly counts phrases, not quality. It is a reporting-practice text audit, not evidence that studies failed to perform checks, and cannot substantiate “typically does not” without stronger ascertainment and validation.
- TPI1 biology is largely a chain of well-handled negative/bounded findings, not validation of the methodological theorem. Patient discovery effects do not replicate arm-for-arm; TPI1 is FDR 0.119, prior-sensitive coloc, genetically unresolved and broadly essential. Keeping this long worked example costs narrative focus without raising the central claim’s evidential tier.
- Gene-attribution calibration (6/10 accepted genes, failures at MC1R/OCA2/TYR/CDKN2A regions) is potentially one of the paper’s most important analyses, but “generally accepted causal gene” selection and success rules must be independently curated, blinded to nominations, and sensitivity-tested because the denominator is only ten.
- Data/code availability still contains a repository DOI placeholder. A paper whose central value is auditability is not submission-ready until complete executable provenance, frozen candidate lists, known-locus lists, pre-registrations and null/failed analyses are deposited.
- Independent BH recalculation clarified (and corrected an earlier exploratory PowerShell miscount): record-level BH gives 21 records / 14 SNPs / 10 genes; collapsing to one minimum-P record per unique SNP before BH gives 15 SNPs / 11 genes; one minimum-P record per gene gives 13 genes. Thus the present list is measurably unit-dependent, although not explosively so. The manuscript needs a justified hierarchical/cluster-aware primary multiplicity scheme and full sensitivity table.
- Visual inspection of the proposed primary figures shows professional execution but rhetorical overreach remains embedded in the graphics: the power-stability figure is titled “Outcome GWAS power sets the reproducibility” although it benchmarks recovery against the study’s own unstable full-power list; the known/novel stratification rests on only 6 versus 4 reference genes with very wide intervals. The figure makes a fragile calibration look like a law.
- Figure/readability architecture remains fragmented: `MANUSCRIPT_NC.md` proposes 7 items, while `FIGURES_plan.md` still describes an older 8-figure scheme and different numbering. The current paper cannot yet be assessed as a coherent submission package because legends, numbering, final composite figures and Methods are not assembled into the reviewed file.

## GB versus NC comparison — initial structural findings
- GB is ~6,049 words/602 lines versus NC ~3,253 words/345 lines. GB is much closer to a full main-text manuscript and includes structured abstract, fuller background, more diagnostic steps, and Methods/Supplementary pointers.
- GB preserves the same core absolute claim in both title (“is a property of the outcome GWAS”) and abstract conclusion. It is slightly less rhetorically binary than NC’s “not by the exposure,” but still exceeds the evidence.
- GB improves transparency around the R13 transfer: explicitly separates phenotype-definition change from power and admits partial confounding. This is superior to the compressed NC treatment.
- GB worsens some claim calibration relative to NC: it says “MR significance is a property of the outcome GWAS alone” while immediately conceding that exposure selects the variant; it retains “the two diseases behave alike once power is equalised,” which is not an equivalence result.
- GB’s structured abstract is informative but overloaded: it tries to carry release trajectory, HCC, eQTLGen, colocalisation, list turnover, glycolysis, compartment attribution and literature audit at once.
- GB is superior in evidence-boundary reporting: it explicitly calls the FinnGen calibration empirical rather than external, limits curves to observed power, states nested-release agreement is an upper bound, records failed genotype×pseudotime evidence, distinguishes 3 instrumentable / 2 analysable / 1 nominal, and gives a compact TPI1 evidence ladder.
- GB better exposes the project’s own selection history (Box 1), but this also creates reputational risk: “five genes successively designated lead candidate” under changing criteria can read as uncontrolled researcher degrees of freedom. It must be framed as provenance/audit evidence, not as a substitute for a fixed confirmatory analysis.
- GB contains an additional functional-state/multiome branch and functional failure-mode branch absent from NC. These improve completeness and intellectual honesty but make the manuscript feel like three papers unless subordinate to the eight-diagnostic framework.
- GB has a proper limitations section and a clearer relation-to-guidance section, both important advantages over NC. However it still has multiple `[ref]` placeholders and an HTML comment saying a self-imposed ground-truth decision remains open, so it is not submission-ready.
- Direct choice: GB should be the base manuscript; NC is a useful compression donor. NC’s concision improves pace but removes qualifications precisely where this controversial negative-methods paper needs them.

---

# Findings — 105a blind coding

## Dataset inventory
- `105a_blind_coding.tsv` has 92 rows: 46 C1 rows and 46 C3 rows, representing the same 46-paper set.
- Each criterion has 31 `(none found)` snippets and 15 nonempty candidate snippets; snippets are leads only and must be checked against full text.
- `coder2` and `coder2_note` are initially blank.

## Rule anchors
- C1=1 only for an explicit comparison of this paper's own significant candidates/signals with previously reported loci for the same outcome.
- C3=1 only when the paper reports how its candidate list changes or persists under a different outcome GWAS/power setting.
- Generic instrument power, F statistics, exposure-QTL sample size, validation cohort size, and literature-reference titles are C3=0.
- The first line of the rules explicitly limits coding to passages in the TSV `evidence` column. Full-text research may verify provenance/context but must not introduce an unshown passage as coding evidence.
- Consequently, `(none found)` rows remain 0 unless the evidence itself is genuinely undecidable; a full-text passage discovered outside the evidence field cannot convert them to 1.

## Preliminary evidence-only adjudication
- Clear C3 positives: row 3 (more gene-cancer pairs attributed to greater outcome-GWAS power), row 18 (analysis repeated with another outcome GWAS and candidate/protein results compared), and row 47 (smaller-case outcome GWAS groups explicitly yielded limited ability to identify causal genes versus other COVID-19 outcomes).
- Other nonempty C3 snippets are generic power, exposure-QTL power, instrument strength, non-GWAS validation, or power-enhancing meta-analysis without reporting candidate-list dependence; they do not satisfy C3.
- The nonempty C1 snippets are reference titles, dataset/method descriptions, candidate associations without a prior-locus comparison, cell-type annotation, or comparisons to prior biology/MR studies. No clear C1 positive has yet been identified.

## Boundary-case checks
- PMC12697907 directly states that more gene-cancer pairs were observed for breast/prostate cancer because the outcome GWAS had greater power; this meets the explicit C3 positive example.
- PMC11443760 repeats the protein/gene analysis in another AVB outcome GWAS and reports which results remain meaningful/significant; this meets C3.
- PMC12257706's displayed evidence names PABPC1 and rs1693551 from the AD GWAS used in the analysis. It does not say this own signal overlaps a prior known-locus list, so the displayed passage is C1=0. The full paper has other literature-comparison text, but the rules prohibit importing it into the evidence-only code.
- PMC12885728 validates miRNA findings with a second exposure miRNA-eQTL cohort while reusing breast-cancer outcomes; the displayed limitation sentence does not establish candidate-list dependence on an outcome GWAS, so C3=0.
- PMC13147455 uses multiple COPD GWAS cohorts/meta-analyses, but the displayed snippets discuss power construction and cross-cohort exposure-protein coverage rather than reporting how a candidate list changes with the outcome GWAS; C3=0.
- PMC12626534 and PMC13166711 discuss proteomic/RNA-seq validation power, not outcome-GWAS-dependent candidate lists; C3=0.
- PMC11867302 explicitly ties the smaller case counts of two COVID-19 outcome datasets to reduced ability to identify causal genes compared with other COVID-19 outcomes; this is C3=1 under the rule's higher-powered-outcome example.

## Source strategy
- NCBI exposes PMC Open Access full text through the official BioC REST API; it accepts batches of PMCIDs and returns structured JSON/XML.
- Three spot-checked PMC pages resolved correctly, including PMCID 11606077, 11443760, and 12257706.
- Use BioC full text as the batch evidence corpus and retain the canonical PMC article URL per paper for audit notes.
- The TSV contains 46 distinct PMCIDs, each paired with one C1 row and one C3 row.
- A live two-PMCID batch test returned two complete BioC collections; batch retrieval is suitable and reduces request count.
- Batch retrieval produced structured full text for 45/46 papers. PMC13448146 (PMID 42563458) is the sole missing BioC record and requires a canonical-page fallback.
- PMC13448146 is live at its canonical PMC URL (HTTP 200; ~301 KB HTML) despite not yet being exposed by BioC or Europe PMC fullTextXML. The canonical page is the fallback source.

---

# Findings — Four-part methodological upgrade

## Project-specific constraints
- The current draft uses `coloc.abf` across full cis windows and already documents prior/window sensitivity. The unresolved risk is the one-causal-variant assumption.
- The available CD4 dynamic-eQTL cohorts contain only about 85–100 donors and do not provide in-sample LD. The existing 525-person proxy-LD SuSiE run produced 19 MC1R credible sets versus 3 in FinnGen in-sample fine-mapping, so an external-LD SuSiE count cannot be treated as truth without a reliability gate.
- FinnGen R8–R12 already form a fixed-endpoint, nested power series with 2,000 R12 down-sampling replicates. Extending this to regional colocalisation is feasible, but the nested releases are correlated and must not be analysed as independent replications.
- S34 already implemented process blinding against the 2019 Open Targets gold-standard snapshot, but failed its preregistered coverage floor: 2 trait-matched and 16 cross-trait evaluable loci. This analysis must not be reopened by lowering the floor; a new benchmark must expand the prospectively fixed outcome/locus universe.

## Verified method anchors
- Wallace 2021 showed that `coloc.susie` compares each detected signal pair and is generally more accurate than conditioning-based alternatives when multiple causal variants exist.
- Official `coloc` documentation requires dense regional coverage and the same SNP set in both traits; `runsusie` requires a signed LD matrix and sample size and exposes convergence information.
- Current official `susieR` diagnostics explicitly warn that external/finite LD mismatch can create spurious extra credible sets; allele checks, kriging, finite-reference correction, mismatch modelling, convergence and reliability diagnostics are required before interpretation.
- FinnGen publicly provides release-specific summary statistics, GRCh38 variant definitions, and release-specific fine-mapping/LD resources; phenotype definitions still require code-level identity checks before a release enters a trajectory.
- The Open Targets 2019 gold-standard snapshot is an external, versioned locus-to-gene resource, but its evidence is mixed. Molecular-QTL/functional-observational assignments are circular for this benchmark and L2G predictions themselves are not truth labels.

## Design consequences
- Multiple-signal analysis should be a gated triangulation: baseline `coloc.abf`, `coloc.susie` only where LD/fit diagnostics pass, a conditioning or masking sensitivity, and an explicit `unresolved` state when they do not.
- The primary multiple-signal universe should be fixed independently of MR significance (all analysable exposure-region records or a preregistered representative subset), otherwise candidate-only analysis introduces selection bias.
- The FinnGen trajectory should use fixed regions, fixed variant intersections, fixed priors, and fixed exposure data across R8–R12. Continuous regional metrics and categorical transitions should be reported; monotonicity of individual loci should not be assumed.
- A stronger attribution benchmark should keep the Open Targets gene field masked while using only coordinates/traits to prospectively select a larger panel of FinnGen outcomes. It should estimate both end-to-end recovery and conditional attribution, with coverage and structural non-instrumentability kept separate.
- The reusable package should be built before the confirmatory runs, with schemas, synthetic truth fixtures, real positive/negative controls, claim-level verdicts, an HTML report, lockfile/container and CI. The package is infrastructure, not independent evidence, until it reproduces the frozen analyses and passes external/synthetic benchmarks.
