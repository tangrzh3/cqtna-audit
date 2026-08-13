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
