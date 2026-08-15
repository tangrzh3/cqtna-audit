# Progress — MANUSCRIPT_NC recent MR benchmark

- 2026-08-13: Initialized manuscript audit and recent-literature benchmark.
- 2026-08-13: Read full 297-line NC draft; completed initial claim/length/limitation audit.
- 2026-08-13: Broad PubMed query was too non-specific (6,571 hits); switching to title-level and journal-filtered searches.
- 2026-08-13: Extracted manuscript claim-calibration issues and began primary-source verification of recent high-impact MR comparators.
- 2026-08-13: Verified three recent high-impact original designs (Nature Medicine, EHJ, Circulation) and one current Nature Genetics methodological benchmark.
- 2026-08-13: Mapped comparator standards to manuscript-specific priorities and verified FinnGen DF13 availability plus cancer-endpoint definition changes.
- 2026-08-13: Completed prioritized three-reviewer assessment and cross-review synthesis.
- 2026-08-13: Began a fresh severe review of the current manuscript version; reset the active plan while preserving prior benchmark findings.
- 2026-08-13: Read title, abstract, introduction and Results through the crossed disease/resource analysis; logged central claim contradictions and unit-of-analysis risks.
- 2026-08-13: Read the remainder of Results and Discussion plus relevant Methods sections; logged multiplicity, comparator-confounding, colocalisation, gene-attribution, literature-audit and reproducibility limitations.
- 2026-08-13: Located the exact primary FDR implementation; began independent recalculation because duplicated gene×profile records share outcome P values. An R PATH check failed and was logged; switching implementation for verification.
- 2026-08-13: Independently verified the primary record-level BH and quantified unit sensitivity: 10 reported genes versus 11 unique-SNP-BH genes and 13 one-minimum-P-per-gene genes.
- 2026-08-13: Visually inspected the main locus-attribution and power-stability figures; logged graphic overclaim and manuscript/figure-plan numbering mismatch.
- 2026-08-13: Final verification passed: manuscript assertions, companion-Methods placeholder, and independent BH sensitivity counts were confirmed fresh (`VERIFICATION_OK`).
- 2026-08-13: Began direct comparison of GB and NC versions.
- 2026-08-13: Compared structure and read GB through the crossed disease/resource analysis; logged completeness gains and persistent/stronger overclaims.
- 2026-08-13: Read GB in full and compared claim flags, limitations, evidence ladder and narrative branches against NC. Preliminary decision: GB is the stronger base, NC the compression template.
- 2026-08-13: Fresh comparison verification passed (`COMPARISON_VERIFIED`); decision and supporting textual differences confirmed.

---

# Progress — 105a blind coding

- 2026-08-14: Confirmed authorization to research online and fill the TSV.
- 2026-08-14: Loaded the coding rules and inventoried 92 rows / 46 papers; no existing coder2 values.
- 2026-08-14: Started literature-search, spreadsheet-editing, and persistent-audit workflows.
- 2026-08-14: Logged and corrected a stale literature-workflow example path using the manifest-declared filename.
- 2026-08-14: Completed required source-tier, deduplication, spreadsheet API, formatting, and scientific-research instruction reads.
- 2026-08-14: Selected Tier-1 PubMed/PMC full text as the controlling evidence source; secondary sources will be used only if PMC retrieval is incomplete.
- 2026-08-14: Verified the official NCBI PMC BioC full-text API and canonical PMC pages through live web search.
- 2026-08-14: Loaded the bundled spreadsheet runtime required for final TSV authoring.
- 2026-08-14: Live-tested BioC JSON retrieval on PMC11606077; structured full text, DOI, title, section metadata, and Unicode text are available.
- 2026-08-14: Confirmed 46 unique PMCIDs and successful multi-ID BioC batch behavior.
- 2026-08-14: Retrieved and locally indexed 45/46 official PMC full texts; generated criterion-specific candidate passages for manual adjudication.
- 2026-08-14: Logged PMC13448146 as the only BioC Open Access gap and started fallback retrieval.
- 2026-08-14: Confirmed canonical PMC13448146 HTML is available and selected it as the fallback full-text source after search indexing/BioC gaps.
- 2026-08-14: Screened all 30 nonempty evidence rows against the strict definitions; identified three clear C3 positives and no clear C1 positive.
- 2026-08-14: Completed 46/46 full-text retrieval, including canonical HTML fallback for PMC13448146.
- 2026-08-14: Resolved key boundary cases by full-text context while preserving the evidence-only coding rule.
- 2026-08-14: Final adjudication set contains 89 zeros, 3 ones, and 0 blanks; positive rows are 3, 18, and 47.
- 2026-08-14: Authored `105a_blind_coding_coded.tsv` through the required spreadsheet runtime and exported a verification XLSX/preview.
- 2026-08-14: Visual QA found fixed-row-height clipping in the support XLSX; changed only row sizing to automatic and rerendered successfully.
- 2026-08-14: Fresh final verification passed (`VERIFICATION_OK`): 92 rows, 0 source-field mismatches, 89 zeros, 3 ones, 0 blanks, 0 blank notes; positive rows 3/18/47.

---

# Progress — Four-part methodological upgrade

- 2026-08-15: Started a preregistration-oriented design for multi-signal coloc, FinnGen release trajectories, blinded gene attribution, and an eight-diagnostic software package.
- 2026-08-15: Loaded the project’s existing review findings and the required multi-source literature-search workflow.
- 2026-08-15: Established three design constraints: nested releases are correlated; gene truth must be independent of eQTL evidence; failed positive-control/QC gates produce an uninformative result.
- 2026-08-15: Verified the current manuscript’s exact colocalisation, release-trajectory, S34 benchmark, and eight-diagnostic definitions.
- 2026-08-15: Verified primary/official sources for `coloc.susie`, coloc data requirements, SuSiE external-LD mismatch diagnostics, FinnGen release data, and Open Targets gold-standard provenance.
- 2026-08-15: Completed a staged design with frozen universes, control gates, estimands, failure labels, and software deliverables for all four upgrades.
