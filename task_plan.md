# Task plan — Compare MANUSCRIPT_GB with MANUSCRIPT_NC

## Goal
Read `manuscript/MANUSCRIPT_GB.md`, compare it directly with the already reviewed `MANUSCRIPT_NC.md`, and decide which is the stronger manuscript and why.

## Phases
- [completed] 1. Map structural and claim differences between GB and NC.
- [completed] 2. Read GB in full and audit its evidence hierarchy.
- [completed] 3. Compare both versions across originality, significance, technical soundness, broad readership and readability.
- [completed] 4. Deliver a clear choice, trade-offs and recommended base version.

## Boundaries
- Compare the manuscripts as written, not hypothetical journal acceptance.
- Reuse only criticisms rechecked against GB; do not assume NC problems persist.
- No manuscript edits are authorized in this turn.

## Errors encountered
| Error | Attempt | Resolution |
|---|---:|---|
| PowerShell rejected a pipeline directly after a `foreach` block in the comparison count command | 1 | Assign the loop output to a variable, then format it; do not repeat the invalid syntax. |
| Verification used an exact HCC phrase that differs across line wrapping/wording in NC | 1 | Verify the distinctive semantic fragment (`behave alike` + `power is equalised`) separately in each file. |

---

# Task plan — Four-part methodological upgrade

## Goal
Design an executable, preregisterable workflow for: (1) conditional/multiple-signal colocalisation, (2) nested FinnGen-release colocalisation trajectories, (3) an independent blinded gene-attribution benchmark, and (4) a reusable implementation of the eight diagnostic checks.

## Phases
- [completed] 1. Reconstruct the manuscript-specific estimands, candidate universes, and current evidence limitations.
- [completed] 2. Verify primary methods and official implementations for multiple-signal colocalisation and fine-mapping.
- [completed] 3. Specify each analysis with frozen inputs, controls, QC gates, estimands, failure criteria, and outputs.
- [completed] 4. Integrate the four analyses into a staged package/repository architecture and manuscript decision tree.
- [completed] 5. Deliver a prioritized implementation sequence and minimum/ideal versions.

## Boundaries
- Design only; do not execute new analyses or edit the manuscript in this turn.
- Treat nested FinnGen releases as correlated observations, not independent replication cohorts.
- Do not treat an external gene-prioritization score as ground truth if it reuses eQTL/colocalisation evidence being benchmarked.
- Real-data positive controls must be supplemented by simulation because causal truth is not fully known at complex loci.
- Any analysis failing its preregistered QC/positive-control gate is labelled uninformative, not negative.

## Errors encountered
| Error | Attempt | Resolution |
|---|---:|---|
| A stale path was used for the literature workflow (`multi-source-search.md`) | 1 | Loaded the manifest-declared `wf1-multi-source-search.md`; no task analysis was affected. |

---

# Task plan — Blind coding of 105a literature set

## Goal
Research the full text of every paper represented in `105a_blind_coding.tsv`, apply only the definitions in `105_CODING_RULES.txt`, and fill `coder2` (with concise audit notes when useful) without changing the source columns.

## Phases
- [completed] 1. Load required literature/spreadsheet instructions and inventory the 46-paper, 92-row dataset.
- [completed] 2. Retrieve authoritative PubMed Central/full-text sources and build a paper-level evidence cache.
- [completed] 3. Adjudicate C1_known_locus and C3_power_stability independently for all papers.
- [completed] 4. Write only coder2/coder2_note to a preserved-copy TSV through the required spreadsheet tooling.
- [completed] 5. Verify row preservation, value domain, coverage, and representative evidence against full text.

## Boundaries
- The coding rules control; article wording, abstracts, reference titles, and generic power statements do not broaden them.
- Use the paper's own outcome and candidate list as the unit of judgment.
- Preserve row order and all source fields exactly.
- Leave blank only when full-text evidence truly cannot resolve the criterion.

## Errors encountered
| Error | Attempt | Resolution |
|---|---:|---|
| Router example referenced a nonexistent `multi-source-search.md` path | 1 | Use the manifest-declared `wf1-multi-source-search.md` path. |
| PMC BioC Open Access returned 45/46 articles; PMC13448146 absent | 1 | Retrieve PMC13448146 from its canonical PMC page or journal/PubMed fallback and record the source limitation. |
| Web search/open did not index or safely open the newly assigned PMC13448146 URL | 1 | Use a direct read-only HTTP request to the canonical PMC page, which returned HTTP 200 and complete HTML. |
