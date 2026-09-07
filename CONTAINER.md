# The container: what it reproduces, and what it does not

**Files**: `container/Dockerfile`, `container/requirements.txt`,
`container/install_r_packages.R`, plus `150a_environment.lock` in the repository
root.

```bash
docker build -f container/Dockerfile -t cqtna-audit:0.3.0 .
docker run --rm -v "$PWD:/repo" -w /repo cqtna-audit:0.3.0 \
    python3 step127_audit_manuscript_numbers.py /repo
```

**This image has not been built yet.** It was written from
`150b_session.txt` and the Methods, on a machine with no Docker. Three
things are most likely to need adjusting on the first build: the base tag
(`bioconductor/bioconductor_docker:RELEASE_3_19`, chosen because Bioc 3.19
is the release carrying the Bioc package versions in the lockfile), the
Ubuntu codename in the Posit snapshot URL in
`container/install_r_packages.R`, and the two external-binary URLs. Fix
them, then record what you built in section 4 rather than deleting this
paragraph.

---

## 1. Why a container and not just the lockfile

`150a_environment.lock` records 226 R packages at their exact versions. It does
not record the operating system, the BLAS, the compiler or the Python
interpreter, and those are not free parameters. This project has already
demonstrated, on its own data, that a result can move without any package
version changing: `susie_rss` returned 0, 1 and 0 credible sets on three
seedless runs of the same region, every one of them reporting
`converged = TRUE` with no warning, because `susie_get_cs` estimates purity from
a random sample of 100 variants and the 1,851-variant set sat on the threshold
(Supplementary S45). A package manifest cannot see that. So the container is the
stronger artefact — but only if its limits are stated rather than implied.

## 2. What it reproduces

Every dependency at the version the analyses ran under: R 4.4.1 with the package
closure in `150a_environment.lock`, Python 3.12 with the five pinned libraries
named in the Methods, and the two external binaries (SMR 1.3.1, PLINK
2.0.0-a.7.2). A reader who builds this image can run any numbered script and
have it resolve.

## 3. What it does NOT reproduce

**The reported numbers were produced on a different machine from the one this
image describes.** From `150b_session.txt`:

| | Authoring machine | This image |
|---|---|---|
| OS | Windows 11 x64 (build 26200) | Linux (Ubuntu, Bioconductor base) |
| Platform | `x86_64-w64-mingw32/x64` | `x86_64-pc-linux-gnu` |
| BLAS/LAPACK | R reference implementation (`Matrix products: default`) | OpenBLAS |
| Python | CPython 3.12.4 | CPython 3.12.x, patch level as shipped |
| Locale | `Chinese (Simplified)_China.utf8`, `LC_NUMERIC=C` | `C.UTF-8` |

Consequences worth naming rather than hoping about:

1. **Floating-point results may differ in the last digits**, and where a
   quantity sits on a threshold, a last-digit difference can change a verdict —
   which is exactly the failure mode S45 documents.
2. **Sort order can differ with the locale.** Any step that orders by a
   character key and then takes the first row is exposed to this.
3. **The external-binary URLs in the Dockerfile are not content-addressed.**
   Upstream can replace the bytes behind a stable URL. Record the checksums of
   what you actually built with.

## 4. The acceptance test, which has not been run

Building the image is not the deliverable; measuring the gap is. Run these
inside the container and record what agrees:

| Group | Check | Result |
|---|---|---|
| 1 (gate) | the four audits: `step127`, `step153`, `step154`, `step155` | not yet run |
| 2 (gate) | `testthat::test_local("cqtna_r")`, 271 assertions | not yet run |
| 3 | `step140` and `step141` — the identity and its decomposition | not yet run |
| 4 | `step101` — carries its own bitwise positive control | not yet run |
| 5 | `step147` — the one S54 expects to move | not yet run |
| 6 | `step151`, `step156 score`, `step158` — should not move at all | not yet run |

Groups 1 and 2 are a gate: if they fail, the environment is not ready and
nothing about individual numbers is interpretable (S54 section 8).

`step127` is the cheapest and the most informative: it re-derives every quoted
number from its table, so a clean run inside the container means the numbers in
the manuscript survive the environment change. `step147` is the one most likely
to move, for the reason in section 1; it already sweeps five seeds per region,
so compare its `n_cs_across_seeds` column rather than only its verdicts.

Fill this table in and keep it in the deposit. **A container with an unrun
acceptance test is a claim, not evidence** — and this project's own discipline
rules say a claim of that shape is what an audit is for.

## 5. The decision this leaves to the authors

Two coherent positions, and they should not be mixed:

**(a) The container is documentation.** The paper's numbers stay as they are,
authored on Windows; the container lets a reader re-run and see whether the
conclusions hold. Section 3 must then appear in the deposit, and section 4 must
be filled in, so the reader knows which numbers were checked across environments
and which were not.

**(b) The container is canonical.** Re-run the whole pipeline inside it and
report the container's numbers as the paper's numbers. This is the stronger
claim and the honest one for a reproducibility-focused paper, but it is not
free: any digit that moves has to be traced through
`manuscript/NUMBER_MIGRATION_fixed_anchor.md`, `step127` has to come back clean
against the new tables, and a verdict that flips has to be reported as having
flipped rather than quietly adopted.

**The author chose (b) on 2026-09-07.** What a moved number means was decided
in advance and is fixed in `manuscript/PREREG_container_canonical.md` (S54),
committed before this image was built: two tiers of acceptance criterion, an
all-or-nothing clause forbidding a mixture of container and authoring-machine
numbers, a migration rule for numbers that move, and a requirement that a
flipped verdict be reported as a flip rather than adopted silently. Run
`step159_container_acceptance.py` inside the image; it executes S54's must-run
list in S54's order and stops at S54's gate.

## 6. What the container still cannot fix

Seeds. `step147` already re-runs each region at five seeds and records
`n_cs_across_seeds` and `seed_stable` in `147a_finemap_decomposition.tsv`,
because `susie_rss` returned 0, 1 and 0 credible sets on three seedless runs of
KANSL1 with no warning at all. That sweep is the check that matters, and no
container makes it unnecessary: **reproducible is not the same as robust**.

What the container changes is whether the sweep still comes out the same. If
`seed_stable` flips for any region inside the image, the finding is
environment-dependent as well as seed-dependent, and that is a result about the
method rather than a build problem. Run `step147` first and compare
`n_cs_across_seeds` column for column.
