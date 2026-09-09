# Deposit gaps: known, recorded, and not closable today

Things the deposit does not contain that it arguably should. Each is reported
by an audit on every run rather than fixed silently, and each says why it
cannot simply be repaired.

---

## 1. `150a_environment.lock` omits two packages the Methods names

`chromVAR 1.26.0` and `org.Hs.eg.db 3.19.1` appear in the Software section and
in `container/install_r_packages.R`'s `need` list, but **not** among the lock's
226 entries. Confirmed by key against every entry, not by a failed lookup.

In `cqtna-audit:0.3.0` both are present **at exactly the stated versions** —
installed by the fallback that handles anything `renv::restore()` did not
restore. ⚠ **That path is not version-pinned**, so they are correct by luck
rather than by record; a rebuild on another day could get different versions
and nothing would say so.

**Why it is not simply fixed**: regenerating the lock means running `step150`,
which rebuilds it from `installed.packages()` on the **current** machine. That
machine has drifted — numpy 2.5.2 against the recorded 2.0.0 (S54 §9.8) — so
regenerating would overwrite the authoring-machine record with a later one and
destroy the baseline every S54 comparison is made against. **The record being
incomplete is better than the record being wrong.**

**Reported by**: `step127` §1n, as `GAP`, every run.

⚠ **Consequence for the closure check**: `step159`'s gate compares installed
packages against the lock, so a package the lock omits is outside its frame.
226/226 was green while these two went unexamined. The closure figure means
"everything the lock records is present at the recorded version", not
"everything the paper uses is pinned".

## 2. PLINK 2.0.0-a.7.2 is no longer obtainable

`s3://plink2-assets/alpha7/` now holds only 2026 builds;
`plink2_linux_x86_64_20241206.zip`, the version the Methods names, has been
deleted upstream (established by listing the bucket, not by reading a 404).

**Not substituted.** Putting a current build behind a version number the paper
states would make the image quietly wrong rather than openly incomplete.
`/opt/binaries_provenance.txt` in the image records SMR's checksum and PLINK's
absence.

**No analysis on S54's must-run list invokes either binary**, so the acceptance
is unaffected; a reader trying to rerun everything is.

## 3. Four analyses cannot run inside the container

`step158` (S53 visibility overlap) needs externally downloaded bulk RNA-seq;
`step120`, `step121`, `step122` point at authoring-machine paths and at
`landi2020_known_loci_grch38.csv`, which is not distributed. These are data
distribution decisions taken before the container existed, not container
defects, and `step159` records them as `not rerunnable` rather than as
failures.

## 4. Audit coverage of the manuscript's arithmetic is partial

`step160` measures it directly. It is not complete and the uncovered values are
listed in `160a_audit_coverage.tsv`, so a reviewer can see exactly which
numbers no audit checks rather than inferring that a green run means every
number was verified. Two errors have already been found inside that uncovered
set (S54 §9.8's wrong-table 13.53, and the Steiger lower bound), which is the
argument for keeping the list visible rather than quoting a coverage headline.
