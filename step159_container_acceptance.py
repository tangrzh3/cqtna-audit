#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 159 -- execute S54: the container acceptance run.

Executes manuscript/PREREG_container_canonical.md (S54), committed 2026-09-07
before the image existed, as amended in S54 section 10.1 (before the image
existed, and before any result was observed).

WHY THE ORDER IS WHAT IT IS -- THIS IS THE POINT OF THE SCRIPT
The first version ran step127 first, against the tables mounted in from the
authoring machine. That proves the repository arrived intact and nothing more:
it cannot show that the container reproduces the manuscript, because the
container had not yet produced anything. Third-party review caught it. The
order now is

  phase 0  snapshot the incoming tables as an immutable baseline
  phase 1  gate: package suite + the four audits, as a MOUNT CHECK only
  phase 2  rerun every analysis that can run here
  phase 3  ACCEPTANCE: rerun the four audits, now against container tables
  phase 4  diff container tables against the phase-0 baseline, cell by cell

Phase 3 is the tier-one criterion in S54 section 2. Phase 1 is not.

A failed step cannot leave tracked files worse than it found them, either:
`run()` diffs tracked files immediately before and after each step and, on
failure, reverts exactly what that step touched (needs `git` in the image;
the acceptance Dockerfile installs it after the expensive R layer so this
does not invalidate that cache). Found necessary the first time this ran:
step158 failed for a legitimate, documented reason -- its external inputs are
not part of the deposit -- and in doing so overwrote a prior run's console
log with a one-line failure message. Nothing was lost from git history, but a
tool meant for strangers to rerun should not depend on someone noticing that
by hand every time.

  docker build -f container/Dockerfile -t cqtna-audit:0.3.0 .
  docker run --rm -v "$PWD:/repo" -w /repo cqtna-audit:0.3.0 \
      python3 step159_container_acceptance.py /repo

WHAT IT DOES NOT DO
It does not adopt anything. S54 section 4 is all-or-nothing and the choice
between canonical and falling back to (a) is a person's, made against this
report and recorded in S54 section 9.

Outputs: 159a_container_acceptance.tsv, 159b_console.log,
         159d_baseline/ (immutable), 159e_table_diff.tsv
"""
import filecmp
import glob
import io
import os
import platform
import shutil
import subprocess
import sys
import time

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(MR, "159d_baseline")
LOG = []

AUDITS = [
    ("step127 numbers", ["python3", "step127_audit_manuscript_numbers.py", MR]),
    ("step153 pre-registrations", ["python3", "step153_prereg_inventory.py", MR]),
    ("step154 cross-references", ["python3", "step154_crossref_audit.py", MR]),
    ("step155 reference order", ["python3", "step155_renumber_references.py", MR]),
]

# S54 section 3, groups 3 to 6, in S54's order.
ANALYSES = [
    ("step140 identity", ["Rscript", "step140_estimator_identity.R", MR]),
    ("step141 decomposition", ["Rscript", "step141_enrichment_decomposition.R", MR]),
    ("step101 effect-size matching", ["python3", "step101_effect_size_matching.py", MR]),
    ("step147 fine-mapping", ["Rscript", "step147_finemap_decomposition.R", MR]),
    ("step151 corpus scan", ["python3", "step151_instrument_count_in_literature.py", MR]),
    ("step156 C6 score", ["python3", "step156_estimator_provenance_coding.py", MR, "score"]),
    ("step158 visibility", ["python3", "step158_visibility_overlap.py", MR, MR]),
]

# S54 section 3 group 7: everything else with its inputs present. Attempted, and
# a failure here is recorded as "not rerunnable" rather than as a discrepancy.
GROUP7_GLOB = "step1[2-5][0-9]_*.py"
GROUP7_SKIP = {  # already run above, or not analyses
    "step151_instrument_count_in_literature.py",
    "step152_set_identity.py", "step153_prereg_inventory.py",
    "step154_crossref_audit.py", "step155_renumber_references.py",
    "step156_estimator_provenance_coding.py", "step157_perturbseq_tpi1_lookup.py",
    "step158_visibility_overlap.py", "step159_container_acceptance.py",
    "step127_audit_manuscript_numbers.py",
    # Queries a LIVE external database (GWAS Catalog) by design -- its own
    # docstring says "the Catalog grows, so some drift is expected". A rerun
    # two weeks later legitimately returns different counts on the SAME
    # machine; that is not a container finding and running it here only
    # produces noise S54 is not asking about. Confirmed by inspection after
    # the first acceptance run showed 68 changed lines here with nothing
    # else to explain them but the calendar.
    "step143_list_provenance.py",
}

# Tier-two tables (S54 section 2) plus the ones the manuscript quotes from.
WATCH = [
    "140a_identity_check.tsv", "140d_threshold.tsv",
    "141a_enrichment_decomposition.tsv", "141b_reviewer_conditioning.tsv",
    "101a_zdist.tsv", "101b_matched.tsv", "101c_stratified.tsv",
    "101d_zonly_model.tsv", "147a_finemap_decomposition.tsv",
    "151a_estimator_usage.tsv", "156c_C6_agreement.tsv", "156e_C6_result.tsv",
    "158a_visibility_2x2.tsv", "158b_visibility_by_decile.tsv",
    "123d_fixed_anchor_full_grid.tsv", "130c_multilist_verdict.tsv",
    "126a_offgrid_attribution.tsv", "85e_matched_background_fixed_anchor.tsv",
]


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def _git(args):
    try:
        p = subprocess.run(["git"] + args, cwd=MR, capture_output=True,
                           text=True, timeout=60)
        return p.returncode, (p.stdout or "")
    except Exception:
        return 1, ""


def run(label, cmd, timeout=7200):
    # A step that cannot succeed here should not be able to leave tracked
    # files worse off than it found them. Discovered the hard way: step158
    # failed for a documented, legitimate reason (its external inputs are not
    # part of the deposit) and, in doing so, overwrote 158c_console.log --
    # 33 lines of a prior successful run's evidence -- with one line saying
    # it could not run this time. Nothing was lost from git history, but a
    # tool meant to be rerun by strangers should not depend on a human
    # noticing that in every future run. `git diff --name-only` before and
    # after brackets exactly what THIS step touched; on failure, revert only
    # that set, never anything a different step is responsible for.
    rc_git, before = _git(["diff", "--name-only"])
    before = set(before.split()) if rc_git == 0 else None

    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=MR, capture_output=True, text=True,
                           timeout=timeout)
        rc, out = p.returncode, (p.stdout or "") + (p.stderr or "")
    except FileNotFoundError as e:
        rc, out = 127, "not found: %s" % e
    except subprocess.TimeoutExpired:
        rc, out = 124, "timed out"
    say("  %-34s %-12s %6.1f s"
        % (label, "ok" if rc == 0 else "FAIL rc=%d" % rc, time.time() - t0))
    if rc != 0:
        for line in [x for x in out.splitlines() if x.strip()][-8:]:
            say("      | " + line[:108])
        if before is not None:
            rc_git, after = _git(["diff", "--name-only"])
            touched = (set(after.split()) - before) if rc_git == 0 else set()
            if touched:
                _git(["checkout", "--"] + sorted(touched))
                say("      reverted (failed step touched tracked files): %s"
                    % ", ".join(sorted(touched)))
    return rc, out


def cellwise_diff(a, b, limit=40):
    """Report changed cells, not just changed files: a moved digit is what S54
    section 5 has to enumerate."""
    try:
        ra = io.open(a, encoding="utf-8", errors="replace").read().split("\n")
        rb = io.open(b, encoding="utf-8", errors="replace").read().split("\n")
    except Exception as e:
        return [("<unreadable>", "", str(e))]
    out = []
    if ra and rb and ra[0] != rb[0]:
        out.append(("<header>", ra[0][:60], rb[0][:60]))
    for i in range(max(len(ra), len(rb))):
        la = ra[i] if i < len(ra) else "<missing row>"
        lb = rb[i] if i < len(rb) else "<missing row>"
        if la == lb:
            continue
        ca, cb = la.split("\t"), lb.split("\t")
        for j in range(max(len(ca), len(cb))):
            x = ca[j] if j < len(ca) else ""
            y = cb[j] if j < len(cb) else ""
            if x != y:
                out.append(("row %d col %d" % (i, j + 1), x[:40], y[:40]))
                if len(out) >= limit:
                    return out
    return out


def main():
    say("container acceptance run (S54, as amended in section 10.1)")
    say("  when     : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    say("  platform : %s" % platform.platform())
    say("  python   : %s" % platform.python_version())
    rows = []

    # ------------------------------------------------------------- phase 0
    say()
    say("=" * 74)
    say("phase 0 -- immutable baseline of the tables mounted in")
    say("=" * 74)
    if os.path.isdir(BASE):
        say("  159d_baseline/ already exists; refusing to overwrite it.")
        say("  A baseline taken after a rerun is not a baseline. Remove it")
        say("  deliberately if you really are starting over.")
        return 2
    os.makedirs(BASE)
    kept = 0
    for f in WATCH:
        src = os.path.join(MR, f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(BASE, f))
            kept += 1
    say("  copied %d of %d watched tables into 159d_baseline/" % (kept, len(WATCH)))

    # ------------------------------------------------------------- phase 1
    say()
    say("=" * 74)
    say("phase 1 -- GATE. Mount check only, NOT the acceptance test")
    say("=" * 74)
    say("  These audits run against the tables that arrived with the repository.")
    say("  Passing here says the mount is intact. It says nothing yet about")
    say("  whether the container reproduces them (S54 section 10.1).")
    gate_bad = False
    rc, _ = run("cqtna testthat", ["Rscript", "-e",
                                   'testthat::test_local("cqtna_r")'])
    rows.append(dict(phase=1, step="cqtna testthat",
                     status="ok" if rc == 0 else "FAIL", rc=rc))
    gate_bad |= rc != 0
    for label, cmd in AUDITS:
        rc, _ = run(label + " (pre)", cmd)
        rows.append(dict(phase=1, step=label + " (pre)",
                         status="ok" if rc == 0 else "FAIL", rc=rc))
        gate_bad |= rc != 0
    if gate_bad:
        say()
        say("  GATE FAILED. The environment or the mount is not ready.")
        say("  S54 section 8 forbids going on to explain individual numbers.")
        write_rows(rows)
        return 2

    # ------------------------------------------------------------- phase 2
    say()
    say("=" * 74)
    say("phase 2 -- rerun the analyses; these produce the container's tables")
    say("=" * 74)
    for label, cmd in ANALYSES:
        rc, _ = run(label, cmd)
        rows.append(dict(phase=2, step=label,
                         status="ok" if rc == 0 else "not rerunnable", rc=rc))

    say()
    say("  group 7 -- every other numbered step whose inputs are present")
    done = set()
    for path in sorted(glob.glob(os.path.join(MR, GROUP7_GLOB))):
        name = os.path.basename(path)
        if name in GROUP7_SKIP or name in done:
            continue
        done.add(name)
        rc, _ = run(name, ["python3", name, MR], timeout=3600)
        rows.append(dict(phase=2, step=name,
                         status="ok" if rc == 0 else "not rerunnable", rc=rc))

    # ------------------------------------------------------------- phase 3
    say()
    say("=" * 74)
    say("phase 3 -- ACCEPTANCE. The same audits, now against container tables")
    say("=" * 74)
    acc_bad = 0
    for label, cmd in AUDITS:
        rc, _ = run(label + " (post)", cmd)
        rows.append(dict(phase=3, step=label + " (post)",
                         status="ok" if rc == 0 else "FAIL", rc=rc))
        acc_bad += rc != 0

    # ------------------------------------------------------------- phase 4
    say()
    say("=" * 74)
    say("phase 4 -- container tables against the phase-0 baseline")
    say("=" * 74)
    moved, diffs = [], []
    for f in WATCH:
        a, b = os.path.join(BASE, f), os.path.join(MR, f)
        if not os.path.exists(a) or not os.path.exists(b):
            continue
        if filecmp.cmp(a, b, shallow=False):
            continue
        moved.append(f)
        for where, was, now in cellwise_diff(a, b):
            diffs.append(dict(table=f, where=where, baseline=was, container=now))
    if not moved:
        say("  every watched table is byte-identical to the baseline.")
    else:
        say("  %d table(s) differ:" % len(moved))
        for f in moved:
            n = sum(1 for d in diffs if d["table"] == f)
            say("    %-44s %d changed cell(s)" % (f, n))
    with io.open(os.path.join(MR, "159e_table_diff.tsv"), "w",
                 encoding="utf-8", newline="\n") as fh:
        fh.write("table\twhere\tbaseline\tcontainer\n")
        for d in diffs:
            fh.write("%s\t%s\t%s\t%s\n"
                     % (d["table"], d["where"], d["baseline"], d["container"]))
    say("  wrote 159e_table_diff.tsv")

    write_rows(rows)
    say()
    say("=" * 74)
    if acc_bad:
        say("VERDICT: %d audit(s) FAILED in phase 3." % acc_bad)
        say("Those failures are the numbers that moved. Enumerate every one in")
        say("S54 section 9.2 and every flipped verdict in 9.3, then apply the")
        say("migration rule in S54 section 5. Do not adopt selectively:")
        say("S54 section 4 is all-or-nothing.")
        return 1
    say("VERDICT: phase 3 audits all pass against container-produced tables,")
    say("and %d watched table(s) differ from the baseline." % len(moved))
    say()
    say("This is still NOT a decision. Read 159e_table_diff.tsv, fill S54")
    say("section 9, and choose canonical or fall back to (a) with reasons.")
    return 0


def write_rows(rows):
    with io.open(os.path.join(MR, "159a_container_acceptance.tsv"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write("phase\tstep\tstatus\trc\n")
        for r in rows:
            f.write("%s\t%s\t%s\t%s\n"
                    % (r["phase"], r["step"], r["status"], r["rc"]))


if __name__ == "__main__":
    rc = main()
    io.open(os.path.join(MR, "159b_console.log"), "w",
            encoding="utf-8", newline="\n").write("\n".join(LOG) + "\n")
    raise SystemExit(rc)
