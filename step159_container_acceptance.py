#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 159 -- execute S54: the container acceptance run.

Executes manuscript/PREREG_container_canonical.md (S54), committed 2026-09-07
before the image existed. S54 fixes the must-run list, its order, the two tiers
of acceptance criterion, the all-or-nothing clause and the stopping rules. This
script runs the list and reports; it decides nothing that S54 has not already
decided.

RUN IT INSIDE THE CONTAINER, FROM THE REPOSITORY ROOT

  docker build -f container/Dockerfile -t cqtna-audit:0.3.0 .
  docker run --rm -v "$PWD:/repo" -w /repo cqtna-audit:0.3.0 \
      python3 step159_container_acceptance.py /repo

ORDER IS A STOPPING RULE, NOT A PREFERENCE
S54 section 8: if the four audits or the package suite fail, the environment is
not ready, and we do not go on to explain individual numbers. That gate is in
the code. Groups 3 to 7 run only if groups 1 and 2 pass.

WHAT IT DOES NOT DO
It does not adopt anything. S54 section 4 forbids taking the container's
numbers where they agree and keeping the authoring machine's where they do not;
the decision between "canonical" and "fall back to (a)" is made by a person,
against this report, and recorded in S54 section 9.

Outputs: 159a_container_acceptance.tsv, 159b_console.log
         plus 159c_snapshot/ -- copies of the result tables this run produced,
         so a later diff against the authoring machine's tables is possible.
"""
import io
import os
import platform
import shutil
import subprocess
import sys
import time

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(MR, "159c_snapshot")
LOG = []

# S54 section 3, in the order S54 fixes. group 1 and 2 are the gate.
GROUPS = [
    (1, "audits -- the whole reconciliation surface at once", [
        ("step127 numbers", ["python3", "step127_audit_manuscript_numbers.py", MR]),
        ("step153 pre-registrations", ["python3", "step153_prereg_inventory.py", MR]),
        ("step154 cross-references", ["python3", "step154_crossref_audit.py", MR]),
        ("step155 reference order", ["python3", "step155_renumber_references.py", MR]),
    ]),
    (2, "the package's own assertions", [
        ("cqtna testthat", ["Rscript", "-e",
                            'testthat::test_local("cqtna_r")']),
    ]),
    (3, "the central claim", [
        ("step140 identity", ["Rscript", "step140_estimator_identity.R", MR]),
        ("step141 decomposition", ["Rscript", "step141_enrichment_decomposition.R", MR]),
    ]),
    (4, "the most sensitive floating-point probe", [
        ("step101 effect-size matching", ["python3", "step101_effect_size_matching.py", MR]),
    ]),
    (5, "the one S54 expects to move", [
        ("step147 fine-mapping", ["Rscript", "step147_finemap_decomposition.R", MR]),
    ]),
    (6, "text and contingency -- should not move at all", [
        ("step151 corpus scan", ["python3", "step151_instrument_count_in_literature.py", MR]),
        ("step156 C6 score", ["python3", "step156_estimator_provenance_coding.py", MR, "score"]),
        ("step158 visibility", ["python3", "step158_visibility_overlap.py", MR, MR]),
    ]),
]

# Tables worth snapshotting for a later diff. S54 section 2 tier two lives here.
SNAPSHOT = [
    "140a_identity_check.tsv", "140d_threshold.tsv",
    "141a_enrichment_decomposition.tsv", "141b_reviewer_conditioning.tsv",
    "101a_zdist.tsv", "101b_matched.tsv", "101c_stratified.tsv",
    "101d_zonly_model.tsv", "147a_finemap_decomposition.tsv",
    "151a_estimator_usage.tsv", "156e_C6_result.tsv",
    "158a_visibility_2x2.tsv", "158b_visibility_by_decile.tsv",
    "123d_fixed_anchor_full_grid.tsv", "130c_multilist_verdict.tsv",
]


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def run(label, cmd):
    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=MR, capture_output=True, text=True,
                           timeout=7200)
        rc, out = p.returncode, (p.stdout or "") + (p.stderr or "")
    except FileNotFoundError as e:
        rc, out = 127, "not found: %s" % e
    except subprocess.TimeoutExpired:
        rc, out = 124, "timed out after 7200 s"
    say("  %-32s %-6s %6.1f s" % (label, "ok" if rc == 0 else "FAIL rc=%d" % rc,
                                  time.time() - t0))
    if rc != 0:
        for line in [x for x in out.splitlines() if x.strip()][-12:]:
            say("      | " + line[:110])
    return rc, out


def main():
    say("container acceptance run (S54), %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    say("  platform : %s" % platform.platform())
    say("  python   : %s" % platform.python_version())
    say("  cwd      : %s" % MR)
    say()
    say("S54 section 8: groups 1 and 2 are a gate. If they fail the environment")
    say("is not ready and this script stops rather than explaining single numbers.")

    rows, gate_failed = [], False
    for num, title, items in GROUPS:
        say()
        say("=" * 74)
        say("group %d -- %s" % (num, title))
        say("=" * 74)
        if gate_failed:
            say("  skipped: the gate failed (S54 section 8).")
            for label, _ in items:
                rows.append(dict(group=num, step=label, status="skipped", rc=""))
            continue
        for label, cmd in items:
            rc, _ = run(label, cmd)
            rows.append(dict(group=num, step=label,
                             status="ok" if rc == 0 else "FAIL", rc=rc))
            if rc != 0 and num <= 2:
                gate_failed = True
        if gate_failed:
            say()
            say("  GATE FAILED in group %d. Stopping." % num)
            say("  Fix the image first. S54 forbids running on and explaining")
            say("  individual numbers from an environment that is not ready.")

    # snapshot whatever exists, so a diff against the authoring machine is possible
    if not os.path.isdir(SNAP):
        os.makedirs(SNAP)
    kept = 0
    for f in SNAPSHOT:
        src = os.path.join(MR, f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(SNAP, f))
            kept += 1
    say()
    say("snapshotted %d of %d result tables into 159c_snapshot/" % (kept, len(SNAPSHOT)))
    say("Diff them against the authoring machine's copies to fill S54 section 9.2")
    say("and 9.3. This script does not diff them itself: the authoring copies are")
    say("what git already holds, so `git status` and `git diff` are the diff.")

    with io.open(os.path.join(MR, "159a_container_acceptance.tsv"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write("group\tstep\tstatus\trc\n")
        for r in rows:
            f.write("%s\t%s\t%s\t%s\n"
                    % (r["group"], r["step"], r["status"], r["rc"]))

    bad = [r for r in rows if r["status"] == "FAIL"]
    say()
    say("=" * 74)
    if gate_failed:
        say("VERDICT: environment not ready. Nothing about individual numbers")
        say("is interpretable from this run (S54 section 8).")
        return 2
    if bad:
        say("VERDICT: gate passed, %d later step(s) failed. Record them in S54"
            % len(bad))
        say("section 9.4 as not rerunnable, with the reason, before deciding.")
        return 1
    say("VERDICT: every step in the must-run list completed.")
    say()
    say("This is NOT yet a decision. S54 section 4 is all-or-nothing: read the")
    say("audits' output and the table diffs, list every number that moved in 9.2")
    say("and every verdict that flipped in 9.3, then choose canonical or fall")
    say("back to (a) -- and record which, with reasons.")
    return 0


if __name__ == "__main__":
    rc = main()
    io.open(os.path.join(MR, "159b_console.log"), "w",
            encoding="utf-8", newline="\n").write("\n".join(LOG) + "\n")
    raise SystemExit(rc)
