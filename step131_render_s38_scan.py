"""Step 131 -- render S38 section 7.1 from 126b rather than by hand.

The tolerance scan is eight cells by five tolerances with two numbers each, and
it has now been regenerated three times: once when the partition changed, once
when the estimand was frozen, and once when the sampler's order dependence was
fixed. Transcribing it by hand each time is how the primary column ended up
updated while the other four columns still carried pre-fix values.

Rewrites the table between the two markers in PREREG_permutation_estimand.md and
leaves the rest of the file alone.

  python step131_render_s38_scan.py [dir]
"""
import io
import os
import sys

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR")
      or os.path.dirname(os.path.abspath(__file__)))

TOLS = ["0.1", "0.25", "0.5", "1", "2"]
START = "| 格子 | 容差 0.10 | 0.25 | 0.50 | **1.00（主）** | 2.00 | 判读（§3）|"
VERDICT_OK = "**在本匹配规则下**富集未被密度零分布解释掉"
VERDICT_NO = "**未通过本密度零分布 → 降为描述性**"
VERDICT_MIX = "**跨 0.05 分界 → 不得引用任一 P**"


def fmt(p):
    return "P≤1e-4" if p <= 1.0001e-4 else f"P={p:.4f}"


def main():
    import csv
    rows = list(csv.DictReader(
        io.open(os.path.join(MR, "126b_permutation_scan.tsv"), encoding="utf-8"),
        delimiter="\t"))
    cells = []
    for r in rows:
        if r["cell"] not in cells:
            cells.append(r["cell"])

    out = [START, "|---|---|---|---|---|---|---|"]
    for c in cells:
        by = {r["tolerance"]: r for r in rows if r["cell"] == c}
        line, ps = f"| {c} |", []
        for t in TOLS:
            r = by.get(t)
            if r is None:
                line += " — |"
                continue
            cov = float(r["matched_fraction"])
            if r["ok"].strip().upper() != "TRUE":
                line += f" NA（{cov * 100:.0f}%）|"
                continue
            p, f = float(r["empirical_p"]), float(r["fold_vs_null"])
            ps.append(p)
            body = f"{f:.2f}×, {fmt(p)}"
            line += f" **{body}** |" if t == "1" else f" {body} |"
        if ps and all(p < 0.05 for p in ps):
            v = VERDICT_OK
        elif ps and all(p >= 0.05 for p in ps):
            v = VERDICT_NO
        else:
            v = VERDICT_MIX
        out.append(line + f" {v} |")

    p = os.path.join(MR, "manuscript", "PREREG_permutation_estimand.md")
    s = io.open(p, encoding="utf-8").read()
    i = s.index(START)
    j = s.index("\n\n", i)
    io.open(p, "w", encoding="utf-8", newline="\n").write(
        s[:i] + "\n".join(out) + s[j:])
    print(f"rendered S38 §7.1: {len(cells)} cells x {len(TOLS)} tolerances")
    for line in out[2:]:
        print("  ", line[:110])


if __name__ == "__main__":
    main()
