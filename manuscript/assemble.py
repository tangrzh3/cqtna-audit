"""Assemble the submission-shaped manuscript from its two hand-edited sources.

WHY THIS EXISTS
Reviewer 3 objected that Methods sat in a separate file behind a placeholder, so
the manuscript was not a complete submission document. The obvious fix -- paste
Methods into the manuscript -- would create two copies of the same text, and this
project has already been bitten twice by exactly that (Supplementary numbering in
Step 82, withdrawn claims surviving in figure scripts in Step 88).

So the two sources stay single-copy and hand-edited:
    MANUSCRIPT_v2_dual_thread.md   main text, with a §5 Methods pointer
    METHODS_draft.md               the full Methods, the only copy
and this script generates the combined document. The output is DERIVED --
never edit MANUSCRIPT_assembled.md by hand; edit a source and re-run.

Usage:  python assemble.py
"""
import io
import os
import re
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
MAIN = os.path.join(HERE, "MANUSCRIPT_v2_dual_thread.md")
METH = os.path.join(HERE, "METHODS_draft.md")
OUT = os.path.join(HERE, "MANUSCRIPT_assembled.md")

BANNER = (
    "<!-- GENERATED FILE - do not edit by hand.\n"
    "     Sources: MANUSCRIPT_v2_dual_thread.md + METHODS_draft.md\n"
    f"     Regenerate with: python assemble.py   (last built {date.today()}) -->\n\n"
)


def main():
    main_text = io.open(MAIN, encoding="utf-8").read()
    methods = io.open(METH, encoding="utf-8").read()

    # the Methods file carries its own title line; drop it and demote its headings
    # one level so they nest under the manuscript's "## 5. Methods"
    lines = methods.split("\n")
    if lines and lines[0].startswith("# "):
        lines = lines[1:]
    body = []
    for l in lines:
        m = re.match(r"^(#{2,3}) (\d+[a-z]?\.\s*)(.*)$", l)
        if m:
            body.append(f"{m.group(1)}# 5.{m.group(2).rstrip('. ')} {m.group(3)}")
        else:
            body.append(l)
    methods_body = "\n".join(body).strip()

    # replace the §5 pointer block with the real thing
    start = main_text.index("## 5. Methods")
    end = main_text.index("## 6. Figure legends")
    pointer = main_text[start:end]
    assembled = (
        main_text[:start]
        + "## 5. Methods\n\n"
        + methods_body
        + "\n\n---\n\n"
        + main_text[end:]
    )

    # references: splice REFERENCES.md in place of the pointer paragraph, so the
    # assembled document is submission-shaped and the bibliography stays single-copy
    refs_path = os.path.join(HERE, "REFERENCES.md")
    if os.path.exists(refs_path):
        rl = io.open(refs_path, encoding="utf-8").read().split("\n")
        if rl and rl[0].startswith("# "):
            rl = rl[1:]
        refs = "\n".join(rl).strip()
        i = assembled.index("## 8. References")
        tail = assembled[i + 20:]
        k = tail.find("\n## ")
        j2 = (i + 20 + k) if k >= 0 else len(assembled)
        assembled = (assembled[:i] + "## 8. References\n\n"
                     + refs + "\n" + assembled[j2:])
    io.open(OUT, "w", encoding="utf-8").write(BANNER + assembled)
    nw = len(assembled.split())
    print(f"wrote {os.path.basename(OUT)}: {len(assembled.splitlines())} lines, ~{nw:,} words")
    print(f"  (replaced a {len(pointer.split())}-word Methods pointer with "
          f"{len(methods_body.split()):,} words of Methods)")


if __name__ == "__main__":
    main()
