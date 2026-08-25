#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Parse MANUSCRIPT_GB.md into a JSON block list for the docx builder."""
import io
import json
import re
import sys

src, dest = sys.argv[1], sys.argv[2]
lines = io.open(src, encoding="utf-8").read().split("\n")

blocks = []
i = 0
n = len(lines)


def flush(buf):
    if buf:
        blocks.append({"t": "p", "text": " ".join(buf).strip()})
    return []


buf = []
while i < n:
    l = lines[i]

    # headings
    m = re.match(r"^(#{1,4})\s+(.*)$", l)
    if m:
        buf = flush(buf)
        blocks.append({"t": "h", "level": len(m.group(1)), "text": m.group(2).strip()})
        i += 1
        continue

    # horizontal rule
    if re.match(r"^---+\s*$", l):
        buf = flush(buf)
        i += 1
        continue

    # table: a line with | followed by a separator row
    if l.strip().startswith("|") and i + 1 < n and re.match(r"^\s*\|[\s:|-]+\|\s*$", lines[i + 1]):
        buf = flush(buf)
        rows = []
        while i < n and lines[i].strip().startswith("|"):
            raw = lines[i].strip().strip("|")
            if not re.match(r"^[\s:|-]+$", raw):
                rows.append([c.strip() for c in raw.split("|")])
            i += 1
        blocks.append({"t": "table", "rows": rows})
        continue

    # bullet list
    if re.match(r"^\s*[-*]\s+", l):
        buf = flush(buf)
        items = []
        while i < n and re.match(r"^\s*[-*]\s+", lines[i]):
            items.append(re.sub(r"^\s*[-*]\s+", "", lines[i]).strip())
            i += 1
        blocks.append({"t": "ul", "items": items})
        continue

    # blank line ends a paragraph
    if not l.strip():
        buf = flush(buf)
        i += 1
        continue

    buf.append(l.strip())
    i += 1

flush(buf)

io.open(dest, "w", encoding="utf-8").write(
    json.dumps(blocks, ensure_ascii=False, indent=1))
print("parsed %d blocks (%d headings, %d tables, %d paragraphs)"
      % (len(blocks),
         sum(1 for b in blocks if b["t"] == "h"),
         sum(1 for b in blocks if b["t"] == "table"),
         sum(1 for b in blocks if b["t"] == "p")))
