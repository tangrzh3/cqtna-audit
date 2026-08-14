from __future__ import annotations

import csv
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


API = (
    "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/"
    "pmcoa.cgi/BioC_json/{ids}/unicode"
)

C1_TERMS = {
    r"\bknown (?:risk )?loc(?:us|i)\b": 5,
    r"\bpreviously (?:reported|identified|established)\b": 4,
    r"\bestablished (?:risk )?loc(?:us|i)\b": 5,
    r"\breported (?:risk )?loc(?:us|i)\b": 4,
    r"\bGWAS catalog\b": 5,
    r"\b(?:overlap|overlapped|overlapping)\b": 2,
    r"\bnovel (?:risk )?loc(?:us|i)\b": 2,
    r"\bGWAS[- ]identified\b": 3,
    r"\brisk (?:gene|variant|SNP)s?\b": 1,
}

C3_TERMS = {
    r"\b(?:another|second|independent|external|alternative|different) (?:outcome )?(?:GWAS|dataset|cohort)\b": 6,
    r"\b(?:validation|replication) (?:outcome )?(?:GWAS|dataset|cohort)\b": 5,
    r"\brepeat(?:ed|ing)? (?:the )?analysis\b": 4,
    r"\b(?:replicat|validat)(?:e|ed|es|ing|ion)\b": 2,
    r"\bgreater power of the outcome GWAS\b": 7,
    r"\bpower of (?:the )?outcome GWAS\b": 6,
    r"\bmore (?:gene|protein|candidate|target)s?\b": 1,
    r"\b(?:still|remained) (?:significant|meaningful|consistent|robust)\b": 3,
    r"\bFinnGen\b": 1,
    r"\bUK Biobank\b": 1,
}


def fetch_batch(pmcids: list[str], attempts: int = 3):
    url = API.format(ids=",".join(pmcids))
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Codex blind-coding audit/1.0"},
    )
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt < attempts:
                time.sleep(attempt * 2)
    raise RuntimeError(f"BioC retrieval failed for {pmcids}: {last_error}")


def as_collections(payload):
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        return [payload]
    return []


def passage_records(collection):
    for document in collection.get("documents", []):
        for passage in document.get("passages", []):
            text = re.sub(r"\s+", " ", passage.get("text", "")).strip()
            if not text:
                continue
            infons = passage.get("infons") or {}
            yield {
                "section": infons.get("section_type", ""),
                "title": infons.get("title", ""),
                "text": text,
                "infons": infons,
            }


def score(text: str, terms: dict[str, int]) -> int:
    return sum(weight for pattern, weight in terms.items() if re.search(pattern, text, re.I))


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: 105a_research_fetch.py INPUT_TSV OUTPUT_DIR")
    source = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    text_dir = output_dir / "fulltext"
    output_dir.mkdir(parents=True, exist_ok=True)
    text_dir.mkdir(parents=True, exist_ok=True)

    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    pmcids = sorted({row["pmc"].strip() for row in rows if row["pmc"].strip()})

    collections = []
    errors = []
    for start in range(0, len(pmcids), 8):
        batch = pmcids[start : start + 8]
        try:
            collections.extend(as_collections(fetch_batch(batch)))
        except RuntimeError as exc:
            errors.append(str(exc))
            for pmcid in batch:
                try:
                    collections.extend(as_collections(fetch_batch([pmcid])))
                except RuntimeError as single_exc:
                    errors.append(str(single_exc))
        time.sleep(0.5)

    articles = {}
    for collection in collections:
        passages = list(passage_records(collection))
        if not passages:
            continue
        infons = passages[0]["infons"]
        pmcid = infons.get("article-id_pmc", "")
        if not pmcid:
            for item in passages:
                pmcid = item["infons"].get("article-id_pmc", "")
                if pmcid:
                    break
        if not pmcid:
            continue
        title = next((p["text"] for p in passages if p["section"] == "TITLE"), "")
        doi = next((p["infons"].get("article-id_doi", "") for p in passages if p["infons"].get("article-id_doi")), "")
        articles[pmcid] = {
            "pmc": pmcid,
            "title": title,
            "doi": doi,
            "url": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/",
            "passages": passages,
        }

    with (output_dir / "fulltexts.json").open("w", encoding="utf-8") as handle:
        json.dump(articles, handle, ensure_ascii=False, indent=2)

    with (output_dir / "article_index.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["pmc", "title", "doi", "url", "passage_count"])
        for pmcid in pmcids:
            article = articles.get(pmcid, {})
            writer.writerow([
                pmcid,
                article.get("title", ""),
                article.get("doi", ""),
                article.get("url", ""),
                len(article.get("passages", [])),
            ])

    with (output_dir / "candidate_passages.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["pmc", "criterion", "rank", "score", "section", "title", "text"])
        for pmcid in pmcids:
            article = articles.get(pmcid)
            if not article:
                continue
            clean_lines = []
            for idx, passage in enumerate(article["passages"], start=1):
                label = passage["section"] or passage["title"] or "TEXT"
                clean_lines.append(f"[{idx:04d}] [{label}] {passage['text']}")
            (text_dir / f"{pmcid}.txt").write_text("\n\n".join(clean_lines), encoding="utf-8")

            for criterion, terms in (("C1_known_locus", C1_TERMS), ("C3_power_stability", C3_TERMS)):
                scored = []
                for idx, passage in enumerate(article["passages"]):
                    value = score(passage["text"], terms)
                    if value:
                        scored.append((value, idx, passage))
                scored.sort(key=lambda item: (-item[0], item[1]))
                for rank, (value, _, passage) in enumerate(scored[:14], start=1):
                    writer.writerow([
                        pmcid,
                        criterion,
                        rank,
                        value,
                        passage["section"],
                        passage["title"],
                        passage["text"],
                    ])

    missing = sorted(set(pmcids) - set(articles))
    summary = {
        "requested": len(pmcids),
        "retrieved": len(articles),
        "missing": missing,
        "errors": errors,
    }
    (output_dir / "retrieval_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False))
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
