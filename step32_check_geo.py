"""Verify candidate GEO series before committing to any download."""
import gzip
import re
import sys
import urllib.request

KEYS = ("!Series_title", "!Series_summary", "!Series_overall_design",
        "!Series_type", "!Series_platform_id", "!Sample_title",
        "!Sample_source_name", "!Sample_characteristics",
        "!Sample_library_strategy", "!Sample_organism")


def check(acc):
    stem = acc[:-3] + "nnn"
    print("=" * 80)
    print(acc)
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{acc}/matrix/"
           f"{acc}_series_matrix.txt.gz")
    try:
        raw = urllib.request.urlopen(url, timeout=120).read()
        txt = gzip.decompress(raw).decode("utf-8", "replace")
    except Exception as ex:
        print(f"  matrix FAIL {type(ex).__name__}: {ex}")
        txt = ""
    for line in txt.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith(KEYS):
            s = line.strip()
            print(s[:900] + ("..." if len(s) > 900 else ""))
    print("\n  -- supplementary files --")
    try:
        idx = urllib.request.urlopen(
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{acc}/suppl/",
            timeout=120).read().decode("utf-8", "replace")
        files = [m for m in re.findall(r'href="([^"]+)"', idx)
                 if not m.startswith(("/", "http")) and m != "../"]
        for f in files[:25]:
            try:
                u = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{acc}/suppl/{f}"
                r = urllib.request.urlopen(urllib.request.Request(u, method="HEAD"),
                                           timeout=60)
                mb = int(r.headers.get("Content-Length", 0)) / 1e6
                print(f"     {f:<60}{mb:9.1f} MB")
            except Exception:
                print(f"     {f}")
        if len(files) > 25:
            print(f"     ... and {len(files)-25} more")
    except Exception as ex:
        print(f"    suppl FAIL {type(ex).__name__}: {ex}")
    print()


for acc in (sys.argv[1:] or ["GSE166188", "GSE282266"]):
    check(acc)
