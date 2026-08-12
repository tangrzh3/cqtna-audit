"""Verify what GSE199994 actually contains before designing anything around it."""
import gzip
import io
import urllib.request
import re

URLS = [
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199994/matrix/GSE199994_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199994/soft/GSE199994_family.soft.gz",
]
KEYS = ("!Series_title", "!Series_summary", "!Series_overall_design",
        "!Series_type", "!Sample_title", "!Sample_source_name",
        "!Sample_characteristics", "!Sample_organism", "!Sample_library_strategy",
        "!Sample_geo_accession", "!Series_supplementary_file",
        "!Sample_supplementary_file")

for url in URLS:
    print("=" * 78)
    print(url)
    try:
        raw = urllib.request.urlopen(url, timeout=120).read()
    except Exception as ex:
        print(f"  FAIL {type(ex).__name__}: {ex}")
        continue
    txt = gzip.decompress(raw).decode("utf-8", "replace")
    for line in txt.splitlines():
        if line.startswith(KEYS):
            s = line.strip()
            print(s[:600] + ("..." if len(s) > 600 else ""))
        if line.startswith("!series_matrix_table_begin"):
            break
    break

# supplementary file listing
print("\n" + "=" * 78)
print("supplementary files on the FTP")
try:
    idx = urllib.request.urlopen(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199994/suppl/",
        timeout=120).read().decode("utf-8", "replace")
    for m in re.findall(r'href="([^"]+)"', idx):
        if not m.startswith("/") and m not in ("../",):
            print("  ", m)
except Exception as ex:
    print(f"  FAIL {type(ex).__name__}: {ex}")
