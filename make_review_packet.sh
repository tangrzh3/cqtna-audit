#!/usr/bin/env bash
# Assemble the third-party review packet from tracked sources.
#
# The packet is a build product. Tracking its contents as well as their originals
# means two copies of every table, drifting apart the moment one is regenerated
# -- the same reason .Rcheck output and the zip are ignored. Only the two
# hand-written READMEs inside the packet are tracked, because nothing else
# produces them.
#
#   bash make_review_packet.sh
# Produces: review_packet_v9/ and review_packet_v9.zip
set -euo pipefail
cd "$(dirname "$0")"

OUT=review_packet_v9
REP="$OUT/reproduce"
mkdir -p "$REP"

# --- manuscript and the decision records the attribution line rests on
cp manuscript/MANUSCRIPT_GB.md \
   manuscript/PREREG_locus_partition.md \
   manuscript/PREREG_permutation_estimand.md \
   manuscript/PREREG_multilist_control.md \
   manuscript/PREREG_noncancer_outcome.md \
   manuscript/SUPP_window_sensitivity.md \
   manuscript/NUMBER_MIGRATION_fixed_anchor.md \
   manuscript/AMENDMENT_S38_permutation_order.md \
   CONVENTION_REVIEW.md "$OUT/"

# --- every table computed on the frozen partition
cp 123d_fixed_anchor_full_grid.tsv \
   125a_mismatch_locus_diagnostics.tsv 125b_mismatch_spread_summary.tsv \
   126a_offgrid_attribution.tsv 126b_permutation_scan.tsv 126c_permutation_primary.tsv \
   128a_known_kb_sweep.tsv 128b_locus_kb_sweep.tsv 128c_distances.tsv \
   129a_selfcheck_partition_comparison.tsv \
   130a_multilist_main.tsv 130b_multilist_sweep.tsv 130c_multilist_verdict.tsv \
   85e_matched_background_fixed_anchor.tsv "$OUT/"

# --- figures
cp figures/Fig1_locus_attribution.pdf figures/Fig2_generality.pdf \
   figures/FigS_window_control.pdf figures/FigS_multilist_control.pdf "$OUT/"

# --- everything needed to recompute those tables from per-record data
cp 123b_ra_records.tsv.gz 123c_eqtlgen_hcc_records.tsv.gz 99d_r13_records.tsv.gz \
   123b_ra_known_positions.tsv \
   landi2020_known_loci_grch38.csv 84a_hcc_known_loci_grch38.csv \
   known_loci_colorectal_grch38.csv known_loci_prostate_grch38.csv \
   13_meta_locus_annotation.tsv 92c_locus_annotated.tsv \
   85a_HCC_high_annotated.tsv 85a_HCC_low_annotated.tsv \
   116b_mr_Nathan_2022.tsv 116b_mr_Randolph_2021.tsv \
   116b_mr_Schmiedel_2018.tsv \
   known_loci_breast_grch38.csv known_loci_lung_grch38.csv \
   96a_attribution_selfcheck.tsv "$REP/"
cp step124_cells.R step124_full_grid.R step125_mismatch_loci.R \
   step126_recompute_on_fixed_anchor.R step127_audit_manuscript_numbers.py \
   step128_window_sensitivity_fixed_anchor.R step129_selfcheck_on_fixed_anchor.R \
   step130_multilist_control.R "$REP/"

# --- the package, rebuilt so the tarball cannot lag its sources
R_BIN="${R_BIN:-/d/R/R-4.4.1/bin/R.exe}"
rm -f cqtna_0.3.0.tar.gz
if [ -x "$R_BIN" ]; then
  "$R_BIN" CMD build cqtna_r >/dev/null
else
  R CMD build cqtna_r >/dev/null
fi
cp cqtna_0.3.0.tar.gz "$OUT/"

rm -f review_packet_v9.zip
if command -v zip >/dev/null 2>&1; then
  zip -qr review_packet_v9.zip "$OUT"
else
  powershell -NoProfile -Command \
    "Compress-Archive -Path 'review_packet_v9\\*' -DestinationPath 'review_packet_v9.zip' -Force"
fi

echo "packet: $(find "$OUT" -type f | wc -l) files"
echo "zip:    $(du -h review_packet_v9.zip | cut -f1)"
echo
echo "The two READMEs in the packet are hand-written and tracked; everything else"
echo "is a copy. Edit the originals, not the copies."

# --- smoke test: every file the reproduce scripts READ must be present.
# --- Matching every filename in the source also caught outputs and regex
# --- debris, which is noise; only read calls are inputs.
missing=0
inputs=$(grep -ohE '(read[.](delim|csv|table)|gzfile|file[.]exists)[(]"[^"]+"' \
           "$REP"/step1*.R "$REP"/step1*.py 2>/dev/null |
         grep -oE '"[^"]+"' | tr -d '"' | grep -E '[.](tsv|csv|gz)$' | sort -u)
for f in $inputs; do
  if [ ! -e "$REP/$f" ]; then
    echo "  MISSING from reproduce/: $f"
    missing=1
  fi
done
if [ "$missing" -eq 0 ]; then
  echo "reproduce/: every input the scripts read is present ($(echo "$inputs" | wc -w) files)"
else
  echo "reproduce/: INCOMPLETE -- the README promises a rebuild that will fail"
  exit 1
fi
