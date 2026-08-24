#!/usr/bin/env bash
# step136 -- S40 前提三：在收到任何外部关联统计之前，冻结并散列全部输入。
#
# S40 §3 要求存放：protocol、可执行容器、renv.lock、CQTNA 源码包、
# 冻结的工具表、已知位点名单、对照名单、检验族清单、SHA-256 清单。
# 本脚本产出其中**本项目自己能产出的部分**，并在 S40_LOCK.md 里
# 逐条写明哪几项产不出来、为什么、以及要谁来补。
#
#   bash step136_lock_manifest.sh
# 产出：136a_lock_manifest.tsv（每行一个文件：类别、路径、字节、SHA-256）
#       136b_environment.tsv（R 版本与实际加载的包版本）
#       136c_console.log
set -euo pipefail
cd "$(dirname "$0")"

MAN=136a_lock_manifest.tsv
printf 'category\tpath\tbytes\tsha256\n' > "$MAN"

add() {                      # add <category> <path...>
  local cat="$1"; shift
  for f in "$@"; do
    if [ ! -f "$f" ]; then
      echo "MISSING: $f" >&2
      printf '%s\t%s\tNA\tMISSING\n' "$cat" "$f" >> "$MAN"
      continue
    fi
    printf '%s\t%s\t%s\t%s\n' "$cat" "$f" \
      "$(wc -c < "$f" | tr -d ' ')" \
      "$(sha256sum "$f" | cut -d' ' -f1)" >> "$MAN"
  done
}

# --- 冻结的记录表（每个格子的暴露×结局记录，工具表即由此决定）
add frozen_records \
  13_meta_locus_annotation.tsv 92c_locus_annotated.tsv \
  85a_HCC_high_annotated.tsv 85a_HCC_low_annotated.tsv \
  123b_ra_records.tsv.gz 123c_eqtlgen_hcc_records.tsv.gz

# --- 已知位点名单：自身名单 + S39 的全部对照名单
add known_lists \
  landi2020_known_loci_grch38.csv 84a_hcc_known_loci_grch38.csv \
  known_loci_colorectal_grch38.csv known_loci_prostate_grch38.csv \
  known_loci_breast_grch38.csv known_loci_lung_grch38.csv \
  123b_ra_known_positions.tsv

# --- 分析代码：共享格子定义 + 冻结分区之后的每一步
add scripts \
  step124_cells.R step124_full_grid.R step125_mismatch_loci.R \
  step126_recompute_on_fixed_anchor.R step127_audit_manuscript_numbers.py \
  step128_window_sensitivity_fixed_anchor.R step129_selfcheck_on_fixed_anchor.R \
  step130_multilist_control.R step131_render_s38_scan.py \
  step132_external_power_gate.R step133_projection_reading_b.R \
  step134_candidate_cohorts.R step135_hcc_gate_distance.R \
  step140_estimator_identity.R step141_enrichment_decomposition.R \
  step145_method_benchmark.R step148_publish_litaudit.py \
  step149_typeI_calibration.R step150_environment_lock.R \
  make_review_packet.sh step136_lock_manifest.sh

# --- R 包源码 tarball
add package cqtna_0.3.0.tar.gz

# --- 决策与预注册记录（检验族由这些定义，不由代码定义）
add decision_records \
  EXTERNAL_VALIDATION_PROTOCOL_v9_3.md \
  manuscript/PREREG_external_validation.md \
  manuscript/PREREG_locus_partition.md \
  manuscript/PREREG_permutation_estimand.md \
  manuscript/PREREG_multilist_control.md \
  manuscript/AMENDMENT_S38_permutation_order.md \
  S40_POWER_GATE_SPEC.md S40_COHORT_SCREEN.md CONVENTION_REVIEW.md

# --- 结果表：锁定当下的值，日后任何改动都能被看见
add result_tables \
  123d_fixed_anchor_full_grid.tsv 126b_permutation_scan.tsv \
  126c_permutation_primary.tsv 130a_multilist_main.tsv 130c_multilist_verdict.tsv \
  132a_power_gate_curve.tsv 132b_exclusion_roster.tsv \
  133a_projection_reading_b.tsv 134a_candidate_evaluation.tsv \
  135a_hcc_gate_distance.tsv

# --- 环境
Rscript -e '
  ps <- c("cqtna","stats","utils","graphics")
  out <- data.frame(component = c("R", ps),
                    version = c(paste(R.version$major, R.version$minor, sep="."),
                                sapply(ps, function(p)
                                  tryCatch(as.character(packageVersion(p)),
                                           error = function(e) "NOT INSTALLED"))),
                    stringsAsFactors = FALSE)
  out <- rbind(out, data.frame(component = "platform", version = R.version$platform))
  write.table(out, "136b_environment.tsv", sep="\t", row.names=FALSE, quote=FALSE)
' > /dev/null
add environment 136b_environment.tsv 150a_environment.lock 150b_session.txt

# --- 清单自身的散列，以及 git 状态
SELF=$(sha256sum "$MAN" | cut -d' ' -f1)
N=$(( $(wc -l < "$MAN") - 1 ))
MISS=$(grep -c 'MISSING' "$MAN" || true)
COMMIT=$(git rev-parse HEAD)
DIRTY=$(git status --porcelain --untracked-files=no | wc -l | tr -d ' ')

{
  echo "files hashed          : $N"
  echo "missing               : $MISS"
  echo "manifest SHA-256      : $SELF"
  echo "git commit            : $COMMIT"
  echo "tracked files modified: $DIRTY"
  echo
  echo "⚠ 这份清单锁的是**输入与代码**，不是一个可执行容器。"
  echo "⚠ 容器镜像与 renv.lock 仍然缺席 —— 见 S40_LOCK.md §2。"
} | tee 136c_console.log
