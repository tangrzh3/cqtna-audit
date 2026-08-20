## Step 130 -- the negative control, run against several unrelated diseases.
##
## Executes manuscript/PREREG_multilist_control.md (S39), frozen before any list
## was scored against any cell.
##
## Why: S36 showed the single mismatched list is not orthogonal for RA -- the
## control fails at six of seven window settings because RA and melanoma really
## do share immune loci. A single list cannot separate "the attribution is
## specific" from "we happened to pick a disease that shares structure". Scoring
## every cell against several pre-selected unrelated lists removes that degree of
## freedom: the reported quantity becomes how many unrelated lists also enrich,
## rather than the verdict of the one list we chose.
##
## Lists are those already in the tree with >= 30 GRCh38-placeable lead SNPs
## (S39 §3). No new data is fetched, so a disappointing result cannot be followed
## by "let us add a few more lists".
##
##   Rscript step130_multilist_control.R
## Output: 130a_multilist_main.tsv    every cell x every list, main setting
##         130b_multilist_sweep.tsv   the same over four window widths
##         130c_multilist_verdict.tsv one row per cell: k, F_own, F_max, verdict

MR <- "D:/R_ex/MR"
setwd(MR)
source("step124_cells.R")   # cells, runs, status_of()

KB       <- 1000
CONV     <- "any_record"
MIN_LEAD <- 30              # S34's floor, reused; see S39 §3
SWEEP    <- c(100, 250, 500, 1000)

## ---------------------------------------------------------------- the lists
## Loaded by rule, not by hand, so the inclusion criterion is executed rather
## than asserted. A list that drops below the floor drops out visibly.
list_files <- c(melanoma   = "landi2020_known_loci_grch38.csv",
                HCC        = "84a_hcc_known_loci_grch38.csv",
                colorectal = "known_loci_colorectal_grch38.csv",
                prostate   = "known_loci_prostate_grch38.csv",
                breast     = "known_loci_breast_grch38.csv",
                lung       = "known_loci_lung_grch38.csv")

known <- list(); sizes <- integer(0)
for (nm in names(list_files)) {
  f <- list_files[[nm]]
  if (!file.exists(f)) { sizes[nm] <- 0L; next }
  d <- tryCatch(read.csv(f, stringsAsFactors = FALSE),
                error = function(e) data.frame())
  d <- d[!is.na(suppressWarnings(as.numeric(d$pos))), , drop = FALSE]
  sizes[nm] <- nrow(d)
  if (nrow(d) < MIN_LEAD) next
  d$source <- nm
  known[[nm]] <- as_cqtna_known(d[, c("chr", "pos", "source")], build = "GRCh38")
}
cat("known-locus lists (floor =", MIN_LEAD, "lead SNPs):\n")
for (nm in names(sizes))
  cat(sprintf("  %-11s %5d  %s\n", nm, sizes[nm],
              if (nm %in% names(known)) "included" else "EXCLUDED (below floor)"))

## RA's reference list is the placed Okada positions, which live in 123b rather
## than in a csv, so it is added here rather than by the loop above.
ra_known <- read.delim("123b_ra_known_positions.tsv", stringsAsFactors = FALSE)
known[["RA"]] <- as_cqtna_known(within(ra_known, source <- "Okada2014"),
                                build = "GRCh38")
cat(sprintf("  %-11s %5d  included (own list for the RA cells)\n",
            "RA", nrow(ra_known)))

## which list is a cell's OWN outcome, and therefore not one of its controls
own_of <- function(cell) {
  if (startsWith(cell, "melanoma")) "melanoma"
  else if (startsWith(cell, "HCC")) "HCC"
  else if (startsWith(cell, "RA")) "RA"
  else NA_character_
}

## S39 §3: the six registered cells and the two power sensitivities, whole
## genome only. The post-hoc MHC-excluded rows are not part of this test.
use <- Filter(function(cl) cl$region == "all_genome", cells)

score <- function(d, locus_kb, known_kb) {
  mr <- as_cqtna_mr(d, locus_kb = locus_kb, build = "GRCh38",
                    locus_method = "fixed_centre")
  out <- list()
  for (nm in names(known)) {
    a <- cqtna_attribution(mr, known[[nm]], known_kb = known_kb,
                           known_from = CONV)
    out[[nm]] <- data.frame(
      list_name = nm, bg_loci = a$background_loci, bg_known = a$background_known,
      sig_loci = a$significant_loci, sig_known = a$significant_known,
      fold = round(a$fold, 2), fisher_p = a$fisher_p_one_sided,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

## ------------------------------------------------------------ A. main setting
rows <- list()
for (cl in use) {
  own <- own_of(cl$name)
  s <- score(cl$d, KB, KB)
  s$cell <- cl$name
  s$role <- cl$role
  s$relation <- ifelse(s$list_name == own, "own outcome", "unrelated control")
  rows[[length(rows) + 1L]] <- s
}
main <- do.call(rbind, rows)
main <- main[, c("cell", "role", "list_name", "relation", setdiff(
  names(main), c("cell", "role", "list_name", "relation")))]
write.table(main, "130a_multilist_main.tsv", sep = "\t", row.names = FALSE,
            quote = FALSE)

## ------------------------------------------------------------- B. the verdict
verdict <- list()
for (cl in use) {
  own <- own_of(cl$name)
  s <- main[main$cell == cl$name, ]
  o <- s[s$list_name == own, ]
  ctl <- s[s$relation == "unrelated control", ]
  hit <- ctl[is.finite(ctl$fisher_p) & ctl$fisher_p < 0.05 &
             is.finite(ctl$fold) & ctl$fold > 1, , drop = FALSE]
  k <- nrow(hit)
  f_own <- o$fold[1]
  ## Two margins, because they answer different questions and only one of them
  ## is safe to quote as a ratio.
  ##   F_max_sig  largest fold among controls that reach P < 0.05. This is what
  ##              the S39 section 5 verdict rule uses. It is 0 when no control
  ##              is significant, so the ratio against it is undefined -- NOT
  ##              "an unbounded margin", which is how an earlier version of the
  ##              figure and text described it. A control with no significant
  ##              enrichment still has a point estimate.
  ##   F_max_all  largest fold among ALL controls, significant or not. Always
  ##              finite, so margin_all is the ratio to report.
  f_max <- if (k) max(hit$fold) else 0
  f_max_all <- suppressWarnings(max(ctl$fold[is.finite(ctl$fold)]))
  if (!is.finite(f_max_all)) f_max_all <- NA_real_
  v <- if (!is.finite(o$sig_loci[1]) || o$sig_loci[1] < 2) "D no power"
       else if (k == 0) "A control clean"
       else if (is.finite(f_own) && f_own > f_max) "B not orthogonal, direction separates"
       else "C attribution claim fails"
  verdict[[length(verdict) + 1L]] <- data.frame(
    cell = cl$name, role = cl$role, own_list = own,
    n_controls = nrow(ctl), k_enriching = k,
    F_own = f_own, F_max_sig = f_max, F_max_all = f_max_all,
    margin_all = round(f_own / f_max_all, 2),
    best_rival_any = ctl$list_name[which.max(ctl$fold)][1],
    most_favourable = ctl$list_name[which.min(ctl$fold)][1],
    which_enrich = if (k) paste(sprintf("%s(%.2f,P=%.3g)", hit$list_name,
                                        hit$fold, hit$fisher_p),
                                collapse = "; ") else "-",
    verdict = v, stringsAsFactors = FALSE)
}
vd <- do.call(rbind, verdict)
write.table(vd, "130c_multilist_verdict.tsv", sep = "\t", row.names = FALSE,
            quote = FALSE)

## -------------------------------------------------------------- C. the sweep
rows <- list()
for (cl in use) for (kb in SWEEP) {
  for (mode in c("known_kb", "locus_kb")) {
    s <- score(cl$d, if (mode == "locus_kb") kb else KB,
               if (mode == "known_kb") kb else KB)
    s$cell <- cl$name; s$swept <- mode; s$kb <- kb
    s$relation <- ifelse(s$list_name == own_of(cl$name), "own outcome",
                         "unrelated control")
    rows[[length(rows) + 1L]] <- s
  }
}
sw <- do.call(rbind, rows)
write.table(sw, "130b_multilist_sweep.tsv", sep = "\t", row.names = FALSE,
            quote = FALSE)

## ------------------------------------------------------------------- report
cat("\n", strrep("=", 116), "\n",
    "A. every cell against every list -- fixed anchor ", KB, " kb, C1\n", sep = "")
for (nm in unique(main$cell)) {
  s <- main[main$cell == nm, ]
  cat("\n  ", nm, "\n", sep = "")
  cat(sprintf("    %-12s %-18s %10s %10s %8s %11s\n",
              "list", "relation", "bg known", "sig known", "fold", "P"))
  for (i in seq_len(nrow(s))) with(s[i, ], cat(sprintf(
    "    %-12s %-18s %4d/%-5d %4d/%-5d %8s %11.3g%s\n",
    list_name, relation, bg_known, bg_loci, sig_known, sig_loci,
    ifelse(is.na(fold), "-", sprintf("%.2f", fold)), fisher_p,
    ifelse(relation == "unrelated control" & is.finite(fisher_p) &
           fisher_p < 0.05 & fold > 1, "  <-- enriches", ""))))
}

cat("\n", strrep("=", 116), "\n", "B. verdict per cell (S39 §5)\n", sep = "")
cat(sprintf("%-28s %5s %7s %9s %8s  %-11s %-11s %s\n", "cell", "k",
            "F_own", "F_max_all", "margin", "best rival", "most fav.", "verdict"))
for (i in seq_len(nrow(vd))) with(vd[i, ], cat(sprintf(
  "%-28s %5d %7s %9s %8s  %-11s %-11s %s\n", cell, k_enriching,
  ifelse(is.na(F_own), "-", sprintf("%.2f", F_own)),
  ifelse(is.na(F_max_all), "-", sprintf("%.2f", F_max_all)),
  ifelse(is.na(margin_all), "-", sprintf("%.2f", margin_all)),
  best_rival_any, most_favourable, verdict)))

cat("\nwrote 130a / 130b / 130c\n")
