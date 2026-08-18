## Semi-external validation: does the tool recover an enrichment we planted, and
## does it stay quiet when there is none?
##
## Every other check in this package compares against the study it came from, so
## none of them can tell whether the machinery is right -- only whether it is
## unchanged. This one knows the truth because it made it up.
##
##   Rscript inst/validate_synthetic.R

suppressMessages(library(cqtna))
set.seed(20260818)

simulate <- function(n_loci = 500, n_sig = 20, frac_known_bg = 0.10,
                     frac_known_sig, records_per_locus = 5, locus_kb = 100) {
  # one locus per chromosome-like block, well separated so nothing chains
  starts <- seq_len(n_loci) * 10e6
  chr <- as.character(rep(seq_len(22), length.out = n_loci))
  rows <- do.call(rbind, lapply(seq_len(n_loci), function(i) {
    data.frame(locus_i = i, chr = chr[i],
               pos = starts[i] + seq_len(records_per_locus) * 1e4,
               stringsAsFactors = FALSE)
  }))
  sig_loci <- sample(n_loci, n_sig)
  rows$p <- ifelse(rows$locus_i %in% sig_loci, 1e-12, runif(nrow(rows), .2, 1))
  rows$gene <- paste0("G", seq_len(nrow(rows)))
  rows$record_id <- seq_len(nrow(rows))

  # Known status is drawn, not planted exactly: a significant locus is known
  # with probability frac_known_sig, a background locus with frac_known_bg.
  # An earlier version fixed the counts, which made the recovered fold
  # deterministic (sd 0.00) and the calibration meaningless -- the simulation
  # has to carry sampling variability or it is testing arithmetic, not a test.
  is_sig <- seq_len(n_loci) %in% sig_loci
  pk <- ifelse(is_sig, frac_known_sig, frac_known_bg)
  known_loci <- which(stats::runif(n_loci) < pk)
  if (!length(known_loci)) known_loci <- sample(n_loci, 1)
  known <- data.frame(chr = chr[known_loci], pos = starts[known_loci] + 2e4,
                      source = "simulated", stringsAsFactors = FALSE)

  list(mr = as_cqtna_mr(rows[, c("record_id", "gene", "chr", "pos", "p")],
                        locus_kb = locus_kb, build = "SIM"),
       known = as_cqtna_known(known, build = "SIM"),
       truth_fold = frac_known_sig / frac_known_bg)
}

report <- function(label, frac_known_sig, reps = 100, n_perm = 100) {
  fold <- p_fisher <- p_perm <- numeric(reps)
  for (i in seq_len(reps)) {
    s <- simulate(frac_known_sig = frac_known_sig)
    a <- cqtna_attribution(s$mr, s$known)
    fold[i] <- a$fold
    p_fisher[i] <- a$fisher_p_one_sided
    p_perm[i] <- cqtna_permutation_control(s$mr, s$known, n_perm = n_perm)$empirical_p
  }
  truth <- simulate(frac_known_sig = frac_known_sig)$truth_fold
  cat(sprintf("%-26s planted %4.2fx | recovered %4.2fx (sd %.2f) | Fisher P<0.05 %5.1f%% | perm P<0.05 %5.1f%%\n",
              label, truth, mean(fold), stats::sd(fold),
              100 * mean(p_fisher < .05), 100 * mean(p_perm < .05)))
}

cat("Synthetic validation -- 500 loci, 20 significant, 10% known background\n")
cat(strrep("-", 118), "\n")
report("null (no enrichment)",        0.10)   # significant loci known at background rate
report("modest (2x)",                 0.20)
report("strong (4x)",                 0.40)
cat(strrep("-", 118), "\n")
cat("Read the null row as calibration: 'Fisher P<0.05' there is the false-positive rate.\n")
cat("Read the other rows as power. 'recovered' should track 'planted'.\n")
