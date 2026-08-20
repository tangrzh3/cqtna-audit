## step134 -- 把 §7 的闸门（读法 B）直接套在具体候选队列上。
##
## step133 扫的是 R 的网格；本脚本改为按每个候选队列自己的 N_eff 算出 g，
## 再逐个报出它会拿到几个显著位点、fold 多少、Fisher P 多少。
## 判据、分区、归属口径全部沿用 step133，未作任何改动。
##
## 候选队列的病例/对照数来自 132b 之外的公开记录，逐条在 134a 里带出处。
##
## 输出：134a_candidate_evaluation.tsv、134b_console.log

suppressMessages(library(cqtna))

arg <- commandArgs(trailingOnly = TRUE)
root <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
        if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else {
          a <- commandArgs(FALSE)
          f <- sub("^--file=", "", a[grep("^--file=", a)])
          if (length(f)) dirname(normalizePath(f)) else getwd()
        }
setwd(root)

FDR <- 0.05; LOCUS_KB <- 1000; KNOWN_KB <- 1000
neff <- function(a, b) 4 / (1 / a + 1 / b)
N_EFF_DISC <- neff(12530, 789099)

fix_gene <- function(g) {
  g <- as.character(g); bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad)); g
}
cd4 <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
cd4$chr <- sub(":.*", "", cd4$SNP); cd4$pos <- as.numeric(sub(".*:", "", cd4$SNP))
d <- data.frame(record_id = seq_len(nrow(cd4)), gene = fix_gene(cd4$SYMBOL),
                chr = cd4$chr, pos = cd4$pos, p = cd4$pval, stringsAsFactors = FALSE)
keep <- is.finite(d$p) & is.finite(d$pos) & d$pos > 0
d <- d[keep, , drop = FALSE]; cd4 <- cd4[keep, , drop = FALSE]
mel <- read.csv("landi2020_known_loci_grch38.csv", stringsAsFactors = FALSE)
mel$source <- "Landi2020"; K <- as_cqtna_known(mel, build = "GRCh38")

vkey <- paste0(cd4$chr, ":", cd4$pos)
VAR <- factor(vkey, levels = unique(vkey[order(
  suppressWarnings(as.numeric(sub("X", "23", cd4$chr))), cd4$pos)]))
vidx <- as.integer(VAR)
Z_DISC <- cd4$beta_outcome / cd4$se_outcome
ZV <- as.numeric(tapply(Z_DISC, vidx, function(z) z[1]))

run_g <- function(g) {
  dd <- d; dd$p <- 2 * stats::pnorm(-abs(g * ZV))[vidx]
  mrx <- as_cqtna_mr(dd, locus_kb = LOCUS_KB, build = "GRCh38",
                     locus_method = "fixed_centre")
  ax <- cqtna_attribution(mrx, K, known_kb = KNOWN_KB, fdr = FDR,
                          known_from = "any_record")
  c(ST = ax$significant_loci, SK = ax$significant_known,
    fold = ax$fold, p = ax$fisher_p_one_sided)
}

## 自校准：g = 1 必须逐项复现发现集。
v <- run_g(1); stopifnot(v[["ST"]] == 8, v[["SK"]] == 4)
cat(sprintf("self-calibration g=1: ST %d SK %d fold %.3f P %.5g\n",
            v[["ST"]], v[["SK"]], v[["fold"]], v[["p"]]))

cand <- data.frame(
  cohort = c("discovery meta (FinnGen R12 + Rashkin)",
             "MVP PheCode 172.1 (dx or hx)",
             "MVP PheCode 172.11",
             "MVP histologically confirmed invasive"),
  definition = c("reference", "phecode", "phecode", "histology"),
  eligible_s2_2 = c(NA, FALSE, FALSE, TRUE),
  cases = c(12530, 45274, 10468, 5364),
  controls = c(789099, 389597, 436283, 436283),
  evidence = c("manuscript Methods",
               "GCST90475577 (Verma 2024, PMID 39024449)",
               "GCST90475578 (Verma 2024, PMID 39024449)",
               "Wheless 2025, PMID 39853431; control pool assumed as GCST90475578"),
  stringsAsFactors = FALSE)
cand$n_eff <- neff(cand$cases, cand$controls)
cand$R <- cand$n_eff / N_EFF_DISC

out <- list()
for (f in c(0.50, 0.75, 1.00)) for (i in seq_len(nrow(cand))) {
  g <- f * sqrt(cand$R[i]); v <- run_g(g)
  out[[length(out) + 1]] <- data.frame(
    cand[i, c("cohort", "definition", "eligible_s2_2", "cases", "controls",
              "evidence")],
    n_eff = round(cand$n_eff[i]), R = cand$R[i], retention = f, g = g,
    sig_loci = v[["ST"]], sig_known = v[["SK"]], fold = v[["fold"]],
    fisher_p = v[["p"]],
    passes = v[["ST"]] >= 8 && is.finite(v[["fold"]]) && v[["fold"]] > 1 &&
             is.finite(v[["p"]]) && v[["p"]] < 0.05)
}
res <- do.call(rbind, out)
res <- res[order(res$retention, -res$R), ]
for (i in seq_len(nrow(res)))
  cat(sprintf("f %.2f | %-38s R %5.2f g %5.3f | loci %3d known %2d fold %6.3f P %9.3g | %s\n",
              res$retention[i], substr(res$cohort[i], 1, 38), res$R[i], res$g[i],
              res$sig_loci[i], res$sig_known[i], res$fold[i], res$fisher_p[i],
              if (res$passes[i]) "PASS" else "fail"))

## 反解：§2.2 合规口径要达到 g = 1（f = 0.50）需要多少病例
for (ctrl in c(436283, 1e9)) {
  target <- 4 * N_EFF_DISC
  need <- 1 / (4 / target - 1 / ctrl)
  cat(sprintf("cases needed for R=4 with %.0f controls: %.0f  (MVP has 5,364 -> %.1fx)\n",
              ctrl, need, need / 5364))
}
write.table(res, "134a_candidate_evaluation.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("wrote 134a_candidate_evaluation.tsv (", nrow(res), " rows)\n", sep = "")
