## step138 -- S41：机制可迁移性检验（六个疾病，各打自己的名单）。
##
## 规则见 manuscript/PREREG_transport_grid.md（S41），冻结于提交 cf59123 / bd5cdfb，
## 均先于本脚本与任何自身名单结果。本脚本不新增任何判据。
##
## 主终点从未计算过：36e 把每个癌种都打在 Landi 的 melanoma 名单上
## （step30e 文档原文），因此它是错配对照表，不是迁移检验。
##
## 输出：138a_transport_main.tsv    每行 × 每份名单的富集
##       138b_transport_verdict.tsv 每行一个判定
##       138c_console.log

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
MIN_SIG_LOCI <- 8            # S41 §4 信息闸门（与 S40 §4 同一门槛）
MIN_LEAD     <- 30           # S41 §5.1 规则 5
MIN_BUILD    <- 0.50         # S41 §5.1 规则 6
N_PERM <- 10000; SEED <- 1; TOL <- 1.00

fix_gene <- function(g) {
  g <- as.character(g); bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad)); g
}
std <- function(gene, chr, pos, p) {
  d <- data.frame(record_id = seq_along(pos), gene = fix_gene(gene),
                  chr = as.character(chr), pos = as.numeric(pos),
                  p = as.numeric(p), stringsAsFactors = FALSE)
  d[is.finite(d$p) & is.finite(d$pos) & d$pos > 0, , drop = FALSE]
}

## ---- 名单（S41 §5：先过闸门，再进面板） -----------------------------------
audit <- read.delim("137a_list_rebuild_audit.tsv", stringsAsFactors = FALSE)
list_files <- c(melanoma = "landi2020_known_loci_grch38.csv",
                HCC = "84a_hcc_known_loci_grch38.csv",
                colorectal = "known_loci_colorectal_grch38.csv",
                prostate = "known_loci_prostate_grch38.csv",
                breast = "known_loci_breast_grch38.csv",
                lung = "known_loci_lung_grch38.csv",
                pancreas = "known_loci_pancreas_grch38.csv")
K <- list(); list_n <- integer(0); list_gate <- character(0)
for (nm in names(list_files)) {
  f <- list_files[[nm]]
  d <- if (file.exists(f)) tryCatch(read.csv(f, stringsAsFactors = FALSE),
                                    error = function(e) data.frame()) else data.frame()
  if (nrow(d)) d <- d[!is.na(suppressWarnings(as.numeric(d$pos))), , drop = FALSE]
  n <- nrow(d); list_n[nm] <- n
  br <- audit$build_rate[match(nm, audit$disease)]
  g <- if (n < MIN_LEAD) sprintf("inconclusive (list): %d below %d", n, MIN_LEAD)
       else if (!is.na(br) && br < MIN_BUILD)
         sprintf("inconclusive (list): build rate %.1f%%", 100 * br)
       else "usable"
  list_gate[nm] <- g
  if (g == "usable") {
    d$source <- nm
    K[[nm]] <- as_cqtna_known(d[, c("chr", "pos", "source")], build = "GRCh38")
  }
}
cat("known-locus lists (S41 section 5 gates):\n")
for (nm in names(list_files))
  cat(sprintf("  %-11s %5d  %s\n", nm, list_n[nm], list_gate[nm]))

## ---- 六行（S41 §3，事前固定，不得增删） ------------------------------------
mel <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
mel$chr <- sub(":.*", "", mel$SNP); mel$pos <- as.numeric(sub(".*:", "", mel$SNP))
cc <- read.delim("36a_crosscancer_MR_all.tsv", stringsAsFactors = FALSE)
cc$chr <- sub(":.*", "", cc$SNP); cc$pos <- as.numeric(sub(".*:", "", cc$SNP))

rows <- list(list(nm = "Melanoma", own = "melanoma",
                  d = std(mel$SYMBOL, mel$chr, mel$pos, mel$pval)))
for (cn in c("Lung", "Colorectal", "Breast", "Prostate", "Pancreas")) {
  s <- cc[cc$cancer == cn, , drop = FALSE]
  rows[[length(rows) + 1]] <- list(nm = cn, own = tolower(cn),
                                   d = std(s$symbol, s$chr, s$pos, s$pval))
}

main <- list(); verd <- list()
for (r in rows) {
  mr <- as_cqtna_mr(r$d, locus_kb = LOCUS_KB, build = "GRCh38",
                    locus_method = "fixed_centre")
  cat(sprintf("\n  %s  (%d records, %d significant records)\n",
              r$nm, nrow(mr), sum(mr$fdr < FDR)))
  cat("    list         relation           bg known    sig known     fold        P\n")
  own_row <- NULL
  for (nm in names(K)) {
    a <- cqtna_attribution(mr, K[[nm]], known_kb = KNOWN_KB, fdr = FDR,
                           known_from = "any_record")
    rel <- if (nm == r$own) "own outcome" else "unrelated control"
    main[[length(main) + 1]] <- data.frame(
      row = r$nm, list_name = nm, relation = rel,
      bg_loci = a$background_loci, bg_known = a$background_known,
      sig_loci = a$significant_loci, sig_known = a$significant_known,
      fold = a$fold, fisher_p = a$fisher_p_one_sided, stringsAsFactors = FALSE)
    if (rel == "own outcome") own_row <- main[[length(main)]]
    cat(sprintf("    %-11s  %-17s %4d/%-6d %3d/%-7d %7s %8s\n", nm, rel,
                a$background_known, a$background_loci,
                a$significant_known, a$significant_loci,
                ifelse(is.finite(a$fold), sprintf("%.2f", a$fold), "NA"),
                ifelse(is.finite(a$fisher_p_one_sided),
                       sprintf("%.3g", a$fisher_p_one_sided), "NA")))
  }

  ## 判定（S41 §7）
  perm_p <- NA_real_; perm_fold <- NA_real_
  if (is.null(own_row)) {
    v <- sprintf("inconclusive (list): own list %s failed its gate", r$own)
  } else if (own_row$sig_loci < MIN_SIG_LOCI) {
    v <- sprintf("inconclusive (power): %d below %d significant loci",
                 own_row$sig_loci, MIN_SIG_LOCI)
  } else {
    pc <- cqtna_permutation_control(mr, K[[r$own]], known_kb = KNOWN_KB,
                                    fdr = FDR, known_from = "any_record",
                                    tolerance = TOL, n_perm = N_PERM, seed = SEED)
    perm_p <- if (is.null(pc$empirical_p)) NA_real_ else pc$empirical_p
    perm_fold <- if (is.null(pc$fold_vs_null)) NA_real_ else pc$fold_vs_null
    edf <- if (is.null(pc$effective_draw_fraction)) NA_real_ else pc$effective_draw_fraction
    ctl <- do.call(rbind, main)
    ctl <- ctl[ctl$row == r$nm & ctl$relation == "unrelated control", ]
    ctl$p_holm <- stats::p.adjust(ctl$fisher_p, "holm")
    hit <- ctl[is.finite(ctl$p_holm) & ctl$p_holm < 0.05 &
               is.finite(ctl$fold) & ctl$fold > 1, , drop = FALSE]
    fmax <- suppressWarnings(max(ctl$fold[is.finite(ctl$fold)]))
    v <- if (!is.finite(perm_p))
           "inconclusive (permutation)"
         else if (nrow(hit) > 0 || !(is.finite(own_row$fold) && own_row$fold >= fmax))
           "void"
         else if (own_row$fold > 1 && own_row$fisher_p < 0.05 && perm_fold > 1)
           "supportive"
         else "negative"
    cat(sprintf("    permutation: empirical P %s, fold vs null %s, effective draws %s\n",
                ifelse(is.finite(perm_p), sprintf("%.4g", perm_p), "NA"),
                ifelse(is.finite(perm_fold), sprintf("%.2f", perm_fold), "NA"),
                ifelse(is.finite(edf), sprintf("%.1f%%", 100 * edf), "NA")))
    if (!is.finite(perm_p) && !is.null(pc$failed_because))
      cat(sprintf("    permutation failed: %s\n", substr(pc$failed_because, 1, 90)))
  }
  verd[[length(verd) + 1]] <- data.frame(
    row = r$nm, own_list = r$own, n_records = nrow(mr),
    sig_loci = if (is.null(own_row)) NA_integer_ else own_row$sig_loci,
    sig_known = if (is.null(own_row)) NA_integer_ else own_row$sig_known,
    bg_known = if (is.null(own_row)) NA_integer_ else own_row$bg_known,
    bg_loci = if (is.null(own_row)) NA_integer_ else own_row$bg_loci,
    fold = if (is.null(own_row)) NA_real_ else own_row$fold,
    fisher_p = if (is.null(own_row)) NA_real_ else own_row$fisher_p,
    perm_p = perm_p, perm_fold = perm_fold,
    verdict = v, stringsAsFactors = FALSE)
  cat(sprintf("    VERDICT: %s\n", v))
}

M <- do.call(rbind, main); V <- do.call(rbind, verd)
## S41 §7：自身名单 P 跨六行 Holm 校正
V$fisher_p_holm <- stats::p.adjust(V$fisher_p, "holm")
write.table(M, "138a_transport_main.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(V, "138b_transport_verdict.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n", strrep("=", 104), "\n", sep = "")
cat("S41 transport grid -- own-list result per disease\n")
cat(sprintf("%-11s %9s %9s %9s %10s %10s  %s\n",
            "outcome", "sig loci", "own fold", "P", "P (Holm)", "perm P", "verdict"))
for (i in seq_len(nrow(V)))
  cat(sprintf("%-11s %9s %9s %9s %10s %10s  %s\n", V$row[i],
              ifelse(is.na(V$sig_loci[i]), "-", V$sig_loci[i]),
              ifelse(is.finite(V$fold[i]), sprintf("%.2f", V$fold[i]), "-"),
              ifelse(is.finite(V$fisher_p[i]), sprintf("%.3g", V$fisher_p[i]), "-"),
              ifelse(is.finite(V$fisher_p_holm[i]), sprintf("%.3g", V$fisher_p_holm[i]), "-"),
              ifelse(is.finite(V$perm_p[i]), sprintf("%.4g", V$perm_p[i]), "-"),
              V$verdict[i]))
cat("\nwrote 138a / 138b\n")
