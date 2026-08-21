## step135 -- §9 的可选 HCC 扩展，离它自己的八位点信息闸门有多远。
##
## §9 允许把独立 HCC 作为第二阶段单独预注册，并要求它"pass its own eight-locus
## information gate"。melanoma 那一格发现集实得 8 个显著位点，恰好卡在闸门上；
## HCC 那一格实得 2 个。本脚本用与 step133/134 相同的确定性投影（读法 B）
## 算出 HCC 要达到 8 个位点需要多大的 g，进而需要多大的 N_eff。
##
## 判据、分区、归属口径、参照名单全部沿用，未作任何改动。
## 输出：135a_hcc_gate_distance.tsv、135b_console.log

suppressMessages(library(cqtna))
arg <- commandArgs(trailingOnly = TRUE)
root <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
        if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else getwd()
setwd(root)

FDR <- 0.05; LOCUS_KB <- 1000; KNOWN_KB <- 1000
neff <- function(a, b) 4 / (1 / a + 1 / b)

fix_gene <- function(g) {
  g <- as.character(g); bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad)); g
}
std <- function(gene, chr, pos, p, bo, so) {
  d <- data.frame(record_id = seq_along(pos), gene = fix_gene(gene),
                  chr = as.character(chr), pos = as.numeric(pos),
                  p = as.numeric(p), bo = bo, so = so, stringsAsFactors = FALSE)
  d[is.finite(d$p) & is.finite(d$pos) & d$pos > 0, , drop = FALSE]
}

mel <- read.csv("landi2020_known_loci_grch38.csv"); mel$source <- "Landi2020"
hcc <- read.csv("84a_hcc_known_loci_grch38.csv");   hcc$source <- "step84"
K_mel <- as_cqtna_known(mel, build = "GRCh38")
K_hcc <- as_cqtna_known(hcc, build = "GRCh38")

cd4 <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
cd4$chr <- sub(":.*", "", cd4$SNP); cd4$pos <- as.numeric(sub(".*:", "", cd4$SNP))
hh <- read.delim("85a_HCC_high_annotated.tsv", stringsAsFactors = FALSE)

cells <- list(
  list(nm = "melanoma x Soskic_CD4", K = K_mel, n_eff = neff(12530, 789099),
       d = std(cd4$SYMBOL, cd4$chr, cd4$pos, cd4$pval,
               cd4$beta_outcome, cd4$se_outcome)),
  list(nm = "HCC_high x Soskic_CD4", K = K_hcc, n_eff = neff(3748, 1861536),
       d = std(hh$gene_id, hh$chr, hh$pos, hh$p_mr, NA, NA)))

## HCC 表里若无 outcome 侧 beta/se，就由 p 反解 |z|（符号不影响双侧 p，
## 而本投影只用 |z|）。melanoma 表两条路都有，用来核对反解是否等价。
z_from_p <- function(p) stats::qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
{
  zt <- abs(cells[[1]]$d$bo / cells[[1]]$d$so); zp <- z_from_p(cells[[1]]$d$p)
  ok <- is.finite(zt) & is.finite(zp) & cells[[1]]$d$p > 1e-290
  stopifnot(max(abs(zt[ok] - zp[ok])) < 1e-6)
  cat("z-from-p check on melanoma: OK (max diff",
      signif(max(abs(zt[ok] - zp[ok])), 3), ")\n")
}

run_cell <- function(cl, g) {
  d <- cl$d
  z <- z_from_p(d$p)
  d$p <- 2 * stats::pnorm(-abs(g * z))
  mrx <- as_cqtna_mr(d[, c("record_id", "gene", "chr", "pos", "p")],
                     locus_kb = LOCUS_KB, build = "GRCh38",
                     locus_method = "fixed_centre")
  ax <- cqtna_attribution(mrx, cl$K, known_kb = KNOWN_KB, fdr = FDR,
                          known_from = "any_record")
  c(ST = ax$significant_loci, SK = ax$significant_known,
    fold = ax$fold, p = ax$fisher_p_one_sided)
}

out <- list()
for (cl in cells) {
  v1 <- run_cell(cl, 1)
  cat(sprintf("\n%s  N_eff %.0f  |  g=1 reproduces: loci %d known %d fold %.3f P %.5g\n",
              cl$nm, cl$n_eff, v1[["ST"]], v1[["SK"]], v1[["fold"]], v1[["p"]]))
  for (g in c(1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4)) {
    v <- run_cell(cl, g)
    R <- g^2                      # f = 1.00 时 R = g^2；f = 0.50 时 R = 4g^2
    out[[length(out) + 1]] <- data.frame(
      cell = cl$nm, n_eff_disc = round(cl$n_eff), g = g,
      R_at_f1 = R, R_at_f0.5 = 4 * R,
      n_eff_needed_f0.5 = round(cl$n_eff * 4 * R),
      sig_loci = v[["ST"]], sig_known = v[["SK"]], fold = v[["fold"]],
      fisher_p = v[["p"]],
      meets_8 = v[["ST"]] >= 8,
      passes = v[["ST"]] >= 8 && is.finite(v[["fold"]]) && v[["fold"]] > 1 &&
               is.finite(v[["p"]]) && v[["p"]] < 0.05)
    cat(sprintf("  g %4.2f | loci %4d known %3d fold %7.3f P %9.3g | N_eff needed at f=0.50: %10.0f | %s\n",
                g, v[["ST"]], v[["SK"]], v[["fold"]], v[["p"]],
                cl$n_eff * 4 * R, if (out[[length(out)]]$passes) "PASS" else "fail"))
  }
}
res <- do.call(rbind, out)
write.table(res, "135a_hcc_gate_distance.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nwrote 135a_hcc_gate_distance.tsv (", nrow(res), " rows)\n", sep = "")
