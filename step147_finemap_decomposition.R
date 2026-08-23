## step147 -- 精细定位版的模块 H 分解：已达全基因组显著的位点与未达的，
##            在结局侧各有几个可分辨信号？
##
## 复用 step60c 的全部约定（1000G EUR 代理 LD、susie_rss、报告 estimate_s），
## 并复用 PREREG_power_trajectory.md §9.1 **已经收紧**的读法：
##
##   >= 2 个 credible set  → 不可采信（代理 LD 可能虚假分裂）
##   恰好 1 个且对 L 稳健  → 可采信，且保守
##   0 个                  → 该功效下无可分辨的结局侧信号
##
## 本脚本不放松任何一条，也不新增判据。
##
## ⚠ 已知的界限，先写在这里：亚阈值位点给出 0 个 credible set 是**预期之内**的，
## 因为精细定位本来就需要强信号。所以这不是独立证据，而是**与模块 H 的分解一致**。
## 不得读成"已证明亚阈值部分为空"。
##
## 输出：147a_finemap_decomposition.tsv、147b_console.log

suppressPackageStartupMessages({library(susieR); library(data.table)})

arg <- commandArgs(trailingOnly = TRUE)
MR <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
      if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else getwd()
REG <- file.path(MR, "regions")
setwd(MR)

N_FG <- 384502      # FinnGen R12: 5,753 / 378,749
N_META <- 801629    # meta: 12,530 / 789,099
srcs <- c(finngen = N_FG, meta = N_META)
SEED <- 1
SEEDS <- 1:5      # 种子敏感性：同一份输入换种子，credible set 数会不会变

## 八个位点 = 模块 H 的显著 bounded locus。GWS 标记来自 140b。
REGIONS <- list(
  list(rg = "MC1R",     gws = TRUE,  note = "16q24.3; two module-H loci fall here"),
  list(rg = "PARP1",    gws = TRUE,  note = ""),
  list(rg = "ZFYVE19",  gws = FALSE, note = ""),
  list(rg = "SMC2",     gws = FALSE, note = ""),
  list(rg = "KANSL1",   gws = FALSE, note = "17q21.31 inversion; LD is unusual here"),
  list(rg = "MDM4",     gws = FALSE, note = ""),
  list(rg = "KIAA0040", gws = FALSE, note = ""))

verdict <- function(ncs) {
  if (!is.finite(ncs)) "inconclusive"
  else if (ncs == 0) "0 - no distinguishable signal at this power"
  else if (ncs == 1) "1 - admissible (conservative)"
  else "not admissible - proxy LD may split spuriously"
}

out <- list()
for (R0 in REGIONS) {
  rg <- R0$rg
  ldf <- file.path(REG, paste0(rg, "_ld.unphased.vcor1"))
  vrf <- paste0(ldf, ".vars")
  if (!file.exists(ldf)) { cat("[no LD]", rg, "\n"); next }
  vars <- fread(vrf, header = FALSE)$V1
  R <- as.matrix(fread(ldf, header = FALSE))
  colnames(R) <- rownames(R) <- vars
  cat(sprintf("\n===== %s (%s) : LD %d x %d =====\n", rg,
              if (R0$gws) "genome-wide significant" else "below genome-wide",
              nrow(R), ncol(R)))

  for (src in names(srcs)) {
    f <- file.path(REG, paste0(rg, "_", src, ".tsv"))
    if (!file.exists(f)) { cat("  [no sumstats]", src, "\n"); next }
    d <- fread(f)
    ## step60a 的区域文件把标准误列叫 `se`，step146 保留了源文件的 `sebeta`。
    ## 统一到 `se`，而不是让两套文件走两条代码路径。
    if (!"se" %in% names(d) && "sebeta" %in% names(d))
      data.table::setnames(d, "sebeta", "se")
    stopifnot(all(c("chrom", "pos", "ref", "alt", "beta", "se", "pval") %in% names(d)))
    d[, id  := paste0(chrom, ":", pos, ":", ref, ":", alt)]
    d[, idf := paste0(chrom, ":", pos, ":", alt, ":", ref)]
    d <- d[!is.na(beta) & !is.na(se) & se > 0]
    d[, key := fifelse(id %in% vars, id, fifelse(idf %in% vars, idf, NA_character_))]
    d[, flip := !is.na(key) & key == idf]
    d <- d[!is.na(key)][!duplicated(key)]
    if (nrow(d) < 50) {
      cat(sprintf("  %-8s only %d aligned variants, skipped\n", src, nrow(d)))
      next
    }
    z <- d$beta / d$se
    z[d$flip] <- -z[d$flip]
    Rs <- R[d$key, d$key, drop = FALSE]

    s_est <- tryCatch(estimate_s_rss(z = z, R = Rs, n = srcs[[src]]),
                      error = function(e) NA_real_)
    ## ⚠ susie_get_cs 对大 credible set 的 purity 用随机抽样估计（n_purity = 100）。
    ## 一个含 1,851 个变异、purity 恰在阈值附近的 set，会因抽样而时进时出：
    ## KANSL1 在无种子的三次运行里给出 0 / 1 / 0 个 credible set。加种子使之可复现，
    ## 但**可复现不等于稳健**，所以下面还要跨种子测一次敏感性。
    cs_at <- function(L, seed = SEED) {
      set.seed(seed)
      fit <- tryCatch(susie_rss(z = z, R = Rs, n = srcs[[src]], L = L),
                      error = function(e) NULL)
      if (is.null(fit))
        return(c(n = NA_real_, hi = NA_real_, size = NA_real_,
                 span_kb = NA_real_, converged = NA_real_))
      conv <- if (is.null(fit$converged)) NA_real_ else as.numeric(fit$converged)
      cs <- tryCatch(susie_get_cs(fit, Xcorr = Rs), error = function(e) NULL)
      if (is.null(cs) || is.null(cs$cs) || !length(cs$cs))
        return(c(n = 0, hi = 0, size = NA_real_, span_kb = NA_real_,
                 converged = conv))
      hi <- if (!is.null(cs$purity))
        sum(cs$purity$min.abs.corr >= 0.5, na.rm = TRUE) else NA_real_
      ## 个数不够。一个含 1,851 个变异、跨 650 kb 的 credible set 与一个含 35 个、
      ## 跨 68 kb 的，在"1 个 credible set"这句话里看不出区别，而它们不是一回事：
      ## 前者是 17q21.31 倒位单体型整块在动，什么也没定位到。
      sz <- vapply(cs$cs, length, integer(1))
      sp <- vapply(cs$cs, function(ix) diff(range(d$pos[ix])), numeric(1))
      c(n = length(cs$cs), hi = hi, size = max(sz), span_kb = max(sp) / 1000,
        converged = conv)
    }
    a10 <- cs_at(10); a20 <- cs_at(20)
    ## 跨种子重跑 L=10，记录 credible set 数的取值范围
    ns <- vapply(SEEDS, function(sd) cs_at(10, sd)[["n"]], numeric(1))
    seed_stable <- length(unique(ns[is.finite(ns)])) <= 1
    stable <- is.finite(a10[["n"]]) && is.finite(a20[["n"]]) &&
              a10[["n"]] == a20[["n"]]
    v <- verdict(a10[["n"]])
    if (a10[["n"]] == 1 && !stable) v <- "1 at L=10 but not stable in L"
    if (!seed_stable) v <- paste0(v, "; NOT STABLE ACROSS SEEDS (",
                                  paste(ns, collapse = "/"), ")")

    cat(sprintf("  %-8s n=%d | minP=%.2e | cs L=10: %s (largest %s var / %s kb, conv %s), L=20: %s | s=%.3f | %s\n",
                src, nrow(d), min(d$pval, na.rm = TRUE), a10[["n"]],
                ifelse(is.finite(a10[["size"]]), a10[["size"]], "-"),
                ifelse(is.finite(a10[["span_kb"]]), round(a10[["span_kb"]]), "-"),
                ifelse(is.finite(a10[["converged"]]), a10[["converged"]] == 1, "?"),
                a20[["n"]], s_est, v))
    out[[length(out) + 1]] <- data.frame(
      region = rg, genome_wide_sig = R0$gws, source = src,
      n_var = nrow(d), min_pval = min(d$pval, na.rm = TRUE),
      n_cs_L10 = a10[["n"]], n_cs_highpurity_L10 = a10[["hi"]],
      largest_cs_size = a10[["size"]], largest_cs_span_kb = a10[["span_kb"]],
      converged_L10 = a10[["converged"]], converged_L20 = a20[["converged"]],
      n_cs_across_seeds = paste(ns, collapse = ","), seed_stable = seed_stable,
      n_cs_L20 = a20[["n"]], stable_in_L = stable,
      estimate_s = s_est, verdict = v, note = R0$note,
      stringsAsFactors = FALSE)
  }
}

res <- do.call(rbind, out)
write.table(res, "147a_finemap_decomposition.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("\n", strrep("=", 92), "\n", sep = "")
cat("outcome-side signals, split by whether the outcome GWAS already reached 5e-8\n\n")
for (g in c(TRUE, FALSE)) {
  cat(if (g) "genome-wide significant loci:\n" else "\nbelow genome-wide significance:\n")
  s <- res[res$genome_wide_sig == g & res$source == "meta", ]
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-9s minP %9.2e  cs %2s (largest %5s var /%5s kb, conv %-5s)  s %.3f  %s\n",
                s$region[i], s$min_pval[i], s$n_cs_L10[i],
                ifelse(is.finite(s$largest_cs_size[i]), s$largest_cs_size[i], "-"),
                ifelse(is.finite(s$largest_cs_span_kb[i]), round(s$largest_cs_span_kb[i]), "-"),
                ifelse(is.finite(s$converged_L10[i]), s$converged_L10[i] == 1, "?"),
                s$estimate_s[i], s$verdict[i]))
}
cat("\nwrote 147a\n")
