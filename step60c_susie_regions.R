# Step 60c  分析 B 主检验：结局侧在目标区域有几个独立信号
#
# 设计见 manuscript/PREREG_power_trajectory.md §3
# 逻辑：两轮之间暴露侧数据完全相同，故 ④b 中 H3/H4 的变化必来自结局侧。
#       若某区域在两个功效下结局侧都只有一个可信信号，则"结局侧多信号未被建模"
#       不能解释该区域的移动。
#
# ⚠ 硬边界（不因任何资源而改变）：eQTL 侧 n=85-100 且无样本内 LD，
#   其 susie 结果不可靠，故本脚本**只做结局侧**。
#
# ⚠ 流程阳性对照：MC1R 区。FinnGen 官方用**样本内 LD** 给出 3 个高纯度 credible set
#   （cs1/cs2/cs3, log10BF 36.6/45.6/19.5）。我们用 1000G 代理 LD 必须复现 >=2 个。
#   不通过 = 代理 LD 不够用 = 整段分析记为无信息（不得当作"未推翻结论"）。

suppressPackageStartupMessages({library(susieR); library(data.table)})
MR  <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR"; REG <- file.path(MR, "regions")
N_FG   <- 384502    # FinnGen R12: 5,753 例 / 378,749 对照
N_META <- 801629    # meta: 12,530 / 789,099

regions <- c("MC1R", "PARP1", "ZFYVE19", "TPI1")
srcs    <- c(finngen = N_FG, meta = N_META)
res <- list()

for (rg in regions) {
  ldf <- file.path(REG, paste0(rg, "_ld.unphased.vcor1"))
  vrf <- paste0(ldf, ".vars")
  if (!file.exists(ldf)) { cat("[缺 LD]", rg, "\n"); next }
  vars <- fread(vrf, header = FALSE)$V1
  R <- as.matrix(fread(ldf, header = FALSE))
  colnames(R) <- rownames(R) <- vars
  cat(sprintf("\n===== %s : LD %d x %d =====\n", rg, nrow(R), ncol(R)))

  for (src in names(srcs)) {
    f <- file.path(REG, paste0(rg, "_", src, ".tsv"))
    if (!file.exists(f)) { cat("  [缺]", src, "\n"); next }
    d <- fread(f)
    d[, id  := paste0(chrom, ":", pos, ":", ref, ":", alt)]
    d[, idf := paste0(chrom, ":", pos, ":", alt, ":", ref)]   # 等位互换写法
    d <- d[!is.na(beta) & !is.na(se) & se > 0]

    # 与 LD 面板对齐；面板 ID 为 chr:pos:ref:alt（ref 为 REF 等位）
    d[, key := fifelse(id %in% vars, id, fifelse(idf %in% vars, idf, NA_character_))]
    d[, flip := !is.na(key) & key == idf]
    d <- d[!is.na(key)]
    d <- d[!duplicated(key)]
    if (nrow(d) < 50) { cat(sprintf("  %-8s 对齐后仅 %d 个变异，跳过\n", src, nrow(d))); next }

    z <- d$beta / d$se
    z[d$flip] <- -z[d$flip]          # 翻转到面板的 REF/ALT 定向
    Rs <- R[d$key, d$key, drop = FALSE]

    # LD 失配诊断（预注册要求随结果同行）
    s_est <- tryCatch(estimate_s_rss(z = z, R = Rs, n = srcs[[src]]),
                      error = function(e) NA_real_)

    fit <- tryCatch(
      susie_rss(z = z, R = Rs, n = srcs[[src]], L = 10,
                estimate_residual_variance = FALSE),
      error = function(e) NULL)
    if (is.null(fit)) { cat(sprintf("  %-8s susie 未收敛\n", src)); next }

    cs <- fit$sets$cs
    pur <- fit$sets$purity
    ncs <- if (is.null(cs)) 0L else length(cs)
    # 高纯度定义与 FinnGen 官方一致：min |r| >= 0.5
    hi <- if (ncs == 0) 0L else sum(pur$min.abs.corr >= 0.5, na.rm = TRUE)
    top <- d[which.min(pval)]
    cat(sprintf("  %-8s n=%d 变异 | minP=%.2e @ %s:%d | credible sets=%d (高纯度 %d) | estimate_s=%.3f\n",
                src, nrow(d), top$pval, top$chrom, top$pos, ncs, hi,
                ifelse(is.na(s_est), NA, s_est)))
    if (ncs > 0) for (k in seq_len(ncs)) {
      idx <- cs[[k]]
      cat(sprintf("       cs%d size=%d  min|r|=%.3f  top PIP var=%s\n",
                  k, length(idx), pur$min.abs.corr[k],
                  d$key[idx[which.max(fit$pip[idx])]]))
    }
    res[[length(res) + 1]] <- data.table(
      region = rg, source = src, n_var = nrow(d),
      min_pval = top$pval, top_pos = top$pos,
      n_cs = ncs, n_cs_highpurity = hi, estimate_s = s_est)
  }
}

if (length(res)) {
  out <- rbindlist(res)
  fwrite(out, file.path(MR, "60a_susie_outcome_signals.tsv"), sep = "\t")
  cat("\n写出 60a_susie_outcome_signals.tsv\n")
  print(out)

  cat("\n===== 阳性对照判定 =====\n")
  mc <- out[region == "MC1R"]
  if (nrow(mc)) {
    pass <- any(mc$n_cs_highpurity >= 2)
    cat(sprintf("MC1R 高纯度 credible set 数: %s -> %s\n",
                paste(mc$source, mc$n_cs_highpurity, sep = "=", collapse = ", "),
                ifelse(pass, "PASS（代理 LD 可用）",
                       "FAIL（代理 LD 不足，本段分析记为无信息）")))
  }
}
