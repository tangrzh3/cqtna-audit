## Modules A-G. Each is usable on its own; cqtna_audit() runs them together.
##
## Everything is computed by independent locus wherever a count carries an
## argument, because a significant locus does not name a gene.

#' A. Locus attribution against the outcome's own known loci
#'
#' Fold enrichment of FDR-significant loci on previously reported loci for this
#' outcome, counted by independent locus, with a one-sided Fisher exact test.
#'
#' Run it a second time with a different disease's list as `known` to get the
#' mismatched negative control. Without that control you cannot distinguish
#' attribution specific to your outcome's genetics from loci that are dense in
#' every disease.
#'
#' @param mr a `cqtna_mr` object from [as_cqtna_mr()].
#' @param known a `cqtna_known` object from [as_cqtna_known()].
#' @param known_kb window, in kb, within which a variant counts as landing on a
#'   known locus.
#' @param fdr FDR threshold.
#' @param label what this reference list is, used in the report.
#' @param known_from how a SIGNIFICANT locus inherits its known/novel status.
#'   The background is always taken over every record at a locus.
#'   \describe{
#'     \item{`"significant_records"`}{(default, and what the source study
#'       published) a significant locus is known if one of ITS significant
#'       records is near a known lead SNP. Numerator and denominator therefore
#'       use different record sets, which is a real inconsistency -- but it is
#'       the conservative direction, and it is the only one of the three that
#'       keeps the mismatched-list control clean on a dense exposure resource.}
#'     \item{`"any_record"`}{consistent with the background, but on a dense
#'       resource single-linkage chains distinct regions into blocks of tens of
#'       megabases, and a block that wide contains a known lead SNP for almost
#'       any disease. In the source study this raised the whole-blood cell from
#'       4.44- to 4.88-fold AND raised its mismatched control from 1.30-fold
#'       (P = 0.41) to 3.90-fold (P = 1.7e-4) -- a negative control that fails,
#'       which voids the cell under that study's own criterion.}
#'     \item{`"lead_variant"`}{consistent and immune to chaining, since each
#'       locus is represented by its most significant record. Raises the folds
#'       and leaves the mismatched control borderline (P = 0.096 on the same
#'       cell).}
#'   }
#'   Run all three and read the mismatched control before choosing. See
#'   [cqtna_locus_spans()] for whether chaining is in play at all.
#' @return a list with the counts, `fold`, `fisher_p_one_sided`, the convention
#'   used, and the significant genes split by whether their locus was known.
#' @export
#' @examples
#' mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38")
#' kn <- as_cqtna_known(cqtna_demo("known"), build = "GRCh38")
#' cqtna_attribution(mr, kn)$fold
cqtna_attribution <- function(mr, known, known_kb = 1000, fdr = 0.05,
                              label = "known-locus list",
                              known_from = c("significant_records", "any_record",
                                             "lead_variant")) {
  known_from <- match.arg(known_from)
  cq_check_build(mr, known)
  cq_validate_window(known_kb, "known_kb")
  cq_validate_fdr(fdr)
  idx <- cq_known_index(known)
  lk <- cq_locus_known(mr$locus, mr$chr, mr$pos, mr$fdr, idx, known_kb,
                       known_from, fdr)
  sel <- mr$fdr < fdr

  bg <- lk$by_locus
  sg <- lk$sig_status
  BT <- length(bg); BK <- sum(bg)
  ST <- length(sg); SK <- sum(sg)
  fold <- if (ST > 0 && BK > 0) (SK / ST) / (BK / BT) else NA_real_
  p <- if (ST > 0) cq_fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK))
       else NA_real_

  # 基因标签由**位点**状态给出，同一位点内所有基因必然一致
  st <- lk$sig_by_record
  structure(list(
    reference = label, known_from = known_from,
    background_known = as.integer(BK), background_loci = as.integer(BT),
    background_pct = if (BT) 100 * BK / BT else NA_real_,
    significant_loci = as.integer(ST), significant_known = as.integer(SK),
    pct_known = if (ST) 100 * SK / ST else NA_real_,
    fold = fold, fisher_p_one_sided = p,
    known_genes = sort(unique(mr$gene[sel & !is.na(st) & st])),
    novel_genes = sort(unique(mr$gene[sel & !is.na(st) & !st])),
    locus_known = bg),
    class = "cqtna_attribution")
}

#' B. The same list under four inference units
#'
#' The FDR-significant list recomputed over record, variant, gene and
#' independent locus, collapsing with Simes. A shortest list is not evidence
#' that FDR is controlled under dependence; the point is to state which unit the
#' conclusions are in.
#'
#' @inheritParams cqtna_attribution
#' @return a data frame, one row per unit.
#' @export
cqtna_unit_sweep <- function(mr, fdr = 0.05) {
  out <- data.frame(
    unit = "record", n_tests = nrow(mr),
    n_significant = sum(mr$fdr < fdr),
    n_genes = length(unique(mr$gene[mr$fdr < fdr])),
    n_loci = length(unique(mr$locus[mr$fdr < fdr])),
    stringsAsFactors = FALSE)
  keys <- list(variant = c("chr", "pos"), gene = "gene", locus = "locus")
  for (unit in names(keys)) {
    k <- keys[[unit]]
    grp <- interaction(mr[k], drop = TRUE, sep = "\r")
    pv <- tapply(mr$p, grp, cq_simes)
    q <- cq_bh(as.numeric(pv))
    hit <- names(pv)[q < fdr]
    sub <- mr[grp %in% hit, , drop = FALSE]
    out <- rbind(out, data.frame(
      unit = unit, n_tests = length(pv), n_significant = length(hit),
      n_genes = length(unique(sub$gene)), n_loci = length(unique(sub$locus)),
      stringsAsFactors = FALSE))
  }
  out
}

#' C. List stability against a second outcome GWAS
#'
#' Locus identity is recomputed on the UNION of both tables' coordinates, so the
#' two share one partition. Per-table locus ids are not comparable: on the demo
#' data both tables carry an internal id 206, one at 16:87.7 Mb and one at
#' 16:89.7 Mb, and a naive set intersection called them the same locus.
#'
#' @param mr,mr2 two `cqtna_mr` objects from the same exposure data and
#'   different outcome GWAS.
#' @param fdr FDR threshold.
#' @param locus_kb window used for the joint clustering. Defaults to the window
#'   both tables were built with, and errors if they disagree rather than
#'   silently picking one.
#' @return a list with Jaccard indices at gene and locus level, what moved, the
#'   shared locus keys, and a `coverage` block separating "lost" from "never
#'   tested in the second outcome".
#' @export
cqtna_stability <- function(mr, mr2, fdr = 0.05, locus_kb = NULL) {
  cq_validate_fdr(fdr)
  ba <- attr(mr, "build"); bb <- attr(mr2, "build")
  if (!identical(ba, bb))
    stop("the two MR tables are on different genome builds (", ba, " and ", bb,
         "); loci cannot be compared.", call. = FALSE)
  cq_check_chr_style(mr$chr, mr2$chr, "mr", "mr2")
  ka <- attr(mr, "locus_kb"); kb2 <- attr(mr2, "locus_kb")
  if (is.null(locus_kb)) {
    if (!identical(ka, kb2))
      stop("the two MR tables were clustered with different locus_kb (", ka,
           " and ", kb2, "). Rebuild both with the same window, or pass ",
           "`locus_kb` to re-cluster.", call. = FALSE)
    locus_kb <- ka
  }
  cq_validate_window(locus_kb, "locus_kb")

  # ⚠ 位点身份必须在**两表坐标的并集**上统一算出来。
  # 各自表内的位点编号不可比：在 demo 数据上，两表都有一个内部 id 206，
  # 一个是 16:87.7 Mb、另一个是 16:89.7 Mb，相距 2 Mb，
  # 而 intersect() 会把它们当成同一个位点。
  n1 <- nrow(mr); n2 <- nrow(mr2)
  chr <- c(mr$chr, mr2$chr); pos <- c(mr$pos, mr2$pos)
  uni <- cq_assign_loci(chr, pos, locus_kb)
  l1 <- uni[seq_len(n1)]; l2 <- uni[n1 + seq_len(n2)]

  s1 <- mr$fdr < fdr; s2 <- mr2$fdr < fdr
  a <- unique(mr$gene[s1]);  b <- unique(mr2$gene[s2])
  la <- unique(l1[s1]);      lb <- unique(l2[s2])
  jac <- function(x, y) if (length(union(x, y)))
    length(intersect(x, y)) / length(union(x, y)) else NA_real_

  # 记录覆盖不一致会让"丢失"混入"从未测过"，所以单独报出来
  v1 <- unique(paste(mr$chr, mr$pos, mr$gene))
  v2 <- unique(paste(mr2$chr, mr2$pos, mr2$gene))
  list(genes_outcome1 = length(a), genes_outcome2 = length(b),
       genes_shared = length(intersect(a, b)),
       jaccard_gene = jac(a, b),
       loci_outcome1 = length(la), loci_outcome2 = length(lb),
       loci_shared = length(intersect(la, lb)),
       jaccard_locus = jac(la, lb),
       lost = sort(setdiff(a, b)), gained = sort(setdiff(b, a)),
       locus_kb = locus_kb,
       shared_locus_keys = sort(intersect(la, lb)),
       coverage = list(
         records_outcome1 = n1, records_outcome2 = n2,
         gene_variant_pairs_shared = length(intersect(v1, v2)),
         gene_variant_pairs_only1 = length(setdiff(v1, v2)),
         gene_variant_pairs_only2 = length(setdiff(v2, v1)),
         genes_never_tested_in_2 = sort(setdiff(a, unique(mr2$gene))),
         genes_never_tested_in_1 = sort(setdiff(b, unique(mr$gene)))))
}

#' D. Instrument attrition at three levels
#'
#' Instrumentable, analysable against this outcome, nominally associated. These
#' are three different statements and collapsing them attributes an outcome-side
#' limitation to the biology of the pathway.
#'
#' @param instruments data frame with `gene`, `instrumentable`, `analysable`,
#'   `associated` (0/1 or logical).
#' @return a list of the three counts.
#' @export
cqtna_ladder <- function(instruments) {
  d <- as.data.frame(instruments, stringsAsFactors = FALSE)
  need <- c("instrumentable", "analysable", "associated")
  missing <- setdiff(need, names(d))
  if (length(missing))
    stop("instruments is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  f <- lapply(need, function(k) cq_as_flag(d[[k]], k))
  names(f) <- need

  # 三级是嵌套的：关联 => 可分析 => 可工具化。违反说明输入表本身错了，
  # 而不是"发现了什么"，所以报错而不是照算。
  bad <- (f$associated & !f$analysable) | (f$analysable & !f$instrumentable)
  if (any(bad)) {
    g <- if ("gene" %in% names(d)) as.character(d$gene)[which(bad)] else which(bad)
    stop("the three levels must be nested (associated => analysable => ",
         "instrumentable).
  violated by: ",
         paste(utils::head(g, 5), collapse = ", "),
         if (sum(bad) > 5) sprintf(" ... and %d more", sum(bad) - 5) else "",
         call. = FALSE)
  }
  list(pathway_genes = nrow(d),
       instrumentable = sum(f$instrumentable),
       analysable = sum(f$analysable),
       associated = sum(f$associated),
       note = paste("report all three; collapsing them attributes an",
                    "outcome-side limit to the biology"))
}

#' E. Compartment attribution
#'
#' For each nominated gene, the ratio of its highest expression in any other
#' cell type to its expression in the cell type the instrument came from. A
#' large ratio means tissue-level data cannot validate a target-cell-specific
#' mechanism.
#'
#' @param expression data frame with `gene`, `cell_type`, `mean_expression`.
#' @param genes character vector of nominated genes.
#' @param target_cell_type the cell type the exposure was measured in.
#' @return a data frame, one row per gene that could be evaluated.
#' @export
cqtna_compartment <- function(expression, genes, target_cell_type) {
  x <- as.data.frame(expression, stringsAsFactors = FALSE)
  need <- c("gene", "cell_type", "mean_expression")
  missing <- setdiff(need, names(x))
  if (length(missing))
    stop("expression is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  rows <- list()
  for (g in genes) {
    s <- x[x$gene == g, , drop = FALSE]
    if (!nrow(s) || is.null(target_cell_type) ||
        !target_cell_type %in% s$cell_type) next
    v <- stats::setNames(as.numeric(s$mean_expression), s$cell_type)
    ref <- v[[target_cell_type]]
    other <- v[names(v) != target_cell_type]
    if (!length(other)) next
    top <- names(other)[which.max(other)]
    # 目标细胞表达为 0 而别处 > 0，比值是 Inf，不是"缺失"。
    # 早先版本在 ref == 0 时返回 NA，把"完全不在目标细胞里表达"
    # 和"算不出来"混成了同一件事。
    ratio <- if (is.na(ref)) NA_real_
             else if (ref == 0 && other[[top]] > 0) Inf
             else if (ref == 0) NA_real_
             else other[[top]] / ref
    flag <- if (is.na(ratio)) "not evaluable"
            else if (is.infinite(ratio)) "absent from target cell type"
            else if (ratio > 2) "dominated by another compartment"
            else "comparable across compartments"
    rows[[length(rows) + 1L]] <- data.frame(
      gene = g, target_cell_type = target_cell_type,
      target_expression = ref, highest_other = top,
      highest_other_expression = other[[top]],
      ratio_other_over_target = ratio, compartment_flag = flag,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

#' F. Distance between eQTL and GWAS peaks
#'
#' @param peaks data frame with `gene`, `eqtl_pos`, `gwas_pos`.
#' @return the same rows with `distance_bp`.
#' @export
cqtna_peak_distance <- function(peaks) {
  x <- as.data.frame(peaks, stringsAsFactors = FALSE)
  need <- c("gene", "eqtl_pos", "gwas_pos")
  missing <- setdiff(need, names(x))
  if (length(missing))
    stop("peaks is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  x$distance_bp <- abs(as.numeric(x$eqtl_pos) - as.numeric(x$gwas_pos))
  x[, c("gene", "eqtl_pos", "gwas_pos", "distance_bp")]
}

#' Locus span diagnostic -- is single-linkage chaining doing the work?
#'
#' Single-linkage clustering joins two variants within `locus_kb` and then keeps
#' going, so in a dense resource a whole chromosome arm can chain into one
#' "locus". When that happens, asking whether ANY record at the locus is near a
#' known lead SNP stops being a question about the signal and becomes a question
#' about how wide the block grew.
#'
#' This was not hypothetical. In the study this tool came from, the CD4 time
#' course gave loci with a median span of 62 kb and none above 5 Mb, while the
#' whole-blood resource on the same outcome gave a median of 1,400 kb, thirty
#' loci above 10 Mb, and one significant locus spanning 30.8 Mb across 588
#' records. Locus-level counting means something very different in those two
#' cells, and nothing in the fold enrichment says so.
#'
#' @inheritParams cqtna_attribution
#' @return a list with the span quantiles, the number of loci far wider than the
#'   clustering window, and the spans of the significant loci.
#' @export
#' @examples
#' mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38")
#' cqtna_locus_spans(mr)$span_quantiles_kb
cqtna_locus_spans <- function(mr, fdr = 0.05) {
  cq_validate_fdr(fdr)
  kb <- as.numeric(attr(mr, "locus_kb"))
  span <- tapply(mr$pos, mr$locus, function(v) (max(v) - min(v)) / 1000)
  n_rec <- tapply(mr$pos, mr$locus, length)
  sig <- unique(as.character(mr$locus[mr$fdr < fdr]))
  wide <- span > 5 * kb
  out <- list(
    locus_kb = kb, n_loci = length(span),
    span_quantiles_kb = stats::quantile(span, c(.5, .9, .99, 1)),
    n_wider_than_5x_window = sum(wide),
    n_wider_than_10x_window = sum(span > 10 * kb),
    significant_span_kb = sort(span[sig]),
    significant_records = sort(n_rec[sig]),
    widest_locus = names(span)[which.max(span)])
  if (any(wide[sig]))
    warning(sum(wide[sig]), " significant locus/loci span more than ", 5 * kb,
            " kb, up to ", formatC(max(span[sig]), format = "d", big.mark = ","),
            " kb. Single-linkage chaining has merged distinct regions, so ",
            "locus-level counts on this resource are counting blocks, not loci. ",
            "Read cqtna_window_sweep()'s locus sweep before quoting a fold.",
            call. = FALSE)
  out
}

CQ_SWEEP_KB <- c(100, 250, 500, 1000)

#' G. Window sensitivity, and the distances behind the binary flag
#'
#' Two different conventions get called "the 1 Mb window" and they carry
#' different weight, so they are swept separately: `known_kb` decides what counts
#' as landing on a known locus and moves the fold directly; `locus_kb` decides
#' how independent loci are defined and moves the denominator.
#'
#' A threshold chosen to flatter a result weakens when tightened. If the fold
#' *rises* as the window narrows, the value you reported is the conservative one.
#'
#' Distances are reported on the same footing as the binary flag. A locus counts
#' as known when ANY of its records lies inside the window, so the headline
#' distance is the minimum over all records at that locus. The minimum over the
#' locus's significant records is reported alongside as
#' `significant_distances_sig_records_bp`: when the two differ, the locus is
#' near a known lead SNP because of a record that is not itself significant, and
#' that is worth seeing rather than hiding behind either convention.
#'
#' @inheritParams cqtna_attribution
#' @param locus_kb the locus-definition window held fixed while `known_kb` is
#'   swept.
#' @param sweep_kb windows to sweep, in kb.
#' @return a list with `known_sweep`, `locus_sweep`, the significant-locus
#'   distances, and `threshold_neighbourhood`: the nearest significant-locus
#'   distance below and above the chosen threshold, with their ratio and its
#'   log10. No verdict is attached. An earlier version printed "falls in a gap"
#'   when the ratio exceeded 4, a cutoff with no basis beyond looking reasonable
#'   on one dataset. Read the two neighbours and decide.
#' @export
cqtna_window_sweep <- function(mr, known, known_kb = 1000, locus_kb = 1000,
                               fdr = 0.05, sweep_kb = CQ_SWEEP_KB,
                               known_from = c("significant_records", "any_record",
                                              "lead_variant")) {
  known_from <- match.arg(known_from)
  cq_check_build(mr, known)
  cq_validate_window(known_kb, "known_kb")
  cq_validate_window(locus_kb, "locus_kb")
  cq_validate_window(sweep_kb, "sweep_kb")
  cq_validate_fdr(fdr)
  idx <- cq_known_index(known)
  dist_bp <- cq_nearest_bp(idx, mr$chr, mr$pos)
  sel <- mr$fdr < fdr

  # 与模块 A 同一口径：位点级 known 状态，分子分母都继承它
  enrich <- function(loci, kb) {
    rec <- dist_bp <= kb * 1000
    bg <- tapply(rec, loci, any)
    sg <- if (known_from == "any_record") bg[unique(as.character(loci[sel]))]
          else if (known_from == "lead_variant") {
            ord <- order(mr$fdr); keep <- !duplicated(loci[ord])
            r <- stats::setNames(rec[ord][keep], as.character(loci[ord][keep]))
            bg <- r[sort(names(r))]
            r[unique(as.character(loci[sel]))]
          } else tapply(rec[sel], loci[sel], any)
    BT <- length(bg); BK <- sum(bg); ST <- length(sg); SK <- sum(sg)
    fold <- if (ST > 0 && BK > 0) (SK / ST) / (BK / BT) else NA_real_
    data.frame(background_loci = BT, background_known = BK,
               significant_loci = ST, significant_known = SK,
               fold = fold,
               fisher_p_one_sided = if (ST > 0)
                 cq_fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK))
                 else NA_real_)
  }

  known_sweep <- do.call(rbind, lapply(sweep_kb, function(kb)
    cbind(known_window_kb = kb, locus_window_kb = locus_kb,
          enrich(mr$locus, kb))))
  locus_sweep <- do.call(rbind, lapply(sweep_kb, function(kb)
    cbind(locus_window_kb = kb, known_window_kb = known_kb,
          enrich(cq_assign_loci(mr$chr, mr$pos, kb), known_kb))))

  # ⚠ 距离必须与二分类同源。位点的 known 状态由**该位点的全部记录**决定
  # （见 cq_locus_known），所以主报的距离也取全部记录的最小值。
  # 另报只用显著记录算的距离：一个位点可能因为某条不显著的记录才靠近
  # 已知 lead SNP，那时两个数会分开，而这正是应该看见的事。
  bg_all <- tapply(dist_bp, mr$locus, min)
  sig_keys <- unique(as.character(mr$locus[sel]))
  sig_d <- sort(bg_all[sig_keys])
  sig_d <- sig_d[is.finite(sig_d)]
  sig_rec <- sort(tapply(dist_bp[sel], mr$locus[sel], min))
  sig_rec <- sig_rec[is.finite(sig_rec)]
  bg_d <- bg_all[is.finite(bg_all)]

  # 只报阈值上下最近的两个距离，不下二元判语。
  # 早先版本用 nearest_above >= 4 * nearest_below 自动印"落在空隙里"，
  # 那个 4 倍没有任何理论依据，却会被当成结论读。
  nb <- sig_d[sig_d <= known_kb * 1000]
  na_ <- sig_d[sig_d > known_kb * 1000]
  neighbourhood <- list(
    threshold_bp = known_kb * 1000,
    nearest_below_bp = if (length(nb)) as.numeric(nb[length(nb)]) else NA_real_,
    nearest_above_bp = if (length(na_)) as.numeric(na_[1]) else NA_real_)
  neighbourhood$ratio_above_over_below <-
    if (all(is.finite(unlist(neighbourhood[c("nearest_below_bp", "nearest_above_bp")]))) &&
        neighbourhood$nearest_below_bp > 0)
      neighbourhood$nearest_above_bp / neighbourhood$nearest_below_bp else NA_real_
  neighbourhood$log10_gap <-
    if (is.finite(neighbourhood$ratio_above_over_below))
      log10(neighbourhood$ratio_above_over_below) else NA_real_

  list(known_sweep = known_sweep, locus_sweep = locus_sweep,
       significant_distances_bp = as.numeric(sig_d),
       significant_distances_sig_records_bp = as.numeric(sig_rec),
       significant_median_bp = if (length(sig_d)) stats::median(sig_d) else NA_real_,
       background_median_bp = if (length(bg_d)) stats::median(bg_d) else NA_real_,
       threshold_neighbourhood = neighbourhood)
}
