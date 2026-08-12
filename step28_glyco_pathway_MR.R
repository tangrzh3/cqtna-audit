## =====================================================================
## Step 28  Is TPI1 the only glycolytic gene with a causal signal, or is
##          there a pathway-level MR result?
##
## WHY THIS MATTERS
## Step 26 showed the FUNCTIONAL phenotype is pathway-wide (Post 20/22 genes
## higher in non-responders, 8 surviving FDR; PGAM1 the strongest at Pre).
## But MR named only TPI1 -- because TPI1 is where a genome-wide significant
## dynamic eQTL instrument exists. That asymmetry is currently written up as a
## limitation ("instrument availability, not biological importance, decides
## which gene an MR study names").
##
## If OTHER glycolytic genes also carry instruments and show concordant MR
## effects, the limitation converts into a strength: several independent loci
## in one pathway pointing the same way is far stronger evidence than one gene.
## If they do not, the limitation stands as written.
## =====================================================================

suppressPackageStartupMessages({library(data.table)})
setwd("D:/R_ex/MR")

GLYCO <- c("SLC2A1", "SLC2A3", "HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM",
           "PFKP", "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "PGK1",
           "PGAM1", "ENO1", "ENO2", "ENO3", "PKM", "PKLR", "LDHA", "LDHB",
           "PGM1", "PFKFB3", "PFKFB4", "SLC16A1", "SLC16A3")

ann <- fread("13_meta_locus_annotation.tsv")
cat("MR records:", nrow(ann), " | unique genes:", uniqueN(ann$SYMBOL), "\n")

g <- ann[SYMBOL %in% GLYCO]
cat("\n=== glycolytic genes with a strict-set instrument ===\n")
if (!nrow(g)) {
  cat("none\n")
} else {
  out <- g[, .(SYMBOL, exposure, cell_type, timepoint, SNP, rsid,
               pval_exposure, OR = round(OR, 3), OR_LCI = round(OR_LCI, 3),
               OR_UCI = round(OR_UCI, 3), pval, FDR = round(FDR, 4),
               category, near_locus, dist_kb)]
  setorder(out, pval)
  print(out, digits = 3)
  fwrite(out, "34a_glycolysis_MR_records.tsv", sep = "\t")

  cat("\n--- summary by gene ---\n")
  s <- g[, .(n_profiles = .N,
             best_p = min(pval), best_FDR = min(FDR),
             OR_at_best = OR[which.min(pval)],
             direction = ifelse(OR[which.min(pval)] > 1, "risk", "protective"),
             profiles = paste(unique(paste0(cell_type, "_", timepoint)), collapse = ", ")),
         by = SYMBOL][order(best_p)]
  print(s, digits = 3)
  fwrite(s, "34b_glycolysis_MR_by_gene.tsv", sep = "\t")

  cat("\nnominal p<0.05:", s[best_p < 0.05, .N], "of", nrow(s), "genes with instruments\n")
  cat("FDR<0.05      :", s[best_FDR < 0.05, .N], "\n")
  cat("risk-direction (OR>1) among nominal:", s[best_p < 0.05 & direction == "risk", .N],
      "of", s[best_p < 0.05, .N], "\n")
}

## which glycolytic genes have NO instrument at all
cat("\n=== glycolytic genes with no strict-set instrument ===\n")
cat(paste(setdiff(GLYCO, unique(ann$SYMBOL)), collapse = ", "), "\n")

## relaxed / sensitivity set
if (file.exists("18_sensitivity_meta.tsv")) {
  sens <- fread("18_sensitivity_meta.tsv")
  nm <- names(sens)
  key <- grep("exposure|outcome|method|^b$|^se$|^pval$|nsnp", nm, value = TRUE, ignore.case = TRUE)
  cat("\n=== sensitivity-set columns:", paste(head(nm, 12), collapse = " | "), "\n")
  sens[, gene := sub("\\|.*", "", exposure)]
  map <- unique(ann[, .(gene = sub("\\|.*", "", exposure), SYMBOL)])
  sens <- merge(sens, map, by = "gene", all.x = TRUE)
  gs <- sens[SYMBOL %in% GLYCO]
  if (nrow(gs)) {
    cat("\n=== glycolytic genes in the multi-SNP sensitivity set ===\n")
    print(gs[order(pval)][1:min(30, .N)], digits = 3)
    fwrite(gs, "34c_glycolysis_sensitivity_MR.tsv", sep = "\t")
  } else cat("\nno glycolytic gene in the sensitivity set\n")
}
