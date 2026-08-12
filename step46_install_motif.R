## Motif-analysis dependencies. TFBSTools, GenomicRanges, Biostrings and
## SummarizedExperiment are already present; these are the rest.
## BSgenome.Hsapiens.UCSC.hg38 is ~800 MB and is the slow part.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("motifmatchr", "JASPAR2020",
                       "BSgenome.Hsapiens.UCSC.hg38", "chromVAR"),
                     ask = FALSE, update = FALSE)
for (p in c("motifmatchr", "JASPAR2020", "BSgenome.Hsapiens.UCSC.hg38",
            "chromVAR", "TFBSTools")) {
  cat(sprintf("%-34s %s\n", p,
      ifelse(requireNamespace(p, quietly = TRUE),
             as.character(packageVersion(p)), "STILL MISSING")))
}
