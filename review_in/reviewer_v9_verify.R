## Third-party reviewer's own verification script, v9 round (2026-08-20).
##
## NOT written by this project. Kept because it is the independent recomputation
## behind Reviewer 1's item 4 -- the reviewer re-derived every Fisher p-value in
## 123d and 130a from the counts and confirmed they match the tables, which is
## what allowed them to separate "no arithmetic error" from the real limitation
## (loci are not independent tests, so the p-values are nominal).
##
## `root` below points at the reviewer's own extraction directory; change it to
## wherever you unpack review_packet_v9 before running.

root <- "C:/Users/a1197/AppData/Local/Temp/review_packet_v9_f5529a566e294d4486f4194903632bab"

d <- read.delim(file.path(root, "123d_fixed_anchor_full_grid.tsv"), check.names = FALSE)
p <- with(d, phyper(sig_known - 1, bg_known, bg_loci - bg_known,
                    sig_loci, lower.tail = FALSE))
f <- with(d, (sig_known / sig_loci) / (bg_known / bg_loci))
cat("123d max |P diff| =", max(abs(p - d$fisher_p), na.rm = TRUE),
    " max |fold diff| =", max(abs(f - d$fold), na.rm = TRUE), "\n")

m <- read.delim(file.path(root, "130a_multilist_main.tsv"), check.names = FALSE)
pm <- with(m, phyper(sig_known - 1, bg_known, bg_loci - bg_known,
                     sig_loci, lower.tail = FALSE))
fm <- with(m, (sig_known / sig_loci) / (bg_known / bg_loci))
fm[m$sig_known == 0] <- 0
cat("130a max |P diff| =", max(abs(pm - m$fisher_p), na.rm = TRUE),
    " max |fold diff| =", max(abs(fm - m$fold), na.rm = TRUE), "\n")

x <- subset(d, role == "main" & analysis == "main" &
              region_filter == "all_genome")
print(x[, c("cell", "fold", "fisher_p", "mismatch_p", "control")],
      row.names = FALSE)
cat("BH over six own-list P:\n")
print(p.adjust(x$fisher_p, "BH"))
cat("BH over five evaluable own-list P:\n")
print(p.adjust(x$fisher_p[x$control == "clean"], "BH"))

# The Fisher implementation compares significant and non-significant loci.
# This confirms the four cells sum back to the disjoint background partition.
stopifnot(all(d$bg_known - d$sig_known >= 0),
          all((d$bg_loci - d$bg_known) -
                (d$sig_loci - d$sig_known) >= 0))
