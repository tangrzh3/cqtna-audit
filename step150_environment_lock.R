## step150 -- 可执行环境记录：把实际用到的每个包的确切版本写成锁文件。
##
## ⚠ 这**不是** `renv` 产出的 `renv.lock`，本项目没有使用 `renv`。
## 它是**从当前活动库导出的、采用 renv.lock 同一 JSON 结构**的环境记录，
## 因此 `renv::restore()` 可以直接读它，但它的来源必须说清楚：
## 它记录的是"现在装着什么"，而现在装着的正是产出当前结果表的那一套
## —— 本轮的结果表（132–149）全部在这个环境下重跑过。
##
## 它**不能**声称描述更早的运行。若要一个真正的 renv.lock，
## 必须用 renv 从头重跑一遍，那是作者的决定（S40_LOCK.md §2）。
##
## 输出：150a_environment.lock（renv.lock 结构的 JSON）
##       150b_session.txt（sessionInfo 全文）
##       150c_console.log

arg <- commandArgs(trailingOnly = TRUE)
MR <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
      if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else getwd()
setwd(MR)

## 分析脚本里直接 library()/require()/:: 到的包
direct <- c("cqtna", "susieR", "data.table", "coloc", "Matrix",
            "arrow", "dplyr", "Seurat", "hdf5r",
            "TFBSTools", "motifmatchr", "JASPAR2020", "Biostrings",
            "BSgenome.Hsapiens.UCSC.hg38",
            "stats", "utils", "graphics", "grDevices", "methods")

ip <- utils::installed.packages()
have <- rownames(ip)
present <- intersect(direct, have)
missing <- setdiff(direct, have)

## 递归依赖闭包：锁文件若只列直接依赖，restore 出来的环境仍然不同
deps <- tryCatch(
  unique(unlist(tools::package_dependencies(present, db = ip, recursive = TRUE))),
  error = function(e) character(0))
all_pkgs <- sort(unique(c(present, intersect(deps, have))))

esc <- function(x) gsub('"', '\\\\"', x)
rec <- vapply(all_pkgs, function(p) {
  v <- as.character(ip[p, "Version"])
  pr <- as.character(ip[p, "Priority"])
  src <- if (!is.na(pr) && pr %in% c("base", "recommended")) "R" else "Repository"
  sprintf('    "%s": {\n      "Package": "%s",\n      "Version": "%s",\n      "Source": "%s"\n    }',
          esc(p), esc(p), esc(v), src)
}, character(1))

json <- c(
  "{",
  '  "R": {',
  sprintf('    "Version": "%s.%s",', R.version$major, R.version$minor),
  '    "Repositories": [',
  '      { "Name": "CRAN", "URL": "https://cloud.r-project.org" }',
  "    ]",
  "  },",
  '  "Packages": {',
  paste(rec, collapse = ",\n"),
  "  }",
  "}")
writeLines(json, "150a_environment.lock", useBytes = TRUE)

sink("150b_session.txt"); print(sessionInfo()); sink()

cat(sprintf("R %s.%s on %s\n", R.version$major, R.version$minor,
            R.version$platform))
cat(sprintf("packages recorded: %d (%d named directly, %d pulled in as dependencies)\n",
            length(all_pkgs), length(present), length(all_pkgs) - length(present)))
if (length(missing))
  cat(sprintf("⚠ named but NOT installed here, so NOT in the lock: %s\n",
              paste(missing, collapse = ", ")))
cat("\nversions of the packages the analysis calls directly:\n")
for (p in present)
  cat(sprintf("  %-30s %s\n", p, as.character(ip[p, "Version"])))
cat("\n⚠ this file is in renv.lock's format but was NOT produced by renv.\n")
cat("⚠ it records what is installed now, which is what produced the current\n")
cat("  result tables; it does not describe any earlier run.\n")
cat("\nwrote 150a_environment.lock / 150b_session.txt\n")
