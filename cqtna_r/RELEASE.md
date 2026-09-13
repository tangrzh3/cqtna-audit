# 冻结投稿版本的操作清单

版本 **0.3.0**，建议标签 **`v0.3.0-paper`**。
**仓库**：https://github.com/tangrzh3/cqtna-audit


---

## 1. 三种"检查通过"不能混写

投稿材料里必须分开陈述，它们的强度完全不同：

| 层级 | 含义 | 本地状态 |
|---|---|---|
| **本地测试通过** | `testthat::test_local()` 在一台机器、一个 R 版本上跑过 | ✅ 见 §2 |
| **`R CMD check` 通过** | 打包、文档、示例、依赖、测试的完整检查 | ✅ 见 §2 |
| **跨平台 CI 通过** | Linux / macOS / Windows × R release / oldrel | ❌ **未做**——需要 GitHub 仓库 |

⚠ 在拿到 CI 记录之前，**不得**在论文里写"tested across platforms"。
现在能写的只有前两条，且应写明 R 版本与操作系统。

## 2. 本地检查记录（2026-08-18）

```
平台        Windows 11, x86_64-w64-mingw32
R           4.4.1 (2024-06-14 ucrt)
testthat    3.3.2

R CMD check cqtna_0.3.0.tar.gz --no-manual   ->  Status: OK   (0 error / 0 warning / 0 note)
testthat::test_local()                        ->  FAIL 0 | WARN 3 | SKIP 0 | PASS 204
```

三条 warning 全部是**包在按设计工作**的证据，不是缺陷：
已知位点名单缺 `source` 列的循环性提醒（×2）、单连锁串联提醒（×1）。

⚠ 冻结前重跑一次并把**完整输出**存进 release 附件，不要只存一行结论。

## 3. 冻结步骤

1. **建 GitHub 仓库**，把 `cqtna_r/` 作为子目录或独立仓库推上去。
2. **替换占位符**——四个文件里都有 `OWNER`，还有作者名与邮箱：
   ```bash
   grep -rn "OWNER\|noreply@example.com\|MR audit project" \
        cqtna_r/DESCRIPTION cqtna_r/README.md cqtna_r/CITATION.cff cqtna_r/.zenodo.json
   ```
3. **等 CI 跑绿**（`.github/workflows/R-CMD-check.yaml` 已备好 4 平台矩阵），
   把日志存下来。
4. **打标签** `v0.3.0-paper` 并发 release。
5. **Zenodo**：在 Zenodo 打开该仓库开关后再发 release，DOI 自动铸出。
6. **把 commit SHA 与 DOI 写进论文** 的 code availability。
7. **存档一份完整运行记录**：输入夹具、`cqtna_report()` 输出、
   `sessionInfo()`、以及跑过的命令。

```r
# 冻结时一并存下
writeLines(capture.output(sessionInfo()), "sessioninfo.txt")
```

## 4. 论文里该怎么写（建议措辞）

> cqtna (v0.3.0, DOI: ⟨Zenodo DOI⟩, commit ⟨SHA⟩) is a research companion
> implementing the automatable diagnostics used in this study. It is not a
> validated general-purpose target-nomination platform. It does not perform the
> Mendelian randomization itself, nor colocalisation with multiple-signal
> modelling, nor cell-level lineage matching, nor code-by-code verification of
> endpoint definitions across GWAS releases, and it does not issue
> target-supported conclusions; its outputs are intended for auditing and
> sensitivity analysis. `R CMD check` passes without errors, warnings or notes on
> R 4.4.1; the regression suite asserts this study's own reported values, and a
> second suite checks the R implementation against the Python reference
> implementation included in the deposit.

⚠ 若届时 CI 已跑绿，可在末句后加"and continuous integration covers Linux, macOS
and Windows on current and previous R releases"——**在此之前不要写**。

## 5. 冻结前仍未解决的两件事

1. **位点归属口径**（见 `CONVENTION_REVIEW_PACKET.md`）——
   已送第三方复核，结论可能要求改位点定义并重算论文数字。
   **在复核结论回来之前不要冻结。**
2. **密度匹配置换的结果与 Fisher 不一致**：CD4 那格 Fisher P = 0.028、
   置换经验 P = 0.111。正文若引 Fisher，须同时报置换结果。
