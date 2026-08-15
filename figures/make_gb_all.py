"""Build the whole GB nine-figure set in one command.

    python figures/make_gb_all.py

九张图各自独立可跑；本脚本只是把它们串起来，保证投稿前能一次重现全套，
并在结尾核对 PDF/PNG 是否都落了盘。

⚠ 这套是 `MANUSCRIPT_GB.md` 的**九图**编号。全文源
`MANUSCRIPT_v2_dual_thread.md` 的**八图**方案用的是
`FigA`–`FigG` / `Fig5_compartment` / `Fig6_glycolysis_TPI1` / `Fig7_crosscancer` /
`Fig9_part2`，两套并存，互不覆盖。对照见 `manuscript/FIGURES_GB_mapping.md`。

⚠ Fig 2 依赖 `step119_grid_rebuild.py`，Fig 6 依赖 `step120_tpi1_window.py`。
   两个 step 的产物（`119a`/`119c`/`120a`）已在仓库里；若重跑数据链须先跑它们。
"""
import importlib
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

FIGURES = [
    ("Fig 1", "make_gb_fig1_locus_attribution", "Fig1_locus_attribution"),
    ("Fig 2", "make_gb_fig2_generality", "Fig2_generality"),
    ("Fig 3", "make_gb_fig3_power_stability", "Fig3_power_stability"),
    ("Fig 4", "make_gb_fig4_coloc_heidi", "Fig4_coloc_heidi"),
    ("Fig 5", "make_gb_fig5_instrument_ladder", "Fig5_instrument_ladder"),
    ("Fig 6", "make_gb_fig6_tpi1_window", "Fig6_tpi1_window"),
    ("Fig 7", "make_gb_fig7_axis_chromatin", "Fig7_axis_chromatin"),
    ("Fig 8", "make_gb_fig8_compartment", "Fig8_compartment"),
    ("Fig 9", "make_gb_fig9_patients", "Fig9_patients"),
]


def main():
    print("=" * 78)
    print("GB nine-figure set")
    print("=" * 78)
    failed = []
    for label, mod_name, stem in FIGURES:
        t0 = time.time()
        print(f"\n--- {label}  ({mod_name}) " + "-" * (44 - len(mod_name)))
        try:
            mod = importlib.import_module(mod_name)
            mod.main()
            print(f"    {time.time() - t0:.1f}s")
        except Exception as exc:                      # noqa: BLE001
            failed.append((label, repr(exc)))
            print(f"    *** FAILED: {exc!r}")

    print("\n" + "=" * 78)
    missing = []
    for label, _, stem in FIGURES:
        for ext in (".pdf", ".png"):
            p = os.path.join(HERE, stem + ext)
            if os.path.exists(p):
                print(f"  {label:<7} {stem + ext:<34} {os.path.getsize(p) / 1024:8.0f} KB")
            else:
                missing.append(stem + ext)
                print(f"  {label:<7} {stem + ext:<34} *** MISSING ***")

    print("=" * 78)
    if failed:
        print(f"{len(failed)} figure(s) failed to build:")
        for label, err in failed:
            print(f"  {label}: {err}")
    if missing:
        print(f"{len(missing)} output file(s) missing.")
    if not failed and not missing:
        print("All nine figures built, PDF and PNG present for each.")
    return 1 if (failed or missing) else 0


if __name__ == "__main__":
    sys.exit(main())
