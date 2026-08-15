"""GB Fig 8 — compartment attribution.

九图口径里唯一**内容完全吻合、只差一个文件名**的一张：
`make_fig5_compartment.py` 的五个面板（a 细胞类型表达 · b 患者内配对 ·
c 预设阳性对照 · d 空间语境 · e bulk 生存是组成）正是 GB 图注要的东西。

因此这里**不重画**，只把同一个 figure 对象另存为 GB 的编号。
这样两套编号永远不会各画各的而悄悄漂移。

⚠ `FigE_spatial_compartment` 是更窄的**纯空间**版本，不是这一张，别选错。
⚠ 导入 `make_fig5_compartment` 会顺带重写一次 `Fig5_compartment.pdf/png`
   （内容与原来完全相同）——全文源八图方案的那张仍然在，未被覆盖成别的东西。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_fig5_compartment as F          # noqa: E402  执行即构图，fig 留在模块层

OUT = os.path.join(r"D:/R_ex/MR", "figures")


def main():
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        F.fig.savefig(os.path.join(OUT, "Fig8_compartment" + ext), **kw)
    n_panels = len(F.fig.axes)
    print(f"Fig8_compartment ok  |  re-saved from make_fig5_compartment "
          f"({n_panels} panels), content identical to Fig5_compartment")


if __name__ == "__main__":
    main()
