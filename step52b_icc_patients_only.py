"""Is the glycolysis ICC a person-level trait, or just healthy-vs-patient?

Step 52 gave ICC 0.104 against a random floor of 0.023, which looks like real
person-level structure. But the two highest-scoring donors are exactly the two
healthy donors (HD1 +0.131, HD2 +0.164) while all eight patients fall between
-0.082 and +0.027 -- and the healthy donors are also the only two females and
were plausibly processed separately. A two-group difference, a sex difference
and a batch difference are all confounded here, and any of them inflates the ICC
without there being a person-level trait.

A biomarker needs variation between individuals WITHIN a clinically homogeneous
group, so the ICC is recomputed among the eight patients alone. The random-module
floor is recomputed the same way, since a smaller, more homogeneous set changes
the floor too.
"""
import numpy as np
import pandas as pd

MR = r"D:/R_ex/MR"
d = pd.read_csv(f"{MR}/52c_cell_scores.tsv.gz", sep="\t")
allicc = pd.read_csv(f"{MR}/52a_icc_all.tsv", sep="\t")


def icc(values, donor):
    df = pd.DataFrame(dict(v=values, d=donor)).dropna()
    grp = df.groupby("d").v
    k = grp.size()
    if len(k) < 3:
        return np.nan
    grand = df.v.mean()
    msb = (k * (grp.mean() - grand) ** 2).sum() / (len(k) - 1)
    msw = df.groupby("d").v.apply(lambda s: ((s - s.mean()) ** 2).sum()).sum()
    msw = msw / (len(df) - len(k))
    k0 = (len(df) - (k ** 2).sum() / len(df)) / (len(k) - 1)
    vb = max((msb - msw) / k0, 0)
    return vb / (vb + msw) if (vb + msw) > 0 else np.nan


pat = d[d.donor.str.startswith("P")]
print(f"patients only: {pat.donor.nunique()} donors, {len(pat):,} cells")
print(f"  all 10 donors      ICC(glycolysis) = "
      f"{icc(d.Glyco.values, d.donor.values):.4f}")
print(f"  8 patients only    ICC(glycolysis) = "
      f"{icc(pat.Glyco.values, pat.donor.values):.4f}")
print(f"  2 healthy donors   mean Glyco = "
      f"{d[d.donor.str.startswith('HD')].Glyco.mean():+.4f}")
print(f"  8 patients         mean Glyco = {pat.Glyco.mean():+.4f}")

print("\nper-donor means, patients only:")
pm = pat.groupby("donor").Glyco.agg(["mean", "size"]).sort_values("mean")
print(pm.round(4).to_string())
print(f"  spread across patients: {pm['mean'].min():+.4f} to "
      f"{pm['mean'].max():+.4f}  (range {pm['mean'].max()-pm['mean'].min():.4f})")

print("\nNOTE: the random-module floor from Step 52 was computed over all 10")
print("donors and is not the right comparator for the patients-only ICC; a")
print("floor recomputed within patients is needed before calling this positive.")
print(f"  10-donor random floor was {allicc[allicc.kind=='random_floor'].icc.mean():.4f}"
      f" (95th pct {np.percentile(allicc[allicc.kind=='random_floor'].icc, 95):.4f})")
