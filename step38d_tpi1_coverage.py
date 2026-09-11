"""Which DICE datasets actually quantify TPI1? Without this, NaN in the D2 table
cannot be told apart from a genuine null."""
import json
import os
import sys
import time
import urllib.request
import urllib.parse
import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")

API = "https://www.ebi.ac.uk/eqtl/api/v2"
UA = {"User-Agent": "Mozilla/5.0"}
GENES = {"TPI1": "ENSG00000111669", "SPSB2": "ENSG00000111671"}


def get(path, **p):
    url = f"{API}/{path}?" + urllib.parse.urlencode(p)
    for a in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=120) as r:
                return json.load(r), None
        except Exception as ex:
            code = getattr(ex, "code", None)
            if a == 4:
                return None, code
            time.sleep(2 * (a + 1))
    return None, None


ds, _ = get("datasets/", size=1000)
dice = [d for d in ds if d["study_label"] == "Schmiedel_2018"
        and d["quant_method"] == "ge"]

rows = []
for d in dice:
    rec = dict(dataset=d["dataset_id"], cell=d["sample_group"],
               condition=d["condition_label"], n=int(d["sample_size"]))
    for g, ensg in GENES.items():
        a, code = get(f"datasets/{d['dataset_id']}/associations",
                      molecular_trait_id=ensg, size=1)
        time.sleep(0.4)
        if a:
            rec[f"{g}_quantified"] = True
            rec[f"{g}_median_tpm"] = a[0].get("median_tpm")
        else:
            rec[f"{g}_quantified"] = False
            rec[f"{g}_median_tpm"] = None
            rec[f"{g}_http"] = code
    rows.append(rec)
    print(f"  {d['sample_group']:<28}{d['condition_label']:<18}"
          f"TPI1={rec['TPI1_quantified']!s:<6}SPSB2={rec['SPSB2_quantified']}",
          flush=True)

t = pd.DataFrame(rows)
t.to_csv(os.path.join(MR, "38c_DICE_gene_coverage.tsv"), sep="\t", index=False)
print("\n" + "=" * 78)
print(f"TPI1  quantified in {t.TPI1_quantified.sum()}/{len(t)} DICE datasets")
print(f"SPSB2 quantified in {t.SPSB2_quantified.sum()}/{len(t)} DICE datasets")
print("=" * 78)
print(t[["cell", "condition", "TPI1_quantified", "TPI1_median_tpm",
         "SPSB2_quantified"]].to_string(index=False))
