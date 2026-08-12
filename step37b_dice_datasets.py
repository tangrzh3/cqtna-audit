"""List the DICE (Schmiedel_2018) datasets in the eQTL Catalogue, plus any other
study offering an activated CD4 T-cell condition -- those are the candidates for
an independent test of the activation-window eQTL."""
import json
import urllib.request

API = "https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size=1000"
req = urllib.request.Request(API, headers={"User-Agent": "Mozilla/5.0"})
d = json.load(urllib.request.urlopen(req, timeout=120))
print(f"total datasets in catalogue: {len(d)}\n")

print("=" * 92)
print("Schmiedel_2018  (DICE)")
print("=" * 92)
for r in d:
    if r["study_label"] == "Schmiedel_2018":
        print(f"  {r['dataset_id']:<12}{r['sample_group']:<38}"
              f"{r['condition_label']:<22}n={r['sample_size']:<5}{r['quant_method']}")

print("\n" + "=" * 92)
print("other datasets with an activated / stimulated T-cell condition")
print("=" * 92)
KEY = ("activ", "stim", "cd3", "cd28", "pha", "iav", "anti")
seen = set()
for r in d:
    blob = f"{r['sample_group']} {r['condition_label']} {r['tissue_label']}".lower()
    if any(k in blob for k in KEY) and ("t-cell" in blob or "t_cell" in blob
                                        or "tcell" in blob or "cd4" in blob or "cd8" in blob):
        k = (r["study_label"], r["sample_group"], r["condition_label"])
        if k in seen:
            continue
        seen.add(k)
        print(f"  {r['study_label']:<20}{r['dataset_id']:<12}{r['sample_group']:<34}"
              f"{r['condition_label']:<20}n={r['sample_size']:<5}{r['quant_method']}")
