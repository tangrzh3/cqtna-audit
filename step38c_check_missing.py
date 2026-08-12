"""Is TPI1 x rs12302749 genuinely absent from DICE stimulated CD4, or did the
query miss it? A missing pair and a null association mean different things."""
import json
import time
import urllib.request
import urllib.parse
import pandas as pd

API = "https://www.ebi.ac.uk/eqtl/api/v2"
UA = {"User-Agent": "Mozilla/5.0"}
NAIVE, STIM = "QTD000479", "QTD000484"
TPI1 = "ENSG00000111669"
RSID, POS = "rs12302749", 6867132


def get(path, **p):
    url = f"{API}/{path}?" + urllib.parse.urlencode(p)
    for a in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=120) as r:
                return json.load(r)
        except Exception as ex:
            if a == 4:
                print(f"   FAIL {p} {type(ex).__name__} {getattr(ex,'code','')}")
                return None
            time.sleep(3 * (a + 1))


for tag, ds in (("naive", NAIVE), ("stim4h", STIM)):
    print("=" * 74)
    print(f"{tag}  ({ds})")
    print("=" * 74)

    a = get(f"datasets/{ds}/associations", rsid=RSID, size=500)
    if a is None:
        continue
    df = pd.DataFrame(a)
    print(f"  variant {RSID}: {len(df)} gene associations returned")
    if len(df):
        print(f"    genes: {sorted(df.molecular_trait_id.unique())[:12]}")
        print(f"    TPI1 present? {TPI1 in set(df.molecular_trait_id)}")
    time.sleep(0.5)

    # is TPI1 quantified in this dataset at all?
    b = get(f"datasets/{ds}/associations", molecular_trait_id=TPI1, size=5)
    time.sleep(0.5)
    if b:
        bb = pd.DataFrame(b)
        print(f"  TPI1 as a trait: {len(bb)} records on first page, "
              f"median_tpm={bb.median_tpm.iloc[0]}")
    else:
        print("  TPI1 as a trait: NOT QUANTIFIED in this dataset")

    # is the variant present at all, regardless of gene?
    c = get(f"datasets/{ds}/associations", pos=f"12:{POS}-{POS}", size=100)
    time.sleep(0.5)
    if c:
        cc = pd.DataFrame(c)
        print(f"  position 12:{POS}: {len(cc)} records | "
              f"rsids={sorted(set(cc.rsid.dropna()))[:5]}")
        print(f"    genes at this position: {sorted(cc.molecular_trait_id.unique())[:12]}")
        if TPI1 in set(cc.molecular_trait_id):
            h = cc[cc.molecular_trait_id == TPI1].iloc[0]
            print(f"    *** TPI1 FOUND by position: p={h.pvalue} beta={h.beta} "
                  f"maf={h.maf} ***")
    else:
        print(f"  position 12:{POS}: no records")

    # what is the strongest TPI1 cis signal in this dataset near the gene?
    d = get(f"datasets/{ds}/associations", molecular_trait_id=TPI1,
            pos=f"12:{POS-250000}-{POS+250000}", size=1000)
    time.sleep(0.5)
    if d:
        dd = pd.DataFrame(d)
        dd["pvalue"] = dd.pvalue.astype(float)
        best = dd.loc[dd.pvalue.idxmin()]
        print(f"  strongest TPI1 signal within +/-250 kb ({len(dd)} variants): "
              f"p={best.pvalue:.3g} at {best.position} ({best.rsid})")
