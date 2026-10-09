"""Analisis E3 (dan E6 pada 3 MS): regret*, harga sosial, sebaran seed, perubahan perangkat.
  python e3_analyze.py e3_results/ieee33_t60_h1.json"""
import json, sys
import numpy as np
from systems import System

fn = sys.argv[1]; rows = json.load(open(fn))
name = fn.split("/")[-1].split("_t")[0]
tau = float(fn.split("_t")[1].split("_h")[0]); thr = float(fn.split("_h")[1].split(".json")[0])
sy = System(name, tau_ms_min=tau, thr_h=thr)
def counts(g):
    x, w = g[:sy.nsec], g[sy.nsec:]
    return f"MS {sum(v==1 for v in x)}, RCS {sum(v==2 for v in x)}, REC {sum(v==3 for v in x)}, tie M/R {sum(v==1 for v in w)}/{sum(v==2 for v in w)}"
out = {}
print(f"{name} (manual {tau:g} menit, ambang {thr:g} jam)")
print("anggaran | regret* | regret per seed med [min-maks] | sebaran Komp maks | harga sosial | gen beda | >ambang E/K")
for n in sorted({r["budget"] for r in rows}):
    R = [r for r in rows if r["budget"] == n]
    xE = min(R, key=lambda r: r["f2"]); xK = min(R, key=lambda r: r["f3"])
    Kref, Eref = xK["f3"], xE["f2"]
    regret_star = (xE["f3"] - Kref) / Kref
    rs = [(r["f3"] - Kref) / Kref for r in R if r["obj"] == "ecost" and r["seed"] != "exact"]
    sp = [(r["f3"] - Kref) / Kref for r in R if r["obj"] == "comp" and r["seed"] != "exact"]
    social = (xK["f2"] - Eref) / Eref
    genes = sum(1 for a, b in zip(xE["g"], xK["g"]) if a != b)
    print(f"  {n:4g} MS | {100*regret_star:6.2f}% | {100*np.median(rs):6.2f}% [{100*min(rs):.2f}-{100*max(rs):.2f}] | "
          f"{100*max(sp):6.2f}% | {100*social:+7.2f}% | {genes:2d} | {100*xE['share_above']:.0f}%/{100*xK['share_above']:.0f}%")
    print(f"           ECOST-opt: {counts(xE['g'])} | Komp-opt: {counts(xK['g'])}")
    ent = dict(regret_star=regret_star, regret_seed=rs, spread_comp=sp, social=social, genes=genes,
               planE=sy.describe(xE["g"]), planK=sy.describe(xK["g"]),
               shareE=xE["share_above"], shareK=xK["share_above"])
    ex = {r["obj"]: r for r in R if r["seed"] == "exact"}
    if ex:   # E6: GA vs exact at the smallest budget
        for k, key in (("ecost", "f2"), ("comp", "f3")):
            best_ga = min((r[key] for r in R if r["obj"] == k and r["seed"] != "exact"), default=np.nan)
            ent[f"gap_ga_exact_{k}"] = (best_ga - ex[k][key]) / ex[k][key]
        print(f"           E6: selisih GA terbaik vs eksak: ECOST {100*ent['gap_ga_exact_ecost']:.3f}%, "
              f"Komp {100*ent['gap_ga_exact_comp']:.3f}%")
    out[str(n)] = ent
json.dump(out, open(fn.replace(".json", "_analysis.json"), "w"), indent=1)
