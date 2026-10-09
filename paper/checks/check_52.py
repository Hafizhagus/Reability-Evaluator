"""Section 5.1/5.2 numbers recomputed from the evaluator and the front files."""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "evaluator"))
from systems import System
D = os.path.join(HERE, "..", "paper_data")
out = {}
for s in ("ieee33", "ieee69", "rbts2", "bawean"):
    S = System(s); net = S.net
    N = np.array(net.N); P = np.array(net.P) if hasattr(net, "P") else None
    info = dict(sections=len(net.sec), cand=S.nsec, ties=S.ntie, customers=int(N.sum()), nvar=S.nvar)
    try:
        info["length_km"] = round(float(sum(sec["length_km"] for sec in net.sec)), 2)
        info["peak_kw"] = round(float(sum(net.P_kw)) if hasattr(net, "P_kw") else float("nan"), 1)
    except Exception as e:
        info["err"] = str(e)
    base = S.evaluate(np.zeros(S.nvar, dtype=int))
    info["empty"] = {k: round(float(v), 4) for k, v in base.items()}
    info["comp_over_ecost_empty"] = round(base["f3"] / base["f2"], 4)
    df = pd.read_csv(os.path.join(D, f"front_{s}.csv")); F = df[["f1", "f2", "f3"]].to_numpy(); X = df.iloc[:, 3:].to_numpy().astype(int)
    iS, iU = int(np.argmin(F[:, 0] + F[:, 1])), int(np.argmin(F[:, 0] + F[:, 2]))
    oS, oU = S.evaluate(X[iS]), S.evaluate(X[iU])
    info["front_size"] = len(F)
    info["soc"] = dict(f1=oS["f1"], f2=oS["f2"], f3=oS["f3"], share=oS["share_above"], plan=S.describe(X[iS]))
    info["util"] = dict(f1=oU["f1"], f2=oU["f2"], f3=oU["f3"], share=oU["share_above"], plan=S.describe(X[iU]))
    info["invest_ratio"] = (oS["f1"] / oU["f1"]) if oU["f1"] > 0 else None
    info["ecost_increase_pct"] = 100 * (oU["f2"] - oS["f2"]) / oS["f2"]
    # front values must equal a fresh evaluation
    dmax = 0.0
    for i in (iS, iU, 0, len(F) // 2, len(F) - 1):
        o = S.evaluate(X[i]); dmax = max(dmax, abs(o["f1"] - F[i, 0]), abs(o["f2"] - F[i, 1]), abs(o["f3"] - F[i, 2]))
    info["front_vs_reeval_maxabs"] = dmax
    out[s] = info
    print(f"\n== {s}: seksi {info['sections']} (kandidat {info['cand']}), tie {info['ties']}, pelanggan {info['customers']}, gen {info['nvar']}, front {info['front_size']} solusi")
    print("   tanpa perangkat baru:", {k: round(v, 3) for k, v in base.items()}, "| kompensasi/ECOST =", info["comp_over_ecost_empty"])
    print(f"   sosial : f1 {oS['f1']/1e3:.2f} f2 {oS['f2']/1e3:.2f} f3 {oS['f3']/1e3:.2f} share {100*oS['share_above']:.1f}%  {S.describe(X[iS])}")
    print(f"   utilitas: f1 {oU['f1']/1e3:.2f} f2 {oU['f2']/1e3:.2f} f3 {oU['f3']/1e3:.2f} share {100*oU['share_above']:.1f}%  {S.describe(X[iU])}")
    if info["invest_ratio"]:
        print(f"   rasio investasi sosial/utilitas {info['invest_ratio']:.2f} | ECOST utilitas lebih tinggi {info['ecost_increase_pct']:.1f}%")
    print(f"   selisih front vs evaluasi ulang (maks) {dmax:.2e}")
json.dump(out, open(os.path.join(HERE, "check_52.json"), "w"), indent=1, default=float)
