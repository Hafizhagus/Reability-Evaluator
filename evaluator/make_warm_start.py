"""Titik awal E5 (warm_start.json).
  python make_warm_start.py <folder front_<sistem>.csv> [--f2-scale 1.0] [--prev e5.json e5bw.json]
Isi per sistem:
  social / utility : optimum sosial (min f1+f2) dan utilitas (min f1+f3) pada front gabungan E2;
  bila --prev diberikan, rencana terbaik kasus dasar dari run satu-objektif sebelumnya ikut
  dimasukkan, dan rencana terbaik tiap varian disimpan di "variants".
--f2-scale: pengali f2 untuk berkas yang dihitung dengan faktor konversi ECOST lain
            (hasil lama dengan faktor 1,6543 -> 164.2/158.4)."""
import argparse, json, os
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("folder"); ap.add_argument("--f2-scale", type=float, default=1.0)
ap.add_argument("--prev", nargs="*", default=[])
a = ap.parse_args()
prev = [r for fn in a.prev for r in json.load(open(fn))]
score = {"social": lambda o: o["f1"] + a.f2_scale * o["f2"], "utility": lambda o: o["f1"] + o["f3"]}
out = {}
for s in ("ieee33", "ieee69", "rbts2", "bawean"):
    df = pd.read_csv(os.path.join(a.folder, f"front_{s}.csv"))
    F = df[["f1", "f2", "f3"]].to_numpy(); F[:, 1] *= a.f2_scale
    X = df.iloc[:, 3:].to_numpy().astype(int)
    i, j = int(np.argmin(F[:, 0] + F[:, 1])), int(np.argmin(F[:, 0] + F[:, 2]))
    W = dict(social=[[int(v) for v in X[i]]], utility=[[int(v) for v in X[j]]], variants={})
    print(f"{s:7s} front: sosial indeks {i} ({(F[i,0]+F[i,1])/1e3:.3f} kUSD/th), utilitas indeks {j} ({(F[j,0]+F[j,2])/1e3:.3f} kUSD/th)")
    cfgs = sorted({json.dumps(r["cfg"], sort_keys=True) for r in prev if r["cfg"].get("name") == s
                   and not ({"load_model", "tr_repair"} & set(r["cfg"]))})
    for ck in cfgs:
        cfg = json.loads(ck); vkey = json.dumps({k: v for k, v in cfg.items() if k != "name"}, sort_keys=True)
        for ob in ("social", "utility"):
            R = [r for r in prev if r["cfg"] == cfg and r["obj"] == ob]
            if not R:
                continue
            g = [int(v) for v in min(R, key=lambda r: score[ob](r["own"]))["g"]]
            if vkey == "{}":
                if g not in W[ob]:
                    W[ob].append(g)
            else:
                W["variants"].setdefault(vkey, {})[ob] = [g]
    print(f"        rencana awal kasus dasar: sosial {len(W['social'])}, utilitas {len(W['utility'])}; varian: {len(W['variants'])}")
    out[s] = W
json.dump(out, open("warm_start.json", "w"))
print("tulis warm_start.json")
