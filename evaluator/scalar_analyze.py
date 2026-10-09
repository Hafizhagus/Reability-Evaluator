"""Analisis E4 dan E5.  python scalar_analyze.py e45_results/e4.json  (atau e5.json)"""
import json, sys
import numpy as np
from systems import System

fn = sys.argv[1]; rows = json.load(open(fn)); exp = "e4" if fn.endswith("e4.json") else "e5"
K = lambda c: json.dumps(c, sort_keys=True)
def score(v, ob):
    if ob == "social": return v["f1"] + v["f2"]
    if ob == "utility": return v["f1"] + v["f3"]
    return v["f2"] if ob.startswith("becost") else v["f3"]
def counts(g, sy):
    x, w = g[:sy.nsec], g[sy.nsec:]
    return f"MS {sum(v==1 for v in x)}/RCS {sum(v==2 for v in x)}/REC {sum(v==3 for v in x)}/tie {sum(v>0 for v in w)}"
out = {}
if exp == "e4":
    print("E4 (RBTS Bus 2): model tingkat beban untuk cek kelayakan vs acuan 4 kelompok")
    print("trafo | objektif | model | salah nilai ECOST rencana terpilih | regret keputusan")
    for tr in sorted({r["cfg"]["tr_repair"] for r in rows}):
        for ob in sorted({r["obj"] for r in rows}):
            R = [r for r in rows if r["cfg"]["tr_repair"] == tr and r["obj"] == ob]
            if not R: continue
            best_ref = min(score(r["atref"], ob) for r in R)
            for lm in ("peak", "avg", "ldc4"):
                Rm = [r for r in R if r["cfg"]["load_model"] == lm]
                if not Rm: continue
                xm = min(Rm, key=lambda r: score(r["own"], ob))
                val_err = (xm["own"]["f2"] - xm["atref"]["f2"]) / xm["atref"]["f2"]
                regret = (score(xm["atref"], ob) - best_ref) / best_ref
                out[f"{tr}|{ob}|{lm}"] = dict(valuation_error=val_err, regret=regret, g=xm["g"])
                print(f"  {tr:4g} h | {ob:8s} | {lm:4s} | {100*val_err:+7.2f}% | {100*regret:6.2f}%")
    # fixed transfer-reliant plan: MS on all non-head main sections + both ties manual
    print("\nRencana tetap 'MS di feeder utama + 2 tie manual': ECOST per model (k$/th)")
    for tr in sorted({r["cfg"]["tr_repair"] for r in rows}):
        vals = {}
        for lm in ("peak", "avg", "ldc4"):
            sy = System("rbts2", tr_repair=tr, load_model=lm)
            g = np.zeros(sy.nvar, int)
            for i, s in enumerate(sy.cand):
                if s in sy.net.sec and False: pass
            main = set(__import__("rbts2_data").build()[1]["main"])
            for i, s in enumerate(sy.cand):
                if s in main: g[i] = 1
            g[sy.nsec:] = 1
            vals[lm] = sy.evaluate(g)["f2"]
        print(f"  trafo {tr:4g} h: puncak {vals['peak']/1e3:8.1f} | rata-rata {vals['avg']/1e3:8.1f} | 4 kelompok {vals['ldc4']/1e3:8.1f} "
              f"| salah nilai puncak {100*(vals['peak']-vals['ldc4'])/vals['ldc4']:+.1f}%")
        out[f"fixed|{tr}"] = vals
else:
    print("E5: sensitivitas solusi representatif (terbaik dari semua seed)")
    for sysn in sorted({r["cfg"]["name"] for r in rows}):
        sy = System(sysn)
        base = {ob: min([r for r in rows if r["cfg"] == {"name": sysn} and r["obj"] == ob],
                        key=lambda r: score(r["own"], ob), default=None) for ob in ("social", "utility")}
        print(f"\n{sysn}")
        for cfg in sorted({K(r["cfg"]) for r in rows if r["cfg"]["name"] == sysn}):
            c = json.loads(cfg); lab = {k: v for k, v in c.items() if k != "name"} or "BASE"
            for ob in ("social", "utility"):
                R = [r for r in rows if K(r["cfg"]) == cfg and r["obj"] == ob]
                if not R: continue
                x = min(R, key=lambda r: score(r["own"], ob)); v = x["own"]
                gd = sum(1 for a, b in zip(x["g"], base[ob]["g"])) if False else \
                     (sum(1 for a, b in zip(x["g"], base[ob]["g"]) if a != b) if base[ob] else -1)
                print(f"  {str(lab):45s} {ob:7s}: f1 {v['f1']/1e3:6.1f} ECOST {v['f2']/1e3:7.1f} KOMP {v['f3']/1e3:6.2f} "
                      f">ambang {100*v['share_above']:3.0f}% | {counts(x['g'], sy)} | gen beda vs base {gd}")
                out[f"{cfg}|{ob}"] = dict(own=v, g=x["g"], genes_vs_base=gd)
json.dump(out, open(fn.replace(".json", "_analysis.json"), "w"), indent=1)
