"""E4 dan E5: optimasi skalar (GA + local search) untuk banyak konfigurasi, paralel.
  python scalar_opt.py --exp e4 --workers 8      (tingkat beban & waktu trafo, RBTS)
  python scalar_opt.py --exp e5 --workers 8      (sensitivitas, tiga sistem)
  python scalar_opt.py --exp e5bw                (sensitivitas, Bawean)
Opsi: --objs social (hanya objektif tertentu), --out <folder>, --seeds 5.
Objektif: social = f1+f2, utility = f1+f3, becost:<n> / bcomp:<n> = ECOST / kompensasi
dengan anggaran setara n MS. Setiap rencana hasil juga dinilai dengan model acuan (v4 base)
agar salah nilai dan regret keputusan bisa dihitung di scalar_analyze.py."""
import argparse, json, os, time
from multiprocessing import Pool
import numpy as np
import e3_budget as E

def spec(exp, seeds, objs=None, all_variants=False):
    """objs: keep only these objectives (None = all). Variants that change only the
    compensation (threshold, tariff) leave f1+f2 untouched, so their societal runs would
    repeat the base case with the same seeds; they are skipped unless all_variants is set."""
    T = []
    if exp == "e4":
        for tr in (10, 25, 50, 100, 200):
            for lm in ("peak", "avg", "ldc4"):
                cfg = dict(name="rbts2", tr_repair=tr, load_model=lm)
                ref = dict(name="rbts2", tr_repair=tr, load_model="ldc4")
                for ob in ("social", "becost:5"):
                    T += [(cfg, ref, ob, s) for s in range(seeds)]
    elif exp in ("e5", "e5bw"):
        for sysn in (("ieee33", "ieee69", "rbts2") if exp == "e5" else ("bawean",)):
            base = dict(name=sysn)
            variants = [dict(tau_ms_min=45.0), dict(tau_ms_min=90.0), dict(cost_ratio=1.3),
                        dict(cost_ratio=2.0), dict(thr_h=2.0), dict(thr_h=6.0),
                        dict(tariff_scale={"residential": 2.0}), dict(tariff_scale={"residential": 0.5}),
                        dict(clpu_a=0.6), dict(clpu_a=1.0), dict()]
            for v in variants:
                cfg = dict(base, **v)
                for ob in ("social", "utility"):
                    if ob == "social" and not all_variants and ({"thr_h", "tariff_scale"} & set(v)):
                        continue
                    T += [(cfg, cfg, ob, s) for s in range(seeds)]
    if objs:
        T = [t for t in T if t[2] in objs]
    return T

WARM = json.load(open("warm_start.json")) if os.path.exists("warm_start.json") else {}


def warm_plans(cfg, ob):
    """Starting plans of an E5 run (see make_warm_start.py): the optimum of the E2 union front
    and, if listed, the best plans of earlier single-objective runs of the base case and of
    this variant. Returns a list of gene vectors without duplicates."""
    W = WARM.get(cfg["name"], {})
    vkey = json.dumps({k: v for k, v in cfg.items() if k != "name"}, sort_keys=True)
    cand = []
    for src in (W.get(ob), W.get("variants", {}).get(vkey, {}).get(ob)):
        if not src:
            continue
        cand += [src] if isinstance(src[0], (int, float)) else list(src)   # one plan or several
    out = []
    for g in cand:
        g = [int(v) for v in g]
        if g not in out:
            out.append(g)
    return out


def polish_free(sy, g, key):
    """Local search without budget: single-gene moves and swap moves (remove one device,
    add one elsewhere), first improvement, until no move improves."""
    g = list(g); n = sy.nvar
    improved = True
    while improved:
        improved = False; cur = key(g)
        moves = []
        for i in range(n):
            for v in range(int(sy.xu[i]) + 1):
                if v != g[i]:
                    h = list(g); h[i] = v; moves.append(h)
        on = [i for i in range(n) if g[i]]; off = [j for j in range(n) if not g[j]]
        for i in on:
            for j in off:
                for v in range(1, int(sy.xu[j]) + 1):
                    h = list(g); h[i] = 0; h[j] = v; moves.append(h)
        for h in moves:
            if key(h) < cur - 1e-9:
                g, improved = h, True
                break
    return g


def value(sy, g, ob):
    o = E.ev(sy, g)
    if ob == "social": return o["f1"] + o["f2"]
    if ob == "utility": return o["f1"] + o["f3"]
    return o["f2"] if ob.startswith("becost") else o["f3"]

def run(args):
    cfg, ref, ob, seed = args
    sy = E.sysget(cfg); t0 = time.perf_counter()
    if ob in ("social", "utility"):
        import random
        key = lambda g: value(sy, g, ob) + (0 if E.ev(sy, g)["rec_ok"] else 1e12)
        rng = random.Random(7919 * seed + (0 if ob == "social" else 1))
        n = sy.nvar
        P = []
        # E5 only: start from the best-known plans of the case, so every variant is at least
        # as good as the base plan under its own model. E4 is not warm-started, to avoid
        # biasing the load-model comparison towards the four-level reference.
        is_e5 = not ({"load_model", "tr_repair"} & set(cfg))
        if is_e5:
            P.extend(warm_plans(cfg, ob))
        while len(P) < 30:
            g = [0] * n
            for _ in range(rng.randint(0, int(0.5 * n))):
                i = rng.randrange(n); g[i] = rng.randint(1, int(sy.xu[i]))
            P.append(g)
        for _ in range(40):
            P.sort(key=key); nxt = [P[0], P[1]]
            while len(nxt) < 30:
                a = min(rng.sample(P, 3), key=key); b = min(rng.sample(P, 3), key=key)
                c = [a[i] if rng.random() < 0.5 else b[i] for i in range(n)]
                for i in range(n):
                    if rng.random() < 1.5 / n:
                        c[i] = 0 if rng.random() < 0.4 else rng.randint(0, int(sy.xu[i]))
                nxt.append(c)
            P = nxt
        g = polish_free(sy, min(P, key=key), key)
    else:
        nms = int(ob.split(":")[1]); B = nms * sy.ms_ann
        g = E.ga(sy, B, "ecost" if ob.startswith("becost") else "comp", seed=31 * seed + nms)
    own = E.ev(sy, g); sref = E.sysget(ref); refv = E.ev(sref, g)
    return dict(cfg=cfg, ref=ref, obj=ob, seed=seed, g=[int(v) for v in g],
                own={k: float(v) for k, v in own.items()}, atref={k: float(v) for k, v in refv.items()},
                secs=time.perf_counter() - t0)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True); ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--workers", type=int, default=os.cpu_count()); ap.add_argument("--out", default="e45_results")
    ap.add_argument("--limit", type=int, default=0)      # for quick tests
    ap.add_argument("--objs", default="")                # e.g. "social" or "utility,becost:5"; empty = all
    ap.add_argument("--all-variants", action="store_true")   # also run the redundant societal variants
    a = ap.parse_args()
    T = spec(a.exp, a.seeds, objs=[o for o in a.objs.split(",") if o], all_variants=a.all_variants)
    if a.limit: T = T[:a.limit]
    os.makedirs(a.out, exist_ok=True); fn = os.path.join(a.out, f"{a.exp}.json")
    rows = json.load(open(fn)) if os.path.exists(fn) else []
    done = {(json.dumps(r["cfg"], sort_keys=True), r["obj"], r["seed"]) for r in rows}
    T = [t for t in T if (json.dumps(t[0], sort_keys=True), t[2], t[3]) not in done]
    print(f"{a.exp}: {len(T)} tugas tersisa", flush=True); t0 = time.perf_counter()
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(run, T):
            rows.append(r); json.dump(rows, open(fn, "w"))
            print(f"  {r['cfg']} {r['obj']} seed {r['seed']} [{time.perf_counter()-t0:.0f} s]", flush=True)
