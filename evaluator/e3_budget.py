"""E3: rencana optimal-ECOST vs optimal-kompensasi pada anggaran setara 3/5/8/10/12 MS.
  python e3_budget.py --system ieee33 --budgets 3,5,8,10,12 --seeds 5 --workers 8
Per anggaran: GA + local search (langkah tunggal dan tukar) untuk tiap objektif dan seed;
pada 3 MS juga optimum eksak (enumerasi) untuk kedua objektif sekaligus (dipakai di E6).
Hasil: e3_results/<system>_t<tau>_h<thr>.json (setiap rencana dengan f1, f2, f3, indeks)."""
import argparse, itertools, json, os, random, time
from multiprocessing import Pool
import numpy as np

_SY = {}
def sysget(cfg):
    k = json.dumps(cfg, sort_keys=True)
    if k not in _SY:
        _SY.clear()                                   # keep one system per worker (bounded memory)
        from systems import System
        s = System(**cfg); s._cache = {}; _SY[k] = s
    return _SY[k]

def ev(sy, g):
    g = tuple(int(v) for v in g)
    if g not in sy._cache:
        if len(sy._cache) > 50000:
            sy._cache.clear()
        sy._cache[g] = sy.evaluate(np.array(g))
    return sy._cache[g]

def f1(sy, g):
    return sy.ev.cost(*sy.decode(g))

def obj(sy, g, B, kind):
    o = ev(sy, g)
    pen = (0 if o["rec_ok"] else 1e12) + (0 if o["f1"] <= B + 1e-6 else 1e12)
    return (o["f2"] if kind == "ecost" else o["f3"]) + pen

def repair(sy, g, B, rng):
    g = list(g)
    while f1(sy, g) > B + 1e-6:
        on = [i for i, v in enumerate(g) if v]; g[rng.choice(on)] = 0
    return g

def polish(sy, g, B, kind):
    best = list(g); vals = lambda i: range(int(sy.xu[i]) + 1)
    improved = True
    while improved:
        improved = False; cur = obj(sy, best, B, kind); moves = []
        for i in range(len(best)):
            for v in vals(i):
                if v != best[i]:
                    h = list(best); h[i] = v; moves.append(h)
        on = [i for i, v in enumerate(best) if v]; off = [j for j, v in enumerate(best) if not v]
        for i in on:
            for j in off:
                for v in list(vals(j))[1:]:
                    h = list(best); h[i] = 0; h[j] = v; moves.append(h)
        for h in moves:
            if f1(sy, h) <= B + 1e-6 and obj(sy, h, B, kind) < cur - 1e-9:
                best, improved = h, True; break
    return best

def ga(sy, B, kind, seed, pop=30, gens=30):
    rng = random.Random(seed); n = sy.nvar
    def rnd():
        g = [0] * n
        for _ in range(rng.randint(1, 8)):
            i = rng.randrange(n); g[i] = rng.randint(1, int(sy.xu[i]))
        return repair(sy, g, B, rng)
    P = [rnd() for _ in range(pop)]; key = lambda g: obj(sy, g, B, kind)
    for _ in range(gens):
        P.sort(key=key); nxt = [P[0], P[1]]
        while len(nxt) < pop:
            a = min(rng.sample(P, 3), key=key); b = min(rng.sample(P, 3), key=key)
            c = [a[i] if rng.random() < 0.5 else b[i] for i in range(n)]
            for i in range(n):
                if rng.random() < 1.5 / n:
                    c[i] = 0 if rng.random() < 0.4 else rng.randint(0, int(sy.xu[i]))
            nxt.append(repair(sy, c, B, rng))
        P = nxt
    return polish(sy, min(P, key=key), B, kind)

def exact_small(sy, B):
    """Every plan within B made of up to three MS-cost units (MS or manual tie) or one
    unit of cost <= B (RCS, REC, remote tie); returns the optimum of both objectives."""
    n = sy.nvar; best = {"ecost": (None, np.inf), "comp": (None, np.inf)}
    cands = [[0] * n]
    for k in (1, 2, 3):
        for comb in itertools.combinations(range(n), k):
            g = [0] * n
            for i in comb: g[i] = 1
            cands.append(g)
    for i in range(n):
        for v in range(2, int(sy.xu[i]) + 1):
            g = [0] * n; g[i] = v; cands.append(g)
    for g in cands:
        if f1(sy, g) > B + 1e-6: continue
        for kind in best:
            v = obj(sy, g, B, kind)
            if v < best[kind][1]: best[kind] = (g, v)
    return {k: v[0] for k, v in best.items()}

def task(args):
    cfg, n, kind, seed = args
    sy = sysget(cfg); B = n * sy.ms_ann; t0 = time.perf_counter()
    if kind == "exact":
        res = exact_small(sy, B); out = []
        for k, g in res.items():
            out.append(dict(budget=n, obj=k, seed="exact", g=g, **ev(sy, g)))
    else:
        g = ga(sy, B, kind, seed=10007 * n + 97 * seed + (0 if kind == "ecost" else 1))
        out = [dict(budget=n, obj=kind, seed=seed, g=g, **ev(sy, g))]
    for o in out:
        o["secs"] = time.perf_counter() - t0
        o.update({k: (float(v) if not isinstance(v, (list, str)) else v) for k, v in o.items()
                  if k not in ("g", "seed", "obj")})
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True); ap.add_argument("--budgets", default="3,5,8,10,12")
    ap.add_argument("--seeds", type=int, default=5); ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--seed-offset", type=int, default=0)   # e.g. 5 -> seeds 5..9 (new, not repeats)
    ap.add_argument("--objs", default="ecost,comp")
    ap.add_argument("--tau", type=float, default=60.0); ap.add_argument("--thr", type=float, default=1.0)
    ap.add_argument("--no-exact", action="store_true"); ap.add_argument("--out", default="e3_results")
    a = ap.parse_args()
    cfg = dict(name=a.system, tau_ms_min=a.tau, thr_h=a.thr)
    budgets = [int(v) for v in a.budgets.split(",")]
    objs = a.objs.split(",")
    tasks = [(cfg, n, k, s) for n in budgets for k in objs
             for s in range(a.seed_offset, a.seed_offset + a.seeds)]
    if not a.no_exact and 3 in budgets:
        tasks = [(cfg, 3, "exact", 0)] + tasks
    os.makedirs(a.out, exist_ok=True)
    fn = os.path.join(a.out, f"{a.system}_t{a.tau:g}_h{a.thr:g}.json")
    rows = json.load(open(fn)) if os.path.exists(fn) else []     # append to earlier runs
    done = {(float(r["budget"]), r["obj"], r["seed"]) for r in rows}
    exact_done = {float(r["budget"]) for r in rows if r["seed"] == "exact"}
    def pending(t):
        if t[2] == "exact":
            return float(t[1]) not in exact_done
        return (float(t[1]), t[2], t[3]) not in done
    tasks = [t for t in tasks if pending(t)]
    print(f"{len(tasks)} tugas tersisa (yang sudah ada di {fn} dilewati)", flush=True)
    t0 = time.perf_counter()
    with Pool(a.workers) as pool:
        for out in pool.imap_unordered(task, tasks):
            rows.extend(out); json.dump(rows, open(fn, "w"))
            print(f"  anggaran {out[0]['budget']:g} MS, {out[0]['obj']}, seed {out[0]['seed']}: selesai "
                  f"[{time.perf_counter()-t0:.0f} s]", flush=True)
    print(f"selesai: {fn}")
