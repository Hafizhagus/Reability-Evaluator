"""Bawean: every single addition (one MS, RCS or REC on one candidate section) against the existing configuration."""
import os, sys, time, json
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "evaluator"))
from systems import System
_S = None
def init():
    global _S
    _S = System("bawean")
def ev(job):
    i, v = job
    g = np.zeros(_S.nvar, dtype=int); g[i] = v
    o = _S.evaluate(g)
    return i, v, float(o["f1"]), float(o["f2"]), float(o["f3"]), float(o["share_above"]), bool(o["rec_ok"]), _S.describe(g)
if __name__ == "__main__":
    S = System("bawean"); z = np.zeros(S.nvar, dtype=int); base = S.evaluate(z)
    jobs = [(i, v) for i in range(S.nsec) for v in (1, 2, 3)]
    t0 = time.time()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4, initializer=init) as pool:
        R = pool.map(ev, jobs, chunksize=4)
    ok = [r for r in R if r[6]]
    bs = min(ok, key=lambda r: r[2] + r[3]); bu = min(ok, key=lambda r: r[2] + r[4])
    zero = sorted((r[2], r[7]) for r in ok if r[5] < 1e-12)
    res = dict(n=len(R), n_rec_violation=len(R) - len(ok), seconds=time.time() - t0,
               base_soc=base["f1"] + base["f2"], base_util=base["f1"] + base["f3"],
               best_soc=(bs[2] + bs[3], bs[7]), best_util=(bu[2] + bu[4], bu[7]),
               n_improve_soc=sum(1 for r in ok if r[2] + r[3] < base["f1"] + base["f2"] - 1e-9),
               n_improve_util=sum(1 for r in ok if r[2] + r[4] < base["f1"] + base["f3"] - 1e-9),
               zero_share=zero[:10], n_zero_share=len(zero),
               min_share_single_ms=min(r[5] for r in ok if r[1] == 1))
    json.dump(res, open(os.path.join(HERE, "check_bawean_single.json"), "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))
