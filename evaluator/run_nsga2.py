"""E2: NSGA-II tiga objektif (f1 biaya, f2 ECOST, f3 kompensasi) untuk satu sistem.

Contoh (di PC, 8 core):
  python run_nsga2.py --system ieee33 --seeds 0-9 --pop 100 --gen 150 --workers 8
Hasil per seed: results/<system>/seed<k>.npz (populasi akhir, front, riwayat front).
"""
import argparse, json, os, time
import numpy as np
from multiprocessing import Pool
from pymoo.core.problem import ElementwiseProblem
try:
    from pymoo.parallelization import StarmapParallelization
except ImportError:                                   # older pymoo
    from pymoo.core.problem import StarmapParallelization
from pymoo.core.sampling import Sampling
from pymoo.core.mutation import Mutation
from pymoo.core.callback import Callback
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.ux import UX
from pymoo.optimize import minimize

_SYS = {}


def get_sys(cfg):
    key = json.dumps(cfg, sort_keys=True)
    if key not in _SYS:
        from systems import System
        _SYS[key] = System(**cfg)
    return _SYS[key]


class Placement(ElementwiseProblem):
    def __init__(self, cfg, **kw):
        sy = get_sys(cfg)
        self.cfg = cfg
        super().__init__(n_var=sy.nvar, n_obj=3, n_ieq_constr=1, xl=0, xu=sy.xu,
                         vtype=int, **kw)

    def _evaluate(self, x, out, *args, **kwargs):
        o = get_sys(self.cfg).evaluate(np.asarray(x, int))
        out["F"] = [o["f1"], o["f2"], o["f3"]]
        out["G"] = [0.0 if o["rec_ok"] else 1.0]


class SparseSampling(Sampling):
    """Sparse initial plans (device density 0-40 %) plus the empty plan."""
    def _do(self, problem, n_samples, **kwargs):
        X = np.zeros((n_samples, problem.n_var), int)
        for i in range(1, n_samples):
            dens = np.random.uniform(0.0, 0.4)
            on = np.random.random(problem.n_var) < dens
            X[i, on] = [np.random.randint(1, problem.xu[j] + 1) for j in np.where(on)[0]]
        return X


class ResetMutation(Mutation):
    """Random reset of each gene with probability 1.5/n, biased towards 'no device'."""
    def _do(self, problem, X, **kwargs):
        X = X.copy()
        pm = 1.5 / problem.n_var
        for i in range(len(X)):
            for j in np.where(np.random.random(problem.n_var) < pm)[0]:
                X[i, j] = 0 if np.random.random() < 0.4 else np.random.randint(0, problem.xu[j] + 1)
        return X


class History(Callback):
    def __init__(self):
        super().__init__(); self.F = []
    def notify(self, algorithm):
        opt = algorithm.opt
        feas = opt.get("CV")[:, 0] <= 0
        self.F.append(opt.get("F")[feas].tolist())


def run(cfg, seed, pop, gen, workers, outdir):
    t0 = time.perf_counter()
    if workers > 1:
        pool = Pool(workers)
        problem = Placement(cfg, elementwise_runner=StarmapParallelization(pool.starmap))
    else:
        pool, problem = None, Placement(cfg)
    algo = NSGA2(pop_size=pop, sampling=SparseSampling(), crossover=UX(),
                 mutation=ResetMutation(), eliminate_duplicates=True)
    cb = History()
    res = minimize(problem, algo, ("n_gen", gen), seed=seed, callback=cb, verbose=False)
    if pool:
        pool.close(); pool.join()
    os.makedirs(outdir, exist_ok=True)
    feas = res.opt.get("CV")[:, 0] <= 0
    np.savez(os.path.join(outdir, f"seed{seed}.npz"), X=res.opt.get("X")[feas],
             F=res.opt.get("F")[feas], popX=res.pop.get("X"), popF=res.pop.get("F"),
             hist=np.array(json.dumps(cb.F)), secs=time.perf_counter() - t0,
             n_eval=res.algorithm.evaluator.n_eval, cfg=json.dumps(cfg))
    print(f"[{cfg['name']}] seed {seed}: front {feas.sum()} solusi, {res.algorithm.evaluator.n_eval} "
          f"evaluasi, {time.perf_counter()-t0:.0f} s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)          # ieee33 | ieee69 | rbts2 | rbts2_tr10
    ap.add_argument("--seeds", default="0-9")
    ap.add_argument("--pop", type=int, default=100)
    ap.add_argument("--gen", type=int, default=150)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--tau", type=float, default=60.0)   # manual switching time (min)
    ap.add_argument("--thr", type=float, default=1.0)    # compensation threshold (h)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    lo, hi = (int(v) for v in a.seeds.split("-")) if "-" in a.seeds else (int(a.seeds),) * 2
    cfg = dict(name=a.system, tau_ms_min=a.tau, thr_h=a.thr)
    tag = f"{a.system}_t{a.tau:g}_h{a.thr:g}"
    for sd in range(lo, hi + 1):
        run(cfg, sd, a.pop, a.gen, a.workers, os.path.join(a.out, tag))
