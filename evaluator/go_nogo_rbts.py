"""Uji go/no-go di RBTS Bus 2 (sistem uji kedua).

Model S : snapshot (beban puncak, tanpa CLPU)
Model F : CLPU + variasi beban (4 tingkat dari kurva lama beban RBTS, Billinton dkk.
          1989 Tabel I; waktu gangguan seragam sepanjang tahun) + pemulihan bertahap
Keputusan: perangkat di 32 seksi saluran (10 seksi utama non-kepala + 22 lateral) dan
2 titik normally-open (dipakai sebagai kandidat tie). PARAMETER CLPU MASIH SEMENTARA.
"""
import json, random, sys, time, pickle
import numpy as np
from network import Network, NONE, MS, RCS, REC, DEVICE_NAME
from reliability import Evaluator, default_params
import rbts2_data as R

# ---- Billinton et al. 1989, Table I (100-point load duration curve). Four OCR
# misreads in the source scan were corrected using monotonicity of the curve.
LDC = [(1.0000, 0.0000), (0.9733, 0.0006), (0.9466, 0.0024), (0.9199, 0.0076), (0.8931, 0.0160),
       (0.8664, 0.0333), (0.8397, 0.0614), (0.8130, 0.1004), (0.7863, 0.1452), (0.7596, 0.1918),
       (0.7329, 0.2339), (0.7061, 0.2773), (0.6794, 0.3300), (0.6527, 0.3934), (0.6260, 0.4591),
       (0.5993, 0.5242), (0.5726, 0.5742), (0.5459, 0.6265), (0.5191, 0.6881), (0.4924, 0.7603),
       (0.4657, 0.8302), (0.4390, 0.8880), (0.4123, 0.9420), (0.3856, 0.9783), (0.3588, 0.9949),
       (0.9933, 0.0002), (0.9666, 0.0008), (0.9399, 0.0034), (0.9132, 0.0081), (0.8865, 0.0189),
       (0.8597, 0.0401), (0.8330, 0.0718), (0.8063, 0.1122), (0.7796, 0.1574), (0.7529, 0.2005),
       (0.7262, 0.2436), (0.6995, 0.2909), (0.6727, 0.3448), (0.6460, 0.4094), (0.6193, 0.4771),
       (0.5926, 0.5390), (0.5659, 0.5869), (0.5392, 0.6415), (0.5125, 0.7043), (0.4857, 0.7810),
       (0.4590, 0.8473), (0.4323, 0.9029), (0.4056, 0.9549), (0.3789, 0.9827), (0.3522, 0.9977),
       (0.9866, 0.0003), (0.9599, 0.0010), (0.9332, 0.0040), (0.9065, 0.0100), (0.8798, 0.0239),
       (0.8531, 0.0464), (0.8264, 0.0823), (0.7996, 0.1254), (0.7729, 0.1704), (0.7462, 0.2114),
       (0.7195, 0.2561), (0.6928, 0.3030), (0.6661, 0.3616), (0.6394, 0.4260), (0.6126, 0.4932),
       (0.5859, 0.5501), (0.5592, 0.5992), (0.5325, 0.6544), (0.5058, 0.7218), (0.4791, 0.7992),
       (0.4523, 0.8599), (0.4256, 0.9159), (0.3989, 0.9647), (0.3722, 0.9867), (0.3455, 0.9991),
       (0.9800, 0.0004), (0.9532, 0.0015), (0.9265, 0.0058), (0.8998, 0.0137), (0.8731, 0.0290),
       (0.8464, 0.0517), (0.8197, 0.0906), (0.7930, 0.1353), (0.7662, 0.1823), (0.7395, 0.2232),
       (0.7128, 0.2670), (0.6861, 0.3163), (0.6594, 0.3769), (0.6327, 0.4420), (0.6060, 0.5089),
       (0.5792, 0.5625), (0.5525, 0.6134), (0.5259, 0.6706), (0.4991, 0.7410), (0.4724, 0.8158),
       (0.4457, 0.8758), (0.4190, 0.9293), (0.3922, 0.9721), (0.3655, 0.9905), (0.3388, 1.0000)]


def ldc_levels(k=4):
    pts = sorted(LDC, key=lambda t: t[1])
    per = np.array([p for _, p in pts]); load = np.array([l for l, _ in pts])
    t = np.linspace(0, 1, 4001)
    lt = np.interp(t, per, load)
    out = []
    for i in range(k):
        m = (t >= i / k) & (t <= (i + 1) / k)
        out.append((float(lt[m].mean()), 1.0 / k))
    return out


CLPU = dict(kmax={"residential": 2.0, "commercial": 2.0, "government": 2.0,
                  "industrial": 1.2}, tc_h=1.0, tdec_h=0.5, step_h=10 / 60)
LEVELS = ldc_levels(4)
MODELS = {"S": dict(), "F": dict(clpu=CLPU, load_levels=LEVELS)}


class Problem:
    def __init__(self, model, cache_file=None, clpu=None):
        kw, self.info = R.build("peak")
        fixed = {s: REC for s in self.info["feeder_head"]}
        self.net = Network(**kw, fixed=fixed, s_max_kva=R.SP_CAPACITY_KVA)
        p = default_params()
        p.update(tau_ms=R.SWITCHING, r_rep=R.LINE_REPAIR)
        p.update(MODELS[model])
        if clpu is not None and model == "F":
            p["clpu"] = clpu
        self.ev = Evaluator(self.net, p, mode="loadflow")
        self.cand = [s for s in self.info["main"] if s not in self.info["feeder_head"]] + \
            list(self.info["lateral"])
        self.nsec, self.ntie = len(self.cand), len(self.net.ties)
        self.cache_file, self.cache = cache_file, {}
        if cache_file:
            try:
                self.cache = pickle.load(open(cache_file, "rb"))
            except Exception:
                pass

    def save(self):
        if self.cache_file:
            pickle.dump(self.cache, open(self.cache_file, "wb"))

    def decode(self, g):
        x = [NONE] * len(self.net.sec)
        for i, s in enumerate(self.cand):
            x[s] = g[i]
        return x, list(g[self.nsec:])

    def cost(self, g):
        g = tuple(g)
        if g not in self.cache:
            x, w = self.decode(g)
            o = self.ev.evaluate(x, w)
            self.cache[g] = o["f1"] + o["ECOST"] + (0 if o["rec_ok"] else 1e6)
        return self.cache[g]

    def describe(self, g):
        out = []
        for i, s in enumerate(self.cand):
            if g[i] != NONE:
                out.append(f"{DEVICE_NAME[g[i]]}@s{self.info['sec_no'][s]}")
        for c in range(self.ntie):
            if g[self.nsec + c]:
                out.append(f"tie{c+1}:{'M' if g[self.nsec+c]==1 else 'R'}")
        return out


def random_genome(pr, rng, density=0.15):
    return [rng.choice([MS, RCS, REC]) if rng.random() < density else NONE
            for _ in range(pr.nsec)] + [rng.choice([0, 1, 2]) for _ in range(pr.ntie)]


def mutate(pr, g, rng):
    g = list(g); pm = 1.5 / len(g)
    for i in range(len(g)):
        if rng.random() < pm:
            g[i] = rng.choice([NONE, NONE, MS, RCS, REC]) if i < pr.nsec else rng.choice([0, 1, 2])
    return g


def ga(pr, seed, pop=40, gens=60, seeds=()):
    rng = random.Random(seed)
    P = [list(s) for s in seeds] + [random_genome(pr, rng) for _ in range(pop - len(seeds))]
    for _ in range(gens):
        P.sort(key=pr.cost)
        nxt = [P[0], P[1]]
        while len(nxt) < pop:
            a = min(rng.sample(P, 3), key=pr.cost); b = min(rng.sample(P, 3), key=pr.cost)
            nxt.append(mutate(pr, [a[i] if rng.random() < 0.5 else b[i] for i in range(len(a))], rng))
        P = nxt
    return min(P, key=pr.cost)


def local_search(pr, g):
    g, best, improved = list(g), pr.cost(g), True
    while improved:
        improved = False
        for i in range(len(g)):
            for v in ([NONE, MS, RCS, REC] if i < pr.nsec else [0, 1, 2]):
                if v == g[i]:
                    continue
                h = list(g); h[i] = v
                if pr.cost(h) < best - 1e-6:
                    g, best, improved = h, pr.cost(h), True
    return g


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "S"
    prS, prF = Problem("S", "rbts_cacheS.pkl"), Problem("F", "rbts_cacheF.pkl")
    t0 = time.perf_counter()
    if stage == "S":
        seeds = []
        g0 = [NONE] * prS.nsec + [0, 0]
        for i, s in enumerate(prS.cand):
            if s in prS.info["main"]:
                g0[i] = MS
        seeds.append(g0)
        xs = [local_search(prS, ga(prS, 10 + r, seeds=seeds, gens=60)) for r in range(3)]
        xS = min(xs, key=prS.cost); prS.save()
        json.dump(list(xS), open("rbts_xS.json", "w"))
        print(f"x_S* C_S={prS.cost(xS)/1e3:.2f} {prS.describe(xS)} ({time.perf_counter()-t0:.0f} s)")
    elif stage == "F":
        xS = json.load(open("rbts_xS.json"))
        xF = local_search(prF, ga(prF, 20, pop=30, gens=int(sys.argv[2]) if len(sys.argv) > 2 else 25,
                                  seeds=[xS]))
        prF.save()
        json.dump(list(xF), open("rbts_xF.json", "w"))
        print(f"x_F* C_F={prF.cost(xF)/1e3:.2f} {prF.describe(xF)} ({time.perf_counter()-t0:.0f} s)")
    elif stage == "final":
        xS, xF = json.load(open("rbts_xS.json")), json.load(open("rbts_xF.json"))
        if prS.cost(xF) < prS.cost(xS):
            xS = local_search(prS, xF)
        res = dict(xS=prS.describe(xS), xF=prS.describe(xF), C_S_of_xS=prS.cost(xS),
                   C_S_of_xF=prS.cost(xF), C_F_of_xS=prF.cost(xS), C_F_of_xF=prF.cost(xF),
                   regret=(prF.cost(xS) - prF.cost(xF)) / prF.cost(xF),
                   genes_different=sum(1 for a, b in zip(xS, xF) if a != b), levels=LEVELS)
        prS.save(); prF.save()
        print(json.dumps(res, indent=1)); json.dump(res, open("rbts_go_nogo.json", "w"), indent=1)
