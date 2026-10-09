"""Regresi otomatis: semua hasil yang sudah tervalidasi harus tetap sama."""
import numpy as np
from network import Network, NONE, MS, RCS, REC, bfs_loadflow, normal_edges
from reliability import Evaluator, default_params
import ieee33_data, rbts2_data as R
ok = True
def check(name, got, exp, tol):
    global ok
    good = abs(got - exp) <= tol
    ok &= good
    print(f"{'OK ' if good else 'GAGAL'} {name}: {got:.4f} (acuan {exp:.4f})")
data = ieee33_data.build(); net = Network(**data)
S = {(d['frm']+1, d['to']+1): i for i, d in enumerate(net.sec)}; T = {(t['a']+1, t['b']+1): i for i, t in enumerate(net.ties)}
def dev(d): 
    x = [NONE]*32
    for q, t in d.items(): x[S[q]] = t
    return x
def tie(d):
    w = [0]*5
    for q, t in d.items(): w[T[q]] = t
    return w
ev = Evaluator(net, default_params(), mode="connectivity")
V = [("V0", dev({}), tie({}), 19.2000, 71.328), ("V1", dev({(2,3):RCS}), tie({}), 17.2076, 64.083),
     ("V2", dev({(2,3):MS}), tie({}), 17.8338, 66.360), ("V3", dev({(2,3):RCS,(6,7):MS}), tie({}), 14.4022, 53.619),
     ("V4", dev({(2,3):REC}), tie({}), 17.1507, 63.876), ("V5", dev({(15,16):RCS,(16,17):RCS}), tie({}), 17.5308, 65.124),
     ("V6", dev({(15,16):RCS,(16,17):RCS}), tie({(18,33):1}), 17.5145, 65.064),
     ("V7", dev({(15,16):RCS,(16,17):RCS}), tie({(18,33):2}), 17.5070, 65.037),
     ("V8", dev({(2,3):RCS,(8,9):RCS,(2,19):RCS,(3,23):RCS,(27,28):RCS}), tie({}), 7.6716, 28.543)]
for name, x, w, sd, ens in V:
    o = ev.evaluate(x, w); check(f"{name} SAIDI", o["SAIDI"], sd, 6e-5); check(f"{name} ENS", o["ENS"], ens, 1.5e-3)
lf = bfs_loadflow(net, range(33), normal_edges(net))
check("L1 rugi kW", lf["loss_kw"], 202.68, 0.02); check("L1 Vmin", float(np.nanmin(lf["V"])), 0.9131, 1e-4)
kw, info = R.build(); rnet = Network(**kw, fixed={s: REC for s in info["feeder_head"]})
p = default_params(); p.update(tau_ms=R.SWITCHING, r_rep=R.LINE_REPAIR); rev = Evaluator(rnet, p, mode="connectivity")
bus_all = [info["lp_bus"][lp] for lp in range(1, 23)]
def saidi(x, w):
    o = rev.evaluate(x, w); N = rnet.N[bus_all]; return float((o["U_j"][bus_all]*N).sum()/N.sum())
xB = [NONE]*len(rnet.sec); xF = list(xB)
for s in info["main"]:
    if s not in info["feeder_head"]: xF[s] = MS
check("RBTS kasus B SAIDI", saidi(xB, [0,0]), 22.50, 0.006)
check("RBTS kasus F SAIDI", saidi(xF, [0,0]), 9.93, 0.006)
check("RBTS kasus D SAIDI", saidi(xF, [1,1]), 6.745, 0.01)
# IEEE 69: load flow base case (Baran & Wu 1989) dan cek tangan konektivitas
import ieee69_data
n69 = Network(**ieee69_data.build())
lf69 = bfs_loadflow(n69, range(n69.n), normal_edges(n69))
check("IEEE69 rugi kW", lf69["loss_kw"], 224.99, 0.05)
check("IEEE69 Vmin", float(np.nanmin(lf69["V"])), 0.9092, 1e-4)
ev69 = Evaluator(n69, default_params(), mode="connectivity")
o69 = ev69.evaluate([NONE]*68, [0]*5)
check("IEEE69 base SAIFI (=0,2 x 68)", o69["SAIFI"], 13.6, 1e-9)
check("IEEE69 base SAIDI (=13,6 x 3)", o69["SAIDI"], 40.8, 1e-9)
# satu RCS di seksi 3-4 (hitung tangan): sisi hulu = bus 1-3 + lateral 28-35 dan 36-46
# (2 + 8 + 11 = 21 seksi); sisi hilir = 47 seksi. Pelanggan hulu: gangguan hilir pulih
# 5 menit, gangguan hulu 3 jam -> U_up = 0,2*(21*3 + 47*5/60); pelanggan hilir: 40,8.
up_bus = [b - 1 for b in [2, 3] + list(range(28, 36)) + list(range(36, 47))]
N_up = sum(n69.N[b] for b in up_bus); N_tot = float(n69.N.sum())
U_up = 0.2 * (21 * 3 + 47 * 5 / 60)
saidi_hand = (N_up * U_up + (N_tot - N_up) * 40.8) / N_tot
x = [NONE]*68; x[2] = RCS
check("IEEE69 RCS di 3-4 vs hitung tangan", ev69.evaluate(x, [0]*5)["SAIDI"], saidi_hand, 1e-9)
# Bawean: load flow terhadap hasil PowerFactory di ekspor (V akar PF 1,002485 pu)
import pandas as pd, bawean_data as BW
bkw, binfo = BW.build(); bnet = Network(**bkw)
blf = bfs_loadflow(bnet, range(bnet.n), normal_edges(bnet))
bn = binfo["raw"]["nodes"]; bn["u_pu"] = pd.to_numeric(bn["u_pu"], errors="coerce")
pfv = dict(zip(bn.terminal, bn.u_pu)); v0 = pfv["Bus2"]
dmax = max(abs((v0 - pfv[t]) - (1.0 - blf["V"][i])) for i, t in enumerate(binfo["order"]) if t != "ROOT" and t in pfv)
check("Bawean rugi saluran kW (PF dikoreksi ke 1,0 pu)", blf["loss_kw"], 33.905 * v0 ** 2, 0.1)
check("Bawean selisih jatuh tegangan maks (pu)", dmax, 0.0, 2e-4)
print("SEMUA REGRESI LOLOS" if ok else "ADA REGRESI GAGAL")
