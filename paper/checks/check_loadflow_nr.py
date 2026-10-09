"""Cek independen load flow: backward-forward sweep evaluator vs Newton-Raphson (pandapower)
untuk base case IEEE 33, IEEE 69, RBTS Bus 2 dan Bawean (data yang sama, pemecah berbeda).
  pip install pandapower ;  python check_loadflow_nr.py"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "evaluator"))
import pandapower as pp
from network import bfs_loadflow, normal_edges
from systems import System
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
for name in ("ieee33", "ieee69", "rbts2", "bawean"):
    sy = System(name); net = sy.net
    lf = bfs_loadflow(net, range(net.n), normal_edges(net))
    p = pp.create_empty_network(sn_mva=10.0)
    b = [pp.create_bus(p, vn_kv=net.V_kv) for _ in range(net.n)]
    pp.create_ext_grid(p, b[net.root], vm_pu=1.0)
    for (u, v, r, x, amp) in normal_edges(net):
        if r == 0 and x == 0:                      # zero-impedance branch (LV transformer in RBTS)
            pp.create_switch(p, b[u], b[v], et="b", closed=True)
            continue
        pp.create_line_from_parameters(p, b[u], b[v], length_km=1.0, r_ohm_per_km=r, x_ohm_per_km=x,
                                       c_nf_per_km=0.0, max_i_ka=amp / 1e3 if np.isfinite(amp) else 10.0)
    for j in range(net.n):
        if net.P[j] or net.Q[j]:
            pp.create_load(p, b[j], p_mw=net.P[j] / 1e3, q_mvar=net.Q[j] / 1e3)
    pp.runpp(p, algorithm="nr", tolerance_mva=1e-9, max_iteration=100, numba=False)
    v = p.res_bus.vm_pu.to_numpy(); loss = float(p.res_line.pl_mw.sum() * 1e3)
    dv = float(np.nanmax(np.abs(v - lf["V"])))
    out[name] = dict(loss_bfs_kw=float(lf["loss_kw"]), loss_nr_kw=loss, vmin_bfs=float(np.nanmin(lf["V"])),
                     vmin_nr=float(v.min()), max_dv_pu=dv)
    print(f"{name:7s}: rugi BFS {lf['loss_kw']:.3f} kW vs NR {loss:.3f} kW | Vmin BFS {np.nanmin(lf['V']):.5f} vs NR {v.min():.5f} pu "
          f"| selisih tegangan maks {dv:.2e} pu")
json.dump(out, open(os.path.join(HERE, "check_loadflow_nr.json"), "w"), indent=1)
