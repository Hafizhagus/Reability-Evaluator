"""Step 3: tests with hand-calculated expectations + IEEE 33 load-flow validation."""
import random
import numpy as np

from network import Network, bfs_loadflow, normal_edges, NONE, MS, RCS, REC
from reliability import Evaluator, default_params
import ieee33_data

TOL = 1e-9
H5 = 5.0 / 60.0


def toy(main_len, lateral=None, ties=()):
    """Toy feeder: buses 0..n on a main line, optional lateral (frm, to).
    Every section 1 km, every load bus 100 kW / 10 customers."""
    secs = [dict(frm=i, to=i + 1, R=0.1, X=0.1, length_km=1.0, ampacity_a=1e6)
            for i in range(main_len)]
    n = main_len + 1
    if lateral:
        secs.append(dict(frm=lateral[0], to=n, R=0.1, X=0.1, length_km=1.0,
                         ampacity_a=1e6))
        n += 1
    tie_list = [dict(a=a, b=b, R=0.1, X=0.1, length_km=1.0, ampacity_a=1e6)
                for a, b in ties]
    P = [0] + [100] * (n - 1)
    return Network(n, secs, tie_list, P, [0] * n, [0] + [10] * (n - 1), V_kv=20.0)


def params():
    p = default_params()
    p.update(lambda_per_km=0.1, r_rep=3.0, tau_ms=1.0, tau_rcs=H5)
    return p


def close(a, b, msg):
    assert abs(a - b) < 1e-9, f"{msg}: got {a}, expected {b}"


def test_no_device():
    net = toy(3, lateral=(1, None))
    ev = Evaluator(net, params())
    out = ev.evaluate([NONE] * 4, [])
    close(out["SAIFI"], 0.4, "no-device SAIFI = total lambda")
    close(out["SAIDI"], 1.2, "no-device SAIDI = total lambda x r_rep")


def test_one_rcs():
    net = toy(3, lateral=(1, None))
    ev = Evaluator(net, params())
    out = ev.evaluate([NONE, RCS, NONE, NONE], [])
    u_up = 0.1 * 3 + 0.1 * H5 + 0.1 * H5 + 0.1 * 3
    close(out["SAIDI"], (2 * u_up + 2 * 1.2) / 4, "one RCS SAIDI")
    close(out["SAIFI"], 0.4, "sectionalizing switch does not change SAIFI")


def test_tie_and_timing():
    # RCS at s1, MS at s2, tie 3-4 remote
    net = toy(3, lateral=(1, None), ties=[(3, 4)])
    ev = Evaluator(net, params())
    out = ev.evaluate([NONE, RCS, MS, NONE], [2])
    u14 = 0.1 * 3 + 0.1 * H5 + 0.1 * H5 + 0.1 * 3
    u2 = 0.1 * 3 + 0.1 * 3 + 0.1 * 1.0 + 0.1 * 3
    u3 = 0.1 * 3 + 0.1 * 1.0 + 0.1 * 3 + 0.1 * 3
    close(out["U_j"][1], u14, "bus1"); close(out["U_j"][4], u14, "bus4")
    close(out["U_j"][2], u2, "bus2"); close(out["U_j"][3], u3, "bus3")


def test_far_terminal_energization():
    # MS at s1, RCS at s2, remote tie 3-4: far end (bus 4) is back only after 1 h
    net = toy(3, lateral=(1, None), ties=[(3, 4)])
    ev = Evaluator(net, params())
    B_out = ev.evaluate([NONE, MS, RCS, NONE], [2])
    # fault on s1: bus 3 restored via tie at max(5 min, 5 min, 1 h) = 1 h
    u3 = 0.1 * 3 + 0.1 * 1.0 + 0.1 * 3 + 0.1 * 3
    close(B_out["U_j"][3], u3, "far-terminal energization dominates")


def test_two_stage():
    # main 0-1-2-3-4, RCS at s1, MS at s3; fault on s3
    net = toy(4)
    ev = Evaluator(net, params())
    x = [NONE, RCS, NONE, MS]
    r, _ = ev.fault_durations(x, [], 3)
    close(r[1], H5, "bus 1 restored by the upstream RCS (stage 1)")
    close(r[2], 1.0, "bus 2 restored by the MS (stage 2)")
    close(r[3], 1.0, "bus 3 restored by the MS (stage 2)")
    close(r[4], 3.0, "bus 4 in the faulted zone")


def test_recloser():
    net = toy(3, lateral=(1, None))
    ev = Evaluator(net, params())
    out = ev.evaluate([NONE, REC, NONE, NONE], [])
    close(out["SAIFI"], (0.2 * 2 + 0.4 * 2) / 4, "recloser SAIFI")
    close(out["SAIDI"], (0.6 * 2 + 1.2 * 2) / 4, "recloser SAIDI")


def test_ieee33_loadflow():
    net = Network(**ieee33_data.build())
    res = bfs_loadflow(net, range(net.n), normal_edges(net))
    assert res["converged"] and res["radial"]
    vmin, bus = np.nanmin(res["V"]), int(np.nanargmin(res["V"])) + 1
    print(f"  IEEE 33 base: loss = {res['loss_kw']:.2f} kW, Vmin = {vmin:.4f} pu at bus {bus}")
    assert abs(res["loss_kw"] - 202.7) < 0.5 and abs(vmin - 0.9131) < 5e-4 and bus == 18


def test_ieee33_base_and_head():
    net = Network(**ieee33_data.build())
    ev = Evaluator(net)
    base = ev.evaluate([NONE] * 32, [0] * 5)
    lam_tot = ev.lam.sum()
    close(base["SAIFI"], lam_tot, "IEEE33 no-device SAIFI")
    close(base["SAIDI"], lam_tot * 3.0, "IEEE33 no-device SAIDI")
    x = [NONE] * 32; x[0] = RCS                      # device on section 1-2 only
    head = ev.evaluate(x, [0] * 5)
    close(head["SAIDI"], base["SAIDI"], "device at section 1-2 gives no benefit")


def test_monotonic(mode, trials=150, seed=1):
    net = Network(**ieee33_data.build())
    ev = Evaluator(net, mode=mode)
    rng, bad = random.Random(seed), 0
    for _ in range(trials):
        x = [rng.choice([NONE, NONE, NONE, MS, RCS, REC]) for _ in range(32)]
        w = [rng.choice([0, 0, 1, 2]) for _ in range(5)]
        empty = [s for s in range(32) if x[s] == NONE]
        if not empty:
            continue
        a = ev.evaluate(x, w)
        x2 = list(x); x2[rng.choice(empty)] = rng.choice([MS, RCS])
        b = ev.evaluate(x2, w)
        if b["SAIDI"] > a["SAIDI"] + 1e-9 or b["ENS"] > a["ENS"] + 1e-9:
            bad += 1
    print(f"  monotonic ({mode}): {bad} violations in {trials} trials")
    return bad


if __name__ == "__main__":
    for t in [test_no_device, test_one_rcs, test_tie_and_timing,
              test_far_terminal_energization, test_two_stage, test_recloser,
              test_ieee33_loadflow, test_ieee33_base_and_head]:
        t()
        print(f"PASS {t.__name__}")
    v1 = test_monotonic("connectivity")
    v2 = test_monotonic("loadflow", trials=60)
    print("PASS monotonic" if v1 == 0 else "CHECK monotonic (connectivity)")
