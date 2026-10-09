"""IEEE 33-bus test system (Baran & Wu, 1989), 12.66 kV.

Bus numbering in code is 0-based: code bus i = IEEE bus i+1 (bus 0 = IEEE bus 1 = GI).
The feeder breaker (CB) sits upstream of IEEE bus 1, so every section is a candidate.

The original IEEE 33 data give only R, X and loads. Reliability studies add their own
assumptions for line length and customers. This file uses:

  * Line length: DEFAULT is a uniform 1 km per section and per tie, matching the
    PowerFactory model of this study (base case SAIFI 6.400, SAIDI 19.200 h/cust.yr
    with the SPLN 59:1985 parameters). The IEEE 33 data give no lengths.
    build(lengths="impedance") instead derives each length from |Z| assuming one
    conductor type, AAAC 150 mm2 with z' = 0.22 + j0.35 ohm/km (total 66.85 km);
    use it for the sensitivity analysis.
  * Ampacity: 412 A for every branch (AAAC 150 mm2 summer rating, ESB Networks
    "Line impedance data for 38 kV"). Independent of the length assumption.
  * Customers per bus: default from the PowerFactory model of this study
    (build(customers="powerfactory")); a literature allocation is available with
    build(customers="literature"). Customer sector per bus comes from the
    literature source either way.
"""

V_KV = 12.66

# (from_ieee, to_ieee, R_ohm, X_ohm)  -- 32 normally closed sections
SECTIONS_IEEE = [
    (1, 2, 0.0922, 0.0470), (2, 3, 0.4930, 0.2511), (3, 4, 0.3660, 0.1864),
    (4, 5, 0.3811, 0.1941), (5, 6, 0.8190, 0.7070), (6, 7, 0.1872, 0.6188),
    (7, 8, 0.7114, 0.2351), (8, 9, 1.0300, 0.7400), (9, 10, 1.0440, 0.7400),
    (10, 11, 0.1966, 0.0650), (11, 12, 0.3744, 0.1238), (12, 13, 1.4680, 1.1550),
    (13, 14, 0.5416, 0.7129), (14, 15, 0.5910, 0.5260), (15, 16, 0.7463, 0.5450),
    (16, 17, 1.2890, 1.7210), (17, 18, 0.7320, 0.5740), (2, 19, 0.1640, 0.1565),
    (19, 20, 1.5042, 1.3554), (20, 21, 0.4095, 0.4784), (21, 22, 0.7089, 0.9373),
    (3, 23, 0.4512, 0.3083), (23, 24, 0.8980, 0.7091), (24, 25, 0.8960, 0.7011),
    (6, 26, 0.2030, 0.1034), (26, 27, 0.2842, 0.1447), (27, 28, 1.0590, 0.9337),
    (28, 29, 0.8042, 0.7006), (29, 30, 0.5075, 0.2585), (30, 31, 0.9744, 0.9630),
    (31, 32, 0.3105, 0.3619), (32, 33, 0.3410, 0.5302),
]

# Candidate tie lines (from_ieee, to_ieee, R_ohm, X_ohm)
TIES_IEEE = [
    (8, 21, 2.0, 2.0), (9, 15, 2.0, 2.0), (12, 22, 2.0, 2.0),
    (18, 33, 0.5, 0.5), (25, 29, 0.5, 0.5),
]

# Peak loads at IEEE bus 1..33 (kW, kvar)
LOADS_IEEE = {
    1: (0, 0), 2: (100, 60), 3: (90, 40), 4: (120, 80), 5: (60, 30), 6: (60, 20),
    7: (200, 100), 8: (200, 100), 9: (60, 20), 10: (60, 20), 11: (45, 30),
    12: (60, 35), 13: (60, 35), 14: (120, 80), 15: (60, 10), 16: (60, 20),
    17: (60, 20), 18: (90, 40), 19: (90, 40), 20: (90, 40), 21: (90, 40),
    22: (90, 40), 23: (90, 50), 24: (420, 200), 25: (420, 200), 26: (60, 25),
    27: (60, 25), 28: (60, 20), 29: (120, 70), 30: (200, 600), 31: (150, 70),
    32: (210, 100), 33: (60, 40),
}

import math

Z_PER_KM = complex(0.22, 0.35)   # AAAC 150 mm2, only for lengths="impedance"
AMPACITY_A = 412.0               # AAAC 150 mm2 summer rating (ESB Networks)

# Customers per load point: exported from the user's PowerFactory IEEE 33 model
# ("Number of connected customers" per Load_i). Note: these are proportional to the
# bus load, about 3 kW per customer (ratio 2.81-3.00), so SAIDI and ENS are collinear
# under this allocation -- see CUSTOMERS_LIT for an alternative.
CUSTOMERS_IEEE = {2: 34, 3: 31, 4: 40, 5: 20, 6: 20, 7: 67, 8: 67, 9: 20, 10: 20,
                  11: 16, 12: 20, 13: 20, 14: 40, 15: 20, 16: 20, 17: 20, 18: 31,
                  19: 31, 20: 31, 21: 31, 22: 31, 23: 31, 24: 140, 25: 140, 26: 20,
                  27: 20, 28: 20, 29: 40, 30: 67, 31: 51, 32: 70, 33: 20}

# Alternative (literature) allocation and the customer sector per bus:
# Jaleel & Abd, Int. J. Intell. Eng. Syst. 14(5), 2021, Appendix Table 4
# (attributed there to Kumar et al., Energies 13(21):5631, 2020).
# Bus 10 is missing in that table -> assumed 10 (same as buses 6-9).
_GROUPS = [((2, 3, 4, 5), 148, "industrial"), ((6, 7, 8, 9), 10, "commercial"),
           ((10,), 10, "commercial"),
           ((11, 12), 132, "commercial"), ((13, 14, 15), 110, "residential"),
           ((16,), 2, "residential"), ((17, 18, 19, 20), 118, "residential"),
           ((21, 22, 23, 24, 25, 26), 126, "residential"),
           ((27, 28, 29, 30, 31), 108, "residential"), ((32, 33), 58, "residential")]
CUSTOMERS_LIT = {b: n for buses, n, _ in _GROUPS for b in buses}
SECTOR_IEEE = {b: t for buses, _, t in _GROUPS for b in buses}


def _length_km(r, x, lengths):
    """1 km per branch by default; 'impedance' scales |Z| by the assumed conductor."""
    return 1.0 if lengths == "uniform" else abs(complex(r, x)) / abs(Z_PER_KM)


def build(customers="powerfactory", lengths="uniform"):
    """Return a dict consumed by network.Network."""
    sections = [dict(frm=f - 1, to=t - 1, R=r, X=x, length_km=_length_km(r, x, lengths),
                     ampacity_a=AMPACITY_A) for f, t, r, x in SECTIONS_IEEE]
    ties = [dict(a=f - 1, b=t - 1, R=r, X=x, length_km=_length_km(r, x, lengths),
                 ampacity_a=AMPACITY_A) for f, t, r, x in TIES_IEEE]
    P = [LOADS_IEEE[i + 1][0] for i in range(33)]
    Q = [LOADS_IEEE[i + 1][1] for i in range(33)]
    table = CUSTOMERS_IEEE if customers == "powerfactory" else CUSTOMERS_LIT
    customers = [table.get(i + 1, 0) for i in range(33)]
    sectors = [SECTOR_IEEE.get(i + 1) for i in range(33)]
    return dict(n_bus=33, sections=sections, ties=ties, P_kw=P, Q_kvar=Q,
                customers=customers, root=0, V_kv=V_KV, sectors=sectors)
