"""IEEE 69-bus test system (Baran & Wu, 1989), 12.66 kV.

Branch impedances (ohm) and peak loads (kW, kvar) parsed from MATPOWER case69.m v2
(ieee69_raw.py), which cites Baran & Wu, IEEE Trans. Power Del. 4(1):725-734, 1989.
Code bus i = IEEE bus i+1; the feeder breaker sits upstream of bus 1.

Assumptions (same conventions as ieee33_data.py, documented in assumptions.md):
  * length: 1 km per section and per tie (default, as for IEEE 33);
    lengths="impedance" derives |Z| / |z'| with AAAC 150 mm2 (0.22 + j0.35 ohm/km)
  * ampacity: 412 A (AAAC 150 mm2)
  * customers: N_j = max(1, round(P_j / 3 kW)) for loaded buses (about the
    3 kW-per-customer ratio of the IEEE 33 PowerFactory model)
  * sector by bus peak load: P >= 200 kW industrial, 60 <= P < 200 commercial,
    P < 60 residential (written rule; no published allocation exists)
  * tie lines (normally open), as commonly used for this system:
    11-43 and 13-21 (0.5 + j0.5 ohm), 15-46 and 27-65 (1.0 + j1.0 ohm),
    50-59 (2.0 + j2.0 ohm)  -- TODO: cite the source used in the paper
"""
from ieee69_raw import SECTIONS_IEEE, LOADS_IEEE

V_KV = 12.66
Z_PER_KM = complex(0.22, 0.35)
AMPACITY_A = 412.0
KW_PER_CUSTOMER = 3.0
TIES_IEEE = [(11, 43, 0.5, 0.5), (13, 21, 0.5, 0.5), (15, 46, 1.0, 1.0),
             (50, 59, 2.0, 2.0), (27, 65, 1.0, 1.0)]


def sector_of(p_kw):
    if p_kw >= 200:
        return "industrial"
    return "commercial" if p_kw >= 60 else "residential"


def _length_km(r, x, lengths):
    return 1.0 if lengths == "uniform" else abs(complex(r, x)) / abs(Z_PER_KM)


def build(lengths="uniform"):
    n = 69
    sections = [dict(frm=f - 1, to=t - 1, R=r, X=x, length_km=_length_km(r, x, lengths),
                     ampacity_a=AMPACITY_A) for f, t, r, x in SECTIONS_IEEE]
    ties = [dict(a=f - 1, b=t - 1, R=r, X=x, length_km=_length_km(r, x, lengths),
                 ampacity_a=AMPACITY_A) for f, t, r, x in TIES_IEEE]
    P = [LOADS_IEEE[i + 1][0] for i in range(n)]
    Q = [LOADS_IEEE[i + 1][1] for i in range(n)]
    customers = [max(1, round(p / KW_PER_CUSTOMER)) if p > 0 else 0 for p in P]
    sectors = [sector_of(p) if p > 0 else None for p in P]
    return dict(n_bus=n, sections=sections, ties=ties, P_kw=P, Q_kvar=Q,
                customers=customers, root=0, V_kv=V_KV, sectors=sectors)
