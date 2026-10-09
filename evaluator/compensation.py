"""f3: expected annual TMP compensation (MEMR Reg. 2/2025, Art. 6A).

Compensation per customer and month = multiplier(excess hours above threshold) x base,
where base = "biaya beban / rekening minimum". Minimum bill modelled as
    base = 40 h x contracted power (kVA) x tariff (Rp/kWh)
with contracted power = peak demand per customer / PF_CONTRACT.

Tariffs (Rp/kWh, non-subsidised, MEMR Reg. 7/2024; unchanged Q1-Q3 2026):
  R-1 1300-2200 VA 1444.70 | R-2 3.5-5.5 kVA 1699.53 | R-3 >=6.6 kVA 1699.53
  B-2 6.6-200 kVA 1444.70  | B-3/TM >200 kVA 1114.74
  I-3/TM >200 kVA 1114.74  | I-4/TT >30 MVA 996.74
  P-1 6.6-200 kVA 1699.53  | P-2/TM >200 kVA 1522.88
Customers below the non-subsidised band of their sector (small business B-1, small
industry I-1/I-2) are priced at the lowest non-subsidised tariff of that sector
(assumption, see assumptions.md).
Exchange rate: 16,959.32 Rp/USD (macro parameter used by MEMR for the Q3-2026 tariff
setting, Feb-Apr 2026 average).
"""
import numpy as np

BRACKETS_H = [0, 2, 4, 8, 16, 40]
MULT = [0.50, 0.75, 1.00, 2.00, 3.00, 5.00]
MIN_HOURS = 40.0
PF_CONTRACT = 0.85
RP_PER_USD = 16959.32


def golongan(sector, kva):
    """(tariff class, Rp/kWh) for a customer of the given sector and contracted kVA."""
    if sector == "residential":
        if kva <= 2.2:
            return "R-1", 1444.70
        return ("R-2", 1699.53) if kva <= 5.5 else ("R-3", 1699.53)
    if sector == "commercial":
        return ("B-3/TM", 1114.74) if kva > 200 else ("B-2", 1444.70)
    if sector == "industrial":
        return ("I-4/TT", 996.74) if kva > 30000 else ("I-3/TM", 1114.74)
    if sector == "government":
        return ("P-2/TM", 1522.88) if kva > 200 else ("P-1", 1699.53)
    raise ValueError(sector)


def monthly_base_usd(net, loads, rp_per_usd=RP_PER_USD, tariff_scale=None, contract_kva=None):
    """Minimum-bill base per customer (USD/month) and tariff class for each load point.
    tariff_scale: optional {sector: factor} for the sector-ratio sensitivity."""
    base, cls = [], []
    for j in loads:
        kva = contract_kva if contract_kva else net.P[j] / net.N[j] / PF_CONTRACT
        g, rp = golongan(net.sector[j], kva)
        f = (tariff_scale or {}).get(net.sector[j], 1.0)
        base.append(MIN_HOURS * kva * rp * f / rp_per_usd)
        cls.append(g)
    return np.array(base), cls


def expected_multiplier(sev, thr_h):
    """Monthly compound Poisson on a 1-minute grid (Panjer recursion), all load points
    at once. sev[i] = {duration_min: annual rate}. Severities are capped at the top
    bracket, so the pmf below the cap is exact.
    Returns E[multiplier], E[monthly duration] (h), P(D > threshold)."""
    cap = int(round((thr_h + BRACKETS_H[-1]) * 60)) + 1
    n = len(sev)
    supp = sorted({min(d, cap) for s in sev for d in s})
    if not supp:
        z = np.zeros(n)
        return z, z, z
    lam = np.zeros(n); ED = np.zeros(n); f = np.zeros((n, len(supp)))
    for i, s in enumerate(sev):
        tot = sum(s.values())
        lam[i] = tot / 12.0
        ED[i] = sum(d * r for d, r in s.items()) / 12.0 / 60.0
        for d, r in s.items():
            if tot > 0:
                f[i, supp.index(min(d, cap))] += r / tot
    g = np.zeros((n, cap)); g[:, 0] = np.exp(-lam)
    coef = lam[:, None] * np.array(supp)[None, :] * f
    for s_ in range(1, cap):
        acc = np.zeros(n)
        for k, y in enumerate(supp):
            if y <= s_:
                acc += coef[:, k] * g[:, s_ - y]
        g[:, s_] = acc / s_
    cdf = np.cumsum(g, axis=1)
    Em = np.zeros(n); prev = 0.0
    for lo, v in zip(BRACKETS_H, MULT):
        Em += (v - prev) * (1.0 - cdf[:, int(round((thr_h + lo) * 60))])
        prev = v
    return Em, ED, 1.0 - cdf[:, int(round(thr_h * 60))]


def f3(net, loads, sev, thr_h=1.0, base=None):
    """Expected annual compensation (USD/yr) and share of customers whose expected
    monthly interruption duration exceeds the threshold."""
    if base is None:
        base, _ = monthly_base_usd(net, loads)
    N = np.array([net.N[j] for j in loads])
    Em, ED, _ = expected_multiplier(sev, thr_h)
    return 12.0 * float((N * base * Em).sum()), float(N[ED > thr_h].sum() / N.sum())
