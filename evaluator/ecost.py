"""Customer interruption cost (ECOST) from sector customer damage functions (SCDF).

SCDF_RBTS below is VERIFIED against Table XV ("Cost of interruption in $/kW",
1987 Canadian dollars) of
  R. Billinton et al., "A reliability test system for educational purposes - basic
  data," IEEE Trans. Power Systems 4(3):1238-1244, 1989.
Values are normalized by the customer's annual PEAK demand, so ECOST must multiply
them by peak load. Residential 1 min and 8 h values are extrapolated in the source.

Conversion to 2025 US dollars:
  factor = CPI_Canada(2025) / CPI_Canada(BASE_YEAR) / (CAD per USD, 2025 average)
         = 164.2 / 68.5 / 1.3978 = 1.7149
  CPI, all items, annual average (2002=100): 1987 = 68.5, 1991 = 82.8 (Statistics Canada,
                  table 18-10-0005-01); 2025 = 164.2 (Statistics Canada, The Daily,
                  19 January 2026, "Consumer Price Index: Annual review, 2025").
  CAD per USD, 2025 annual average = 1.3978 (Bank of Canada, annual average exchange rates).

History: until 6 October 2026 this file used CPI 2025 = 158.4 (a secondary figure that turned
out to be too low), i.e. a factor of 1.6543. f2 is proportional to the factor, so every f2
value computed with the old file is the present value divided by 164.2/158.4 = 1.036616.
"""
import numpy as np

DURATIONS_MIN = [1.0, 20.0, 60.0, 240.0, 480.0]

# $/kW of annual peak demand, 1987 Canadian dollars (Billinton et al. 1989, Table XV)
SCDF_RBTS = {
    "residential": [0.001, 0.093, 0.482, 4.914, 15.690],
    "commercial":  [0.381, 2.969, 8.552, 31.317, 83.008],
    "industrial":  [1.625, 3.868, 9.085, 25.163, 55.808],
    "government":  [0.044, 0.369, 1.492, 6.558, 26.040],
}

BASE_YEAR = 1987            # stated in Billinton et al. 1989 (1987 Cdn $ base)
CPI_CANADA = {1987: 68.5, 1991: 82.8, 2025: 164.2}
CAD_PER_USD_2025 = 1.3978


def conversion_factor(base_year=BASE_YEAR):
    return CPI_CANADA[2025] / CPI_CANADA[base_year] / CAD_PER_USD_2025


def damage(sector, r_hours, table=SCDF_RBTS, factor=None):
    """Cost per kW (USD 2025) for an interruption of r hours; linear interpolation,
    linear extrapolation beyond 8 h with the last slope, 0 for r <= 0."""
    if r_hours <= 0:
        return 0.0
    f = conversion_factor() if factor is None else factor
    t, c = DURATIONS_MIN, table[sector]
    m = r_hours * 60.0
    if m <= t[0]:
        val = c[0] * m / t[0]
    elif m >= t[-1]:
        val = c[-1] + (c[-1] - c[-2]) / (t[-1] - t[-2]) * (m - t[-1])
    else:
        val = float(np.interp(m, t, c))
    return f * val
