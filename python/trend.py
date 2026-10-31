"""Trend statistics of the annual mapped area (Section 4.5 of the paper).

Mann-Kendall test with continuity correction (no ties in areas), Kendall's tau, Theil-Sen slope with
the normal-approximation interval of Gilbert (1987): the M1-th and (M2 + 1)-th of the N ordered
pairwise slopes, M1 = (N - C)/2 and M2 = (N + C)/2 rounded inward, C = 1.96 sqrt(Var S). Period means and the series maximum.
Usage: python trend.py figure_data/fig6_annual_series.csv
"""
import sys, math
import numpy as np
import pandas as pd
from scipy.stats import norm

TRANSITION_YEARS = (1999, 2000, 2012, 2013, 2014, 2022, 2023)   # Table S19


def mk_sen(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(y)
    s = sum(np.sign(y[j] - y[i]) for i in range(n) for j in range(i + 1, n))
    var = n * (n - 1) * (2 * n + 5) / 18
    z = (s - 1) / math.sqrt(var) if s > 0 else (s + 1) / math.sqrt(var) if s < 0 else 0.0
    p = 2 * norm.sf(abs(z))
    tau = s / (n * (n - 1) / 2)
    slopes = np.sort([(y[j] - y[i]) / (x[j] - x[i]) for i in range(n) for j in range(i + 1, n)])
    N = len(slopes)
    c = 1.96 * math.sqrt(var)
    lo, hi = slopes[math.ceil((N - c) / 2) - 1], slopes[math.floor((N + c) / 2 + 1) - 1]
    return dict(n=n, tau=tau, z=z, p=p, sen=float(np.median(slopes)), lo=float(lo), hi=float(hi))


def show(label, r):
    print(f"{label:52s} n={r['n']:2d}  tau={r['tau']:+.3f}  Z={r['z']:+.2f}  p={r['p']:.4f}  "
          f"Sen={r['sen']:+.1f} km2/yr  95% CI {r['lo']:+.1f} to {r['hi']:+.1f}")


if __name__ == '__main__':
    a = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else 'figure_data/fig6_annual_series.csv')
    full = a[a.coverage_pct >= 99]
    core = full[full.year >= 1998]
    show('1998-2025, coverage >= 99 % (primary)', mk_sen(core.year, core.primary_km2))
    show('1989-2025, incomplete early years included', mk_sen(a.year, a.primary_km2))
    nt = core[~core.year.isin(TRANSITION_YEARS)]
    show('1998-2025 without the seven sensor-transition years', mk_sen(nt.year, nt.primary_km2))
    for y0, y1 in ((1998, 2009), (2010, 2019), (2020, 2025)):
        m = a[(a.year >= y0) & (a.year <= y1)].primary_km2.mean()
        print(f'mean mapped area {y0}-{y1}: {m:.0f} km2')
    i = a.primary_km2.idxmax()
    print(f"maximum: {a.primary_km2[i]:.0f} km2 in {int(a.year[i])}")
