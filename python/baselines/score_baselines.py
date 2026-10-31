"""Tune and score the LandTrendr and CCDC baselines against the anniversary rule.

Tuning: one magnitude threshold per variant, chosen on the 400 calibration points (2025) by maximising
area-weighted F1, the criterion used for the rule's thresholds. Calibration strata: preliminary-map jhum
(100 points, 723.9 km2) and background (300 points, rest of the 21,174.4 km2 frame).
Scoring: the 520-point 2025 sample (five strata) and the 160-point 2010 sample (three strata) with the
stratified estimators of the paper (Olofsson et al., 2014; Stehman, 2014); 4,000-resample stratified
bootstrap for intervals and for paired differences with the rule; exact McNemar on discordant points.
"""
import os, json, csv, math, random, collections, sys
import numpy as np

DD = os.environ.get('JHUM_DATA', 'jhum-mapper-data') + '/'
EE = ''
SEED, NBOOT = 20261026, 4000

cal = {r['point_id']: r for r in csv.DictReader(open(DD + 'calibration_points_400_2025.csv'))}
v25 = {r['point_id']: r for r in csv.DictReader(open(DD + 'validation_points_520_2025.csv'))}
v10 = {r['point_id']: r for r in csv.DictReader(open(DD + 'validation_points_160_2010.csv'))}
A25 = {int(r['stratum']): float(r['area_km2']) for r in csv.DictReader(open(DD + 'stratum_areas_2025.csv'))}
a10 = next(csv.DictReader(open(DD + 'stratum_areas_2010.csv')))
A10 = {1: float(a10['area_stratum1_jhum_km2']), 2: float(a10['area_stratum2_near_km2']), 3: float(a10['area_stratum3_far_km2'])}
ACAL = {1: 723.9, 0: 21174.4 - 723.9}


# ------------------------------------------------------------------ detectors
def lt_loss(arr, Y, mode='segment'):
    """Loss magnitude (NBR x 1000) attributed to clearing year Y by a LandTrendr fit, or 0."""
    if not arr:
        return None
    yrs, src, fit, vtx = arr
    if Y not in yrs:
        return None
    if mode == 'segment':                      # getChangeMap logic: loss segment with yod == Y
        vy = [y for y, v in zip(yrs, vtx) if v == 1]
        for s, e in zip(vy, vy[1:]):
            if s + 1 == Y:
                fs, fe = fit[yrs.index(s)], fit[yrs.index(e)]
                return max(fe - fs, 0.0)        # -NBR x1000 increases with loss
        return 0.0
    i = yrs.index(Y)                           # annual fitted difference
    return max(fit[i] - fit[i - 1], 0.0) if i > 0 else 0.0


def ccdc_loss(rec, Y):
    """Largest confirmed NBR drop (positive number) with a break inside the Y clearing season, or 0."""
    if rec is None or rec.get('tBreak') is None:
        return None
    lo, hi = (Y - 1) + 9 / 12, Y + 9 / 12          # 1 Oct (Y-1) to 1 Oct (Y)
    best = 0.0
    for tb, cp, mg in zip(rec['tBreak'], rec['changeProb'], rec['NBR_magnitude']):
        if tb and lo <= tb < hi and cp >= 1 and mg < 0:
            best = max(best, -mg)
    return best


# ------------------------------------------------------------------ estimators
def estimate(rows, areas):
    """rows: list of (stratum, map, ref). Returns OA, UA, PA, F1, mapped and adjusted area."""
    B = sum(areas.values()); by = collections.defaultdict(list)
    for s, m, y in rows:
        by[s].append((m, y))
    tp = fp = fn = 0.0
    for s, v in by.items():
        W, n = areas[s] / B, len(v)
        for m, y in v:
            if m and y: tp += W / n
            elif m: fp += W / n
            elif y: fn += W / n
    ua = tp / (tp + fp) if tp + fp else float('nan'); pa = tp / (tp + fn) if tp + fn else float('nan')
    f1 = 2 * ua * pa / (ua + pa) if ua == ua and pa == pa and ua + pa else float('nan')
    return dict(OA=1 - fp - fn, UA=ua, PA=pa, F1=f1, mapped=B * (tp + fp), adj=B * (tp + fn))


def boot_paired(rowsA, rowsB, areas, keys=('UA', 'PA', 'F1', 'OA')):
    """rowsA/rowsB aligned point by point. Bootstrap CIs for B and for A - B."""
    random.seed(SEED)
    by = collections.defaultdict(list)
    for (s, ma, y), (_, mb, _) in zip(rowsA, rowsB):
        by[s].append((ma, mb, y))
    eb, ed = collections.defaultdict(list), collections.defaultdict(list)
    for _ in range(NBOOT):
        ra, rb = [], []
        for s, v in by.items():
            for _ in range(len(v)):
                ma, mb, y = v[random.randrange(len(v))]
                ra.append((s, ma, y)); rb.append((s, mb, y))
        a, b = estimate(ra, areas), estimate(rb, areas)
        for k in keys:
            if b[k] == b[k]: eb[k].append(b[k])
            if a[k] == a[k] and b[k] == b[k]: ed[k].append(a[k] - b[k])
    q = lambda v: (np.percentile(v, 2.5), np.percentile(v, 97.5)) if v else (float('nan'),) * 2
    return {k: q(eb[k]) for k in keys}, {k: q(ed[k]) for k in keys}


def mcnemar(rowsA, rowsB):
    b = sum(1 for (_, ma, y), (_, mb, _) in zip(rowsA, rowsB) if ma == y and mb != y)
    c = sum(1 for (_, ma, y), (_, mb, _) in zip(rowsA, rowsB) if mb == y and ma != y)
    n = b + c
    p = min(1.0, 2 * sum(math.comb(n, k) for k in range(min(b, c) + 1)) / 2 ** n) if n else 1.0
    return b, c, p


# ------------------------------------------------------------------ tuning
def tune(score_cal):
    """score_cal: {pid: loss magnitude}. Grid over thresholds, area-weighted F1 on calibration points."""
    vals = sorted({v for v in score_cal.values() if v and v > 0})
    grid = sorted(set([0.0] + vals))
    best = (-1, None)
    for t in grid:
        rows = [(int(cal[p]['stratum']), int(v is not None and v > t), int(cal[p]['ref_jhum']))
                for p, v in score_cal.items() if p in cal]
        f1 = estimate(rows, ACAL)['F1']
        if f1 == f1 and f1 > best[0]:
            best = (f1, t)
    return best


def main(variants):
    out = {}
    rule25 = [(int(v25[p]['stratum']), int(v25[p]['map_landsat']), int(v25[p]['ref_jhum'])) for p in sorted(v25)]
    rule10 = [(int(v10[p]['stratum']), int(v10[p]['map_jhum']), int(v10[p]['ref_jhum'])) for p in sorted(v10)]
    out['rule'] = {'2025': estimate(rule25, A25), '2010': estimate(rule10, A10)}
    for name, (fn, src) in variants.items():
        data = json.load(open(EE + src))
        scal = {p: fn(data.get(p), 2025) for p in cal}
        f1c, t = tune(scal)
        res = {'threshold': t, 'cal_F1': f1c, 'missing': sum(v is None for v in scal.values())}
        for Y, ref, areas, rule_rows in (('2025', v25, A25, rule25), ('2010', v10, A10, rule10)):
            rows = []
            for p in sorted(ref):
                sc = fn(data.get(p), int(Y))
                rows.append((int(ref[p]['stratum']), int(sc is not None and sc > t), int(ref[p]['ref_jhum'])))
            e = estimate(rows, areas)
            ci, dci = boot_paired(rule_rows, rows, areas)
            e['ci'] = ci; e['diff_rule_minus_baseline'] = {k: (estimate(rule_rows, areas)[k] - e[k], dci[k]) for k in dci}
            e['mcnemar_rule_correct_baseline_wrong, baseline_correct_rule_wrong, p'] = mcnemar(rule_rows, rows)
            e['n_mapped_points'] = sum(m for _, m, _ in rows)
            res[Y] = e
        out[name] = res
        print(name, 'thr', round(t, 4), 'calF1', round(f1c, 3), '| 2025 UA/PA/F1',
              *(round(res['2025'][k], 3) for k in ('UA', 'PA', 'F1')), '| 2010', *(round(res['2010'][k], 3) for k in ('UA', 'PA', 'F1')),
              flush=True)
    json.dump(out, open(EE + 'baseline_scores.json', 'w'), indent=1, default=float)
    return out


if __name__ == '__main__':
    import os
    V = {}
    for ps in ('standard', 'fast', 'fast8', 'fast10', 'fast12', 'fast16'):
        for end in (2025, 2026):
            for mode in ('segment', 'annual'):
                f = f'lt_{ps}_{end}.json'
                if os.path.exists(EE + f):
                    V[f'LandTrendr {ps} {mode}, record to {end}'] = ((lambda m: (lambda a, Y: lt_loss(a, Y, m)))(mode), f)
    for mo in (6, 3):
        for end in ('2025', '2026'):
            f = f'ccdc_{mo}_{end}.json'
            if os.path.exists(EE + f):
                V[f'CCDC minObs {mo}, record to {end}'] = (ccdc_loss, f)
    main(V)
