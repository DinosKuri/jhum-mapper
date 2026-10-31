"""LandTrendr baseline, run at the 1,080 reference points.

Input: one composite per clearing year Y = 25th percentile of per-scene NBR over 1 Nov (Y-1) - 31 May (Y),
the same harmonized collection and window as the post-season composite of the anniversary rule.
NBR is multiplied by -1000 so that vegetation loss is a positive change (LT-GEE convention).

Parameter sets:
  standard  Kennedy et al. (2018) Table 1 / LT-GEE examples
  fast      for fast-recovering disturbance (Kennedy et al., 2010): recovery threshold and
            one-year-recovery prevention switched off, despiking off, pval 0.1
  fast8 / fast10 / fast12 / fast16   the fast set with 8, 10, 12 and 16 segments, because a jhum pixel
            can be cleared and regrow several times in 37 years and each cycle needs two segments
Record ends: 2025 (same information as the rule) and 2026 (one confirming season after 2025).

Writes lt_<set>_<end>.json: per point, the LandTrendr array (year, source, fitted, vertex).
Usage: python run_landtrendr.py
"""
import os, sys, json, time, concurrent.futures as cf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ee
from ee_common import init, landsat, nbr, load_points, points_fc

PARAMS = {
    'standard': dict(maxSegments=6, spikeThreshold=0.9, vertexCountOvershoot=3, preventOneYearRecovery=True,
                     recoveryThreshold=0.25, pvalThreshold=0.05, bestModelProportion=0.75, minObservationsNeeded=6),
    'fast10': dict(maxSegments=10, spikeThreshold=1.0, vertexCountOvershoot=3, preventOneYearRecovery=False,
                   recoveryThreshold=1.0, pvalThreshold=0.1, bestModelProportion=0.75, minObservationsNeeded=6),
    'fast8': dict(maxSegments=8, spikeThreshold=1.0, vertexCountOvershoot=3, preventOneYearRecovery=False,
                  recoveryThreshold=1.0, pvalThreshold=0.1, bestModelProportion=0.75, minObservationsNeeded=6),
    'fast12': dict(maxSegments=12, spikeThreshold=1.0, vertexCountOvershoot=3, preventOneYearRecovery=False,
                   recoveryThreshold=1.0, pvalThreshold=0.1, bestModelProportion=0.75, minObservationsNeeded=6),
    'fast16': dict(maxSegments=16, spikeThreshold=1.0, vertexCountOvershoot=3, preventOneYearRecovery=False,
                   recoveryThreshold=1.0, pvalThreshold=0.1, bestModelProportion=0.75, minObservationsNeeded=6),
    'fast': dict(maxSegments=6, spikeThreshold=1.0, vertexCountOvershoot=3, preventOneYearRecovery=False,
                 recoveryThreshold=1.0, pvalThreshold=0.1, bestModelProportion=0.75, minObservationsNeeded=6),
}


def annual(LS, y):
    col = LS.filterDate(ee.Date.fromYMD(y - 1, 11, 1), ee.Date.fromYMD(y, 6, 1)).map(nbr)
    img = col.reduce(ee.Reducer.percentile([25])).rename('NBR').multiply(-1000).toFloat()
    return img.set('system:time_start', ee.Date.fromYMD(y, 3, 1).millis())


def run(pset, end, chunk):
    fc = points_fc(chunk)
    LS = landsat(fc.geometry().bounds())
    col = ee.ImageCollection([annual(LS, y) for y in range(1989, end + 1)])
    lt = ee.Algorithms.TemporalSegmentation.LandTrendr(timeSeries=col, **PARAMS[pset]).select('LandTrendr')
    res = lt.reduceRegions(collection=fc, reducer=ee.Reducer.first(), scale=30, crs='EPSG:32646', tileScale=16).getInfo()['features']
    return {f['properties']['pid']: f['properties'].get('first') for f in res}


def safe(pset, end, chunk):
    for attempt in range(3):
        try:
            return run(pset, end, chunk)
        except Exception:
            time.sleep(5)
    out = {}
    for p in chunk:
        try:
            out.update(run(pset, end, [p]))
        except Exception as e:
            print('failed', p['pid'], str(e)[:80], flush=True)
    return out


def main(psets=('standard', 'fast')):
    init()
    pts = load_points()
    chunks = [pts[i:i + 8] for i in range(0, len(pts), 8)]
    for pset in psets:
        for end in (2025, 2026):
            out, t = {}, time.time()
            with cf.ThreadPoolExecutor(6) as ex:
                for r in ex.map(lambda c: safe(pset, end, c), chunks):
                    out.update(r)
            json.dump(out, open(f'lt_{pset}_{end}.json', 'w'))
            print(pset, end, len(out), 'points', round(time.time() - t), 's', flush=True)


if __name__ == '__main__':
    main(tuple(sys.argv[1:]) or ('standard', 'fast'))
