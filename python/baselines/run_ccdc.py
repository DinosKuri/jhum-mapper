"""CCDC baseline, run at the 1,080 reference points.

Input: every clear scene of the harmonized Landsat 5/7/8/9 collection (same preprocessing as the rule),
bands GREEN, RED, NIR, SWIR1, SWIR2 and NBR, reflectance scaled 0-1.
Parameters follow the Earth Engine defaults and the gee-ccdc-tools tutorial (Arevalo et al., 2020):
breakpointBands GREEN/RED/NIR/SWIR1/SWIR2, tmaskBands GREEN/SWIR1, chiSquareProbability 0.99,
minNumOfYearsScaler 1.33, lambda 0.002 (20 scaled for 0-1 reflectance), dateFormat 1 (fractional year).
Variants: minObservations 6 (default) and 3 (cf. Chen et al., 2023); record end 31 May 2025
(same information as the rule) and 30 September 2026 (all data now available).

Writes ccdc_<minobs>_<end>.json: per point, tBreak, changeProb and NBR_magnitude arrays.
"""
import os, sys, json, time, concurrent.futures as cf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ee
from ee_common import init, landsat, load_points, points_fc

ENDS = {'2025': '2025-06-01', '2026': '2026-10-01'}


def add_nbr(i):
    return i.addBands(i.normalizedDifference(['NIR', 'SWIR2']).rename('NBR'))


def run(minobs, end, chunk):
    fc = points_fc(chunk)
    col = landsat(fc.geometry().bounds()).filterDate('1988-01-01', ENDS[end]).map(add_nbr) \
        .select(['GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'NBR'])
    c = ee.Algorithms.TemporalSegmentation.Ccdc(
        collection=col, breakpointBands=['GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'], tmaskBands=['GREEN', 'SWIR1'],
        minObservations=minobs, chiSquareProbability=0.99, minNumOfYearsScaler=1.33, dateFormat=1,
        **{'lambda': 0.002}, maxIterations=25000)
    img = c.select(['tStart', 'tEnd', 'tBreak', 'changeProb', 'NBR_magnitude'])
    res = img.reduceRegions(collection=fc, reducer=ee.Reducer.first(), scale=30, crs='EPSG:32646', tileScale=16).getInfo()['features']
    return {f['properties']['pid']: {k: f['properties'].get(k) for k in ('tStart', 'tEnd', 'tBreak', 'changeProb', 'NBR_magnitude')}
            for f in res}


def safe(minobs, end, chunk):
    for attempt in range(3):
        try:
            return run(minobs, end, chunk)
        except Exception as e:
            err = e
            time.sleep(5)
    out = {}
    for p in chunk:                                # fall back to one point at a time
        try:
            out.update(run(minobs, end, [p]))
        except Exception as e:
            print('failed', p['pid'], str(e)[:80], flush=True)
    return out


def main(variants):
    init()
    pts = load_points()
    chunks = [pts[i:i + 5] for i in range(0, len(pts), 5)]
    for minobs, end in variants:
        out, t = {}, time.time()
        with cf.ThreadPoolExecutor(8) as ex:
            futs = [ex.submit(safe, minobs, end, c) for c in chunks]
            for fu in cf.as_completed(futs):
                out.update(fu.result())
        json.dump(out, open(f'ccdc_{minobs}_{end}.json', 'w'))
        print(minobs, end, len(out), 'points', round(time.time() - t), 's', flush=True)


if __name__ == '__main__':
    main([(6, '2026'), (6, '2025'), (3, '2026'), (3, '2025')])
