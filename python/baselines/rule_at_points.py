"""Anniversary rule (v6) at the 520 validation points of 2025, with both branches kept separately.

Writes rule_points_2025.json: per point, wide-window NBR_post (pw) and dNBR (dw), narrow-window
NBR_post (pn) and dNBR (dn), and the rule class (jh). The class reproduces map_landsat of the data
deposit at every point. Usage: python rule_at_points.py
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ee
from ee_common import init, landsat, nbr, water_free, load_points, points_fc


def pct(col):
    return col.map(nbr).reduce(ee.Reducer.percentile([25])).rename('NBR')


def main():
    init()
    pts = [p for p in load_points() if p['set'] == 'v25']
    fc = points_fc(pts)
    LS = landsat(fc.geometry().bounds())
    wpre, wpost = LS.filterDate('2023-11-01', '2024-06-01'), LS.filterDate('2024-11-01', '2025-06-01')
    npre, npost = LS.filterDate('2024-03-01', '2024-06-01'), LS.filterDate('2025-03-01', '2025-06-01')
    pw = pct(wpost); dw = pct(wpre).subtract(pw)
    pn = pct(npost); dn = pct(npre).subtract(pn)
    A = pw.lte(0.26).And(dw.gte(0.20))
    B = pn.lte(0.14).And(dn.gte(0.46))
    img = ee.Image.cat([pw.rename('pw'), dw.rename('dw'), pn.rename('pn'), dn.rename('dn'),
                        A.Or(B).rename('jh')]).updateMask(water_free())
    feats = img.reduceRegions(collection=fc, reducer=ee.Reducer.first(), scale=30, crs='EPSG:32646').getInfo()['features']
    out = {f['properties']['pid']: {k: f['properties'].get(k) for k in ('pw', 'dw', 'pn', 'dn', 'jh')} for f in feats}
    json.dump(out, open('rule_points_2025.json', 'w'), indent=0)
    print(len(out), 'points')


if __name__ == '__main__':
    main()
