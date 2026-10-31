"""Shared Earth Engine setup for the baseline comparison.

Preprocessing is a line-for-line port of gee/00_final_map_mizoram.js (jhum-mapper v1.2.0):
Collection 2 Level-2 surface reflectance, QA_PIXEL cloud / dilated cloud / shadow mask,
TM and ETM+ harmonized to OLI with Roy et al. (2016) OLS coefficients.
"""
import json, csv, os
import ee

# Earth Engine credentials: set EE_PROJECT, and optionally EE_SERVICE_ACCOUNT_KEY (path to a key file).
# Reference data: set JHUM_DATA to the folder of the data deposit (Khiangte, 2026b).
KEY = os.environ.get('EE_SERVICE_ACCOUNT_KEY')
DD = os.environ.get('JHUM_DATA', 'jhum-mapper-data') + '/'


def init():
    if KEY:
        info = json.load(open(KEY))
        ee.Initialize(ee.ServiceAccountCredentials(info['client_email'], KEY), project=info['project_id'])
    else:
        ee.Initialize(project=os.environ['EE_PROJECT'])


NAMES = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
SLOPES = [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]
INTER = [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]
OLD = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']
NEW = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']


def prep(img, bands, harmonize):
    img = ee.Image(img)
    qa = img.select('QA_PIXEL')
    clear = (qa.bitwiseAnd(1 << 1).eq(0)
             .And(qa.bitwiseAnd(1 << 3).eq(0))
             .And(qa.bitwiseAnd(1 << 4).eq(0)))
    out = img.select(bands, NAMES).multiply(0.0000275).add(-0.2).updateMask(clear)
    if harmonize:
        out = out.multiply(SLOPES).add(INTER)
    return ee.Image(out.toFloat().copyProperties(img, ['system:time_start']))


def landsat(region):
    def c(cid, bands, h):
        return ee.ImageCollection(cid).filterBounds(region).map(lambda i: prep(i, bands, h))
    return (c('LANDSAT/LT05/C02/T1_L2', OLD, True)
            .merge(c('LANDSAT/LE07/C02/T1_L2', OLD, True))
            .merge(c('LANDSAT/LC08/C02/T1_L2', NEW, False))
            .merge(c('LANDSAT/LC09/C02/T1_L2', NEW, False)))


def nbr(i):
    i = ee.Image(i)
    return i.normalizedDifference(['NIR', 'SWIR2']).rename('NBR').copyProperties(i, ['system:time_start'])


WATER_FREE = None


def water_free():
    return ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0).lt(50)


def load_points():
    """All 1,080 reference points with a 'set' property: cal (400), v25 (520), v10 (160)."""
    pts = []
    for r in csv.DictReader(open(DD + 'calibration_points_400_2025.csv')):
        pts.append(dict(pid=r['point_id'], set='cal', lon=float(r['lon']), lat=float(r['lat'])))
    for r in csv.DictReader(open(DD + 'validation_points_520_2025.csv')):
        pts.append(dict(pid=r['point_id'], set='v25', lon=float(r['lon']), lat=float(r['lat'])))
    for r in csv.DictReader(open(DD + 'validation_points_160_2010.csv')):
        pts.append(dict(pid=r['point_id'], set='v10', lon=float(r['lon']), lat=float(r['lat'])))
    return pts


def points_fc(pts):
    return ee.FeatureCollection([ee.Feature(ee.Geometry.Point([p['lon'], p['lat']]), {'pid': p['pid'], 'set': p['set']})
                                 for p in pts])
