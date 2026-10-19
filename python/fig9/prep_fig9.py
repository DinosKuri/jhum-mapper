"""Clean the v1 Earth Engine exports and build 1 km summaries for Figure 9.

Fallow-age fix: the v1 export wrote 'never cleared' and 'outside the state'
as 0, the same value as 'cleared in 2025'. A pixel was ever cleared if and
only if jhum_count > 0, so never-cleared pixels are set to nodata there.
"""
import sys, numpy as np, rasterio
from rasterio.warp import reproject, Resampling
from affine import Affine

IN, OUT = sys.argv[1], sys.argv[2]
rd = lambda n: rasterio.open(f'{IN}/{n}')
with rd('FIG9_jhum_2025.tif') as r:
    jhum = r.read(1); prof = r.profile; tr = r.transform; crs = r.crs
count = rd('FIG9_jhum_count_1989_2025.tif').read(1)
fallow = rd('FIG9_fallow_age_1989_2025.tif').read(1).astype('float32')
cycle = rd('FIG9_cycle_mean_1998_2025.tif').read(1)

ever = count > 0
fallow[~ever] = np.nan
assert np.all(fallow[jhum == 1] == 0) and np.all(jhum[fallow == 0] == 1)
print('2025 mapped, km2: %.1f   ever cleared 1989-2025, km2: %.0f' % ((jhum == 1).sum() * 9e-4, ever.sum() * 9e-4))

# inside-state mask from the export footprint (jhum raster carries 0 inside, 0 outside too),
# so use the geoBoundaries outline for the denominator
import geopandas as gpd
from rasterio.features import geometry_mask
miz = gpd.read_file(sys.argv[3]).to_crs(crs)
inside = ~geometry_mask(miz.geometry, jhum.shape, tr)

def to_km(a, how=Resampling.average, res=1000.0):
    W = int(np.ceil(a.shape[1] * 30 / res)); H = int(np.ceil(a.shape[0] * 30 / res))
    dst = np.full((H, W), np.nan, 'float32')
    reproject(a.astype('float32'), dst, src_transform=tr, src_crs=crs,
              dst_transform=Affine(res, 0, tr.c, 0, -res, tr.f), dst_crs=crs,
              resampling=how, src_nodata=np.nan, dst_nodata=np.nan)
    return dst, Affine(res, 0, tr.c, 0, -res, tr.f)

j = np.where(inside, jhum.astype('float32'), np.nan)
e = np.where(inside, ever.astype('float32'), np.nan)
jp, ktr = to_km(j); ep, _ = to_km(e)
fm, _ = to_km(fallow)                                  # mean over cleared pixels only
fmed, _ = to_km(fallow, Resampling.med)
jp *= 100; ep *= 100
print('1 km cells with any 2025 jhum: %d; percentiles of jhum%% where >0: %s' %
      ((jp > 0).sum(), np.nanpercentile(jp[jp > 0], [25, 50, 75, 90, 99]).round(2)))
print('ever-cleared %% percentiles: %s' % np.nanpercentile(ep, [10, 25, 50, 75, 90]).round(1))
print('mean fallow age percentiles: %s' % np.nanpercentile(fm, [10, 25, 50, 75, 90]).round(1))

p = dict(driver='GTiff', dtype='float32', crs=crs, compress='deflate', nodata=np.nan)
with rasterio.open(f'{OUT}/state_1km.tif', 'w', height=jp.shape[0], width=jp.shape[1], count=4,
                   transform=ktr, **p) as o:
    for i, (a, n) in enumerate([(jp, 'jhum_pct_2025'), (ep, 'ever_cleared_pct'),
                                (fm, 'fallow_age_mean'), (fmed, 'fallow_age_median')], 1):
        o.write(a, i); o.set_band_description(i, n)
with rasterio.open(f'{OUT}/fallow_age_30m_clean.tif', 'w', height=fallow.shape[0], width=fallow.shape[1],
                   count=1, transform=tr, **p) as o:
    o.write(fallow, 1); o.set_band_description(1, 'years since last clearing; nodata = never cleared or outside')
