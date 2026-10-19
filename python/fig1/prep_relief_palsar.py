"""Relief for Figure 1: ALOS PALSAR 12.5 m DEM inside its footprint, Copernicus GLO-30 beyond.
Output grid: 25 m, EPSG:32646, state bounds + 7 km. Also reports elevation range inside the state
from the native 12.5 m PALSAR grid.
Usage: python prep_relief_palsar.py <palsar_dem.tif> <copernicus_dem_30m.tif> <mizoram.gpkg> <out_dir>
"""
import sys, numpy as np, geopandas as gpd, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from rasterio.features import geometry_mask

PAL, COP, MIZ, OUT = sys.argv[1:5]
# The PALSAR RTC DEM is in WGS84 ellipsoid heights. Over Mizoram, Copernicus GLO-30 (EGM2008, i.e.
# above mean sea level) minus PALSAR is +48.5 m (median; 47.9-49.2 m across the state), the local
# geoid separation. Adding it converts PALSAR to orthometric height.
GEOID_OFFSET = float(sys.argv[5]) if len(sys.argv) > 5 else 48.5
CRS, RES, PAD = 'EPSG:32646', 25.0, 7000
state = gpd.read_file(MIZ).to_crs(CRS).geometry.iloc[0]

# native-resolution statistics inside the state
with rasterio.open(PAL) as r:
    z = r.read(1); z[(z == r.nodata) | (z < -60)] = np.nan; z += GEOID_OFFSET
    ins = ~geometry_mask([state], z.shape, r.transform)
    v = z[ins & np.isfinite(z)]
    print('PALSAR in state: n=%d  min %.0f  p0.1 %.0f  median %.0f  p99.9 %.0f  max %.0f  nodata share %.4f'
          % (v.size, v.min(), np.percentile(v, 0.1), np.median(v), np.percentile(v, 99.9), v.max(),
             1 - v.size / ins.sum()))
    ptr, pcrs = r.transform, r.crs
xmin, ymin, xmax, ymax = state.bounds
xmin, ymin, xmax, ymax = xmin - PAD, ymin - PAD, xmax + PAD, ymax + PAD
W, H = int(np.ceil((xmax - xmin) / RES)), int(np.ceil((ymax - ymin) / RES))
tr = from_origin(xmin, ymax, RES, RES)

pal = np.full((H, W), np.nan, 'float32')
reproject(z, pal, src_transform=ptr, src_crs=pcrs, dst_transform=tr, dst_crs=CRS,
          resampling=Resampling.average, src_nodata=np.nan, dst_nodata=np.nan)
del z
cop = np.full((H, W), np.nan, 'float32')
with rasterio.open(COP) as r:
    reproject(rasterio.band(r, 1), cop, dst_transform=tr, dst_crs=CRS, resampling=Resampling.bilinear,
              dst_nodata=np.nan)
dem = np.where(np.isfinite(pal), pal, cop)
src = np.where(np.isfinite(pal), 1, np.where(np.isfinite(cop), 2, 0)).astype('uint8')
print('grid %d x %d; PALSAR share %.3f, Copernicus share %.3f' % (W, H, (src == 1).mean(), (src == 2).mean()))

def hs(zz, az, alt=45, zf=1.4):
    zz = np.where(np.isnan(zz), np.nanmean(zz), zz) * zf
    gy, gx = np.gradient(zz, RES)
    slope = np.arctan(np.hypot(gx, gy)); aspect = np.arctan2(-gx, gy)
    a, b = np.radians(360 - az + 90), np.radians(alt)
    return np.clip(np.sin(b) * np.cos(slope) + np.cos(b) * np.sin(slope) * np.cos(a - aspect), 0, 1)
shade = 0.50 * hs(dem, 315) + 0.20 * hs(dem, 270) + 0.20 * hs(dem, 360) + 0.10 * hs(dem, 225)

prof = dict(driver='GTiff', height=H, width=W, count=1, dtype='float32', crs=CRS, transform=tr,
            compress='deflate', nodata=np.nan, tiled=True)
for name, arr in (('dem_25m', dem), ('hillshade_25m', shade)):
    with rasterio.open(f'{OUT}/{name}.tif', 'w', **prof) as o:
        o.write(arr.astype('float32'), 1)
