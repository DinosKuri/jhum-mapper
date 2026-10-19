"""Prepare terrain and boundary layers for Figure 1 (EPSG:32646)."""
import sys, glob, numpy as np, geopandas as gpd, rasterio
from rasterio.merge import merge
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.transform import from_origin

SRC, OUT = sys.argv[1], sys.argv[2]
CRS = 'EPSG:32646'; RES = 60.0

# --- boundaries ---------------------------------------------------------
a1 = gpd.read_file(f'{SRC}/geoBoundaries-IND-ADM1.geojson')
a2 = gpd.read_file(f'{SRC}/geoBoundaries-IND-ADM2.geojson')
miz = a1[a1.shapeName == 'Mizoram'].to_crs(CRS)
d = a2.to_crs(CRS)
d = d[d.representative_point().within(miz.geometry.iloc[0])].copy()
d['name'] = d.shapeName.replace({'Saiha': 'Siaha'})
print('districts:', len(d), sorted(d.name))
print('state area km2:', round(miz.area.sum() / 1e6, 1))
miz.to_file(f'{OUT}/mizoram.gpkg'); d[['name', 'geometry']].to_file(f'{OUT}/districts.gpkg')

ne = a1[a1.shapeName.isin(['Assam', 'Manipur', 'Tripura', 'Meghālaya', 'Nāgāland',
                           'Arunāchal Pradesh', 'West Bengal', 'Mizoram'])]
ne.to_file(f'{OUT}/ind_ne_states.gpkg')
bgd = gpd.read_file(f'{SRC}/geoBoundaries-BGD-ADM0.geojson')
mmr = gpd.read_file(f'{SRC}/geoBoundaries-MMR-ADM1.geojson').dissolve()
bgd.to_file(f'{OUT}/bgd.gpkg'); mmr.to_file(f'{OUT}/mmr.gpkg')
# India outline for the locator inset of Figure 1
a1.dissolve().to_file(f'{OUT}/india.gpkg')

# --- DEM mosaic to UTM 46N at 60 m ---------------------------------------
xmin, ymin, xmax, ymax = miz.total_bounds
pad = 14000
xmin, ymin, xmax, ymax = xmin - pad, ymin - pad, xmax + pad, ymax + pad
W = int(np.ceil((xmax - xmin) / RES)); H = int(np.ceil((ymax - ymin) / RES))
dst_tr = from_origin(xmin, ymax, RES, RES)

srcs = [rasterio.open(f) for f in sorted(glob.glob(f'{SRC}/Copernicus_DSM_*.tif'))]
mos, tr = merge(srcs, bounds=(91.9, 21.7, 93.8, 24.8))
dem = np.full((H, W), np.nan, 'float32')
reproject(mos[0].astype('float32'), dem, src_transform=tr, src_crs=srcs[0].crs,
          dst_transform=dst_tr, dst_crs=CRS, resampling=Resampling.bilinear,
          src_nodata=srcs[0].nodata, dst_nodata=np.nan)
dem[dem < -50] = np.nan

# --- multi-directional hillshade (softer than single-azimuth) ------------
def hs(z, az, alt=45, zf=1.5):
    zz = np.where(np.isnan(z), np.nanmean(z), z) * zf
    gy, gx = np.gradient(zz, RES)
    slope = np.arctan(np.hypot(gx, gy)); aspect = np.arctan2(-gx, gy)
    az, alt = np.radians(360 - az + 90), np.radians(alt)
    return np.clip(np.sin(alt) * np.cos(slope) +
                   np.cos(alt) * np.sin(slope) * np.cos(az - aspect), 0, 1)
shade = (0.50 * hs(dem, 315) + 0.20 * hs(dem, 270) + 0.20 * hs(dem, 360) + 0.10 * hs(dem, 225))

prof = dict(driver='GTiff', height=H, width=W, count=1, dtype='float32', crs=CRS,
            transform=dst_tr, compress='deflate', nodata=np.nan)
with rasterio.open(f'{OUT}/dem_60m.tif', 'w', **prof) as o: o.write(dem, 1)
with rasterio.open(f'{OUT}/hillshade_60m.tif', 'w', **prof) as o: o.write(shade.astype('float32'), 1)

# elevation range inside the state
from rasterio.features import geometry_mask
inside = ~geometry_mask(miz.geometry, (H, W), dst_tr)
print('elevation in state: min %.0f  max %.0f  median %.0f' %
      (np.nanmin(dem[inside]), np.nanmax(dem[inside]), np.nanmedian(dem[inside])))
