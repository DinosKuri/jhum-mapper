"""Replace geoBoundaries districts/state with the official Mizoram district layer (LGD-coded)."""
import sys, geopandas as gpd, pandas as pd
from shapely.ops import unary_union
SRC, *OUTS = sys.argv[1:]
d = gpd.read_file(SRC).to_crs('EPSG:32646')
d['geometry'] = d.geometry.buffer(0)
fix = {'SAIHA': 'Siaha'}
d['name'] = d.DISTRICT.map(lambda s: fix.get(s, s.title()))
d['lgd_code'] = d.DIST_LGD
state = unary_union(d.geometry.buffer(5)).buffer(-5)          # close hairline gaps between districts
parts = list(state.geoms) if state.geom_type == 'MultiPolygon' else [state]
holes = sum(len(p.interiors) for p in parts)
state = max(parts, key=lambda p: p.area)                     # drop zero-area slivers
print('state parts:', len(parts), ' interior holes:', holes, ' part areas km2:', [round(p.area/1e6, 2) for p in parts][:5])
st = gpd.GeoDataFrame({'name': ['Mizoram']}, geometry=[state], crs='EPSG:32646')
print('state area km2 (equal-area EPSG:6933):', round(st.to_crs(6933).area.iloc[0] / 1e6, 1))
print(d[['name', 'lgd_code']].assign(km2=d.to_crs(6933).area.div(1e6).round(1)).to_string(index=False))
for o in OUTS:
    d[['name', 'lgd_code', 'geometry']].to_file(f'{o}/districts.gpkg')
    st.to_file(f'{o}/mizoram.gpkg')
