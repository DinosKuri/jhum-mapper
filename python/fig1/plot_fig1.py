"""Figure 1 (v3): study area and validation design — WGS 84 (decimal degrees), 170 mm wide.

a, b  Mizoram, hillshaded relief tinted by elevation (ALOS PALSAR inside the state, Copernicus GLO-30 outside),
      validation points by stratum for 2025 and 2010
c     India locator (outline only), Mizoram in dark red; legend stacked below
Text in Times New Roman (Liberation Serif where Times is not installed).

Usage: python plot_fig1.py <layers_dir> <points_dir> <out_prefix>

<points_dir>/fig1_validation_points.csv      2025 sample: lon, lat, stratum (A-E or 1-5)
<points_dir>/fig11_validation_2010_points.csv 2010 sample: lon, lat, stratum (1-3)
Stratum counts and n in the titles and legend are read from these files.
"""
import sys
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.warp import reproject, Resampling
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from rasterio.features import geometry_mask

LAY, PTS, OUT = sys.argv[1:4]
MM = 1 / 25.4
PT_PER_MM = 72 / 25.4
installed = {f.name for f in font_manager.fontManager.ttflist}
SERIF = 'Times New Roman' if 'Times New Roman' in installed else 'Liberation Serif'
mpl.rcParams.update({'font.family': SERIF, 'font.size': 8, 'axes.linewidth': 0.6, 'pdf.fonttype': 42, 'lines.dashed_pattern': (3.5, 2.0)})

FS_HEAD = 9.5            # panel letter and panel title (same size, both bold)
FS_LEG_H = 8.5           # legend block headings
FS_TXT = 8.0             # legend text, axis labels
FS_DIST = 6.5            # district names
FS_NEIGH = 6.5           # neighbouring units (caps)
INK, GREY, DGREY = '#000000', '#6E6E6E', '#4D4D4D'
HALO = [pe.withStroke(linewidth=1.4, foreground='white')]
# 2025: dual-map design (Landsat and Sentinel-2 maps), five strata
STRATA25 = {'A': ('#B03A2E', 'Jhum on both maps'), 'B': ('#1E9E8A', 'Landsat only'),
            'C': ('#7D4BA0', 'Sentinel-2 only'), 'D': ('#E8A21F', 'Neither, within 300 m'),
            'E': ('#2C6FA6', 'Neither, beyond 300 m')}
# 2010: single-map design, three strata
STRATA10 = {1: ('#B03A2E', 'Mapped jhum'), 2: ('#E8A21F', 'Within 300 m of jhum'), 3: ('#2C6FA6', 'Beyond 300 m')}
LW_STATE = 0.5 * PT_PER_MM          # 0.5 mm
LW_DIST = 0.3 * PT_PER_MM           # 0.3 mm
LW_THIN = 0.45
DOT = (1.4 * PT_PER_MM) ** 2        # 1.4 mm diameter, as marker area (pt²)
Z_LO, Z_HI = 40, 2000

# ------------------------------------------------------------------ vectors (all to WGS 84)
G = 'EPSG:4326'
miz = gpd.read_file(f'{LAY}/mizoram.gpkg').to_crs(G); STATE = miz.geometry.iloc[0]
dist = gpd.read_file(f'{LAY}/districts.gpkg').to_crs(G)
bgd = gpd.read_file(f'{LAY}/bgd.gpkg').to_crs(G); mmr = gpd.read_file(f'{LAY}/mmr.gpkg').to_crs(G)
ne = gpd.read_file(f'{LAY}/ind_ne_states.gpkg').to_crs(G)
ne_other = ne[ne.shapeISO != 'IN-MZ']
india = gpd.read_file(f'{LAY}/india.gpkg').to_crs(G)
def pts(path):
    d = pd.read_csv(path)
    d.columns = [c.strip().lower() for c in d.columns]
    lon = next(c for c in ('lon', 'longitude', 'x') if c in d.columns)
    lat = next(c for c in ('lat', 'latitude', 'y') if c in d.columns)
    return gpd.GeoDataFrame(d, geometry=gpd.points_from_xy(d[lon], d[lat]), crs=G)
p25 = pts(f'{PTS}/fig1_validation_points.csv'); p10 = pts(f'{PTS}/fig11_validation_2010_points.csv')
p25['stratum'] = p25.stratum.astype(str).str.strip().str.upper().replace({'1': 'A', '2': 'B', '3': 'C', '4': 'D', '5': 'E'})
p10['stratum'] = p10.stratum.astype(int)
bad = set(p25.stratum) - set(STRATA25)
assert not bad, f'2025 file has strata {bad}; expected A-E (or 1-5). Is this the 520-point table?'
assert set(p10.stratum) <= set(STRATA10), '2010 file must have strata 1-3'
N25 = p25.stratum.value_counts().to_dict(); N10 = p10.stratum.value_counts().to_dict()
print('2025 points', len(p25), N25, '| 2010 points', len(p10), N10)

# map window (degrees); aspect corrected for latitude so the map is not stretched
LON0, LON1, LAT0, LAT1 = 92.18, 93.50, 21.88, 24.58
KLAT = 1 / np.cos(np.radians((LAT0 + LAT1) / 2))       # data aspect for plate carrée display

# ------------------------------------------------------------------ relief, resampled to a WGS 84 grid
RES = 0.0005                                           # ≈ 55 m; print resolution is ≈ 80 m per dot
W = int(round((LON1 - LON0) / RES)); H = int(round((LAT1 - LAT0) / RES))
dst_tr = rasterio.Affine(RES, 0, LON0, 0, -RES, LAT1)
def warp(path, fill):
    out = np.full((H, W), np.nan, 'float32')
    with rasterio.open(path) as r:
        reproject(rasterio.band(r, 1), out, dst_transform=dst_tr, dst_crs=G,
                  resampling=Resampling.average, dst_nodata=np.nan)
    return np.where(np.isfinite(out), out, fill).astype('float32')
dem = warp(f'{LAY}/dem_25m.tif', Z_LO); shade = warp(f'{LAY}/hillshade_25m.tif', 0.8)
inside = ~geometry_mask([STATE], (H, W), dst_tr)
TINT = LinearSegmentedColormap.from_list('tint', [      # muted sage green (40 m) -> soft buff-brown (2,000 m)
    (0.00, '#b9cba9'), (0.25, '#cdd6ae'), (0.50, '#e0dbb6'), (0.75, '#dcc6a0'), (1.00, '#c9a985')])
LUT = TINT(np.linspace(0, 1, 256))[:, :3].astype('float32')
zi = (np.clip((dem - Z_LO) / (Z_HI - Z_LO), 0, 1) * 255).astype('uint8')
s = shade[..., None]
rgb_in = LUT[zi] * (0.55 + 0.55 * s)                    # relief shading inside the state
rgb_out = (0.86 + 0.10 * s) * np.array([0.93, 0.93, 0.93], 'float32')   # flat light grey, faint relief
RGB = (np.clip(np.where(inside[..., None], rgb_in, rgb_out), 0, 1) * 255).astype('uint8')
EXT = (LON0, LON1, LAT0, LAT1)
del dem, shade, zi, s, rgb_in, rgb_out

# ------------------------------------------------------------------ layout (mm)
FIG_W = 190.0
M_L, M_R, M_T, M_B = 2.0, 2.0, 2.0, 2.0
LAT_LAB = 5.5                 # room for latitude labels left of panel a
LON_LAB = 5.0                 # room for longitude labels below the maps
HEAD_H = 6.5                  # panel letter + title band above the maps
GAP_AB, GAP_C = 4.0, 6.0
COL_W = 46.0
MW = (FIG_W - M_L - LAT_LAB - GAP_AB - GAP_C - COL_W - M_R) / 2
MH = MW * (LAT1 - LAT0) * KLAT / (LON1 - LON0)
FIG_H = M_B + LON_LAB + MH + HEAD_H + M_T
fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))
def ax_mm(x, y, w, h): return fig.add_axes([x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H])
def ftxt(x, y, s, **kw): return fig.text(x / FIG_W, y / FIG_H, s, **kw)
XA = M_L + LAT_LAB; XB = XA + MW + GAP_AB; MAP_Y = M_B + LON_LAB; TOP_Y = MAP_Y + MH
COL_X = XB + MW + GAP_C

def heading(x, letter, title):
    t = ftxt(x, TOP_Y + 2.2, letter, fontsize=FS_HEAD, fontweight='bold', va='bottom', ha='left')
    ftxt(x + 3.6, TOP_Y + 2.2, title, fontsize=FS_HEAD, fontweight='bold', va='bottom', ha='left')

def graticule(ax, lat_labels):
    kw = dict(color='#9a9a9a', lw=0.4, ls=(0, (1, 1.6)), zorder=4)
    for lat in np.arange(22.0, 24.51, 0.5):
        ax.plot([LON0, LON1], [lat, lat], **kw)
        if lat_labels:
            ax.text(LON0 - 0.025, lat, f'{lat:.1f}°N', rotation=90, ha='right', va='center', fontsize=FS_TXT)
    for lon in (92.5, 93.0):
        ax.plot([lon, lon], [LAT0, LAT1], **kw)
        ax.text(lon, LAT0 - 0.04, f'{lon:.1f}°E', ha='center', va='top', fontsize=FS_TXT)

# district label positions (lon, lat, rotation), placed inside each district clear of the outline
DPLACE = {'Kolasib': (92.73, 24.15, 0), 'Mamit': (92.43, 23.80, 0), 'Aizawl': (92.79, 23.62, 0),
          'Saitual': (93.09, 23.62, 0), 'Khawzawl': (93.14, 23.42, 90), 'Champhai': (93.31, 23.47, 90),
          'Serchhip': (92.93, 23.17, 0), 'Hnahthial': (92.99, 22.83, 0), 'Lunglei': (92.62, 22.85, 0),
          'Lawngtlai': (92.73, 22.32, 0), 'Siaha': (92.98, 22.30, 0)}
NEIGH = [('ASSAM', 92.45, 24.47, 0), ('MANIPUR', 93.30, 24.38, 0),
         ('BANGLADESH', 92.27, 22.95, 90), ('MYANMAR', 93.40, 22.55, 90)]


# ---- district names: for each district, the spot inside it with the fewest validation points (both years)
from shapely.geometry import box as sbox
ALLPTS = np.vstack([np.c_[p25.geometry.x, p25.geometry.y], np.c_[p10.geometry.x, p10.geometry.y]])
MM_X = MW / (LON1 - LON0); MM_Y = MH / (LAT1 - LAT0)            # mm per degree on the page
CH_MM, TXT_H_MM, DOT_R_MM = FS_DIST * 0.47 / PT_PER_MM, FS_DIST * 0.95 / PT_PER_MM, 0.8
def place_label(name, poly, rot):
    w_mm, h_mm = len(name) * CH_MM, TXT_H_MM
    if rot: w_mm, h_mm = h_mm, w_mm
    hw, hh = w_mm / 2 / MM_X, h_mm / 2 / MM_Y
    rx, ry = DOT_R_MM / MM_X, DOT_R_MM / MM_Y
    inner = poly.buffer(-0.012)
    c = SEED.get(name) or poly.representative_point()
    x0, y0, x1, y1 = poly.bounds; best = None
    for x in np.arange(x0, x1, 0.008):
        for y in np.arange(y0, y1, 0.008):
            if not inner.contains(sbox(x - hw, y - hh, x + hw, y + hh)): continue
            if np.hypot((x - c.x) * MM_X, (y - c.y) * MM_Y) > 14: continue
            n = np.sum((np.abs(ALLPTS[:, 0] - x) < hw + rx) & (np.abs(ALLPTS[:, 1] - y) < hh + ry))
            score = n * 10 + np.hypot((x - c.x) * MM_X, (y - c.y) * MM_Y) * 0.8
            if best is None or score < best[0]: best = (score, x, y, n)
    if best is None: return (c.x, c.y, rot)
    print(f'  {name:10s} points under label: {best[3]}')
    return (best[1], best[2], rot)
from shapely.geometry import Point
SEED = {'Aizawl': Point(92.76, 23.68)}
FORCE_ROT = {'Khawzawl': (90,), 'Champhai': (90,)}
def best_of(name, poly):
    out = []
    for rot in FORCE_ROT.get(name, (0, 90)):
        x, y, _ = place_label(name, poly, rot)
        w_mm, h_mm = len(name) * CH_MM, TXT_H_MM
        if rot: w_mm, h_mm = h_mm, w_mm
        hw, hh = w_mm / 2 / MM_X, h_mm / 2 / MM_Y
        n = np.sum((np.abs(ALLPTS[:, 0] - x) < hw + DOT_R_MM / MM_X) & (np.abs(ALLPTS[:, 1] - y) < hh + DOT_R_MM / MM_Y))
        out.append((n + (2.5 if rot else 0), x, y, rot))      # horizontal preferred on a tie
    n, x, y, rot = min(out)
    print(f'  -> {name:10s} rot {rot:2d}, points under label {int(n)}')
    return (x, y, rot)
LABELS = {r.name: best_of(r.name, r.geometry) for r in dist.itertuples()}

def map_panel(ax, gdf, lat_labels, strata):
    ax.imshow(RGB, extent=EXT, interpolation='antialiased', zorder=1)
    for g in (bgd, mmr): g.boundary.plot(ax=ax, lw=LW_THIN, color='#8c8c8c', zorder=3)          # international
    ne_other.boundary.plot(ax=ax, lw=LW_THIN, color='#8c8c8c', ls='--', zorder=3)     # neighbouring state
    dist.boundary.plot(ax=ax, lw=LW_DIST, color=GREY, zorder=5)                                 # district: one grey line
    miz.boundary.plot(ax=ax, lw=LW_STATE, color=INK, zorder=7, joinstyle='round')               # state
    graticule(ax, lat_labels)
    for name, (x, y, rot) in LABELS.items():
        ax.text(x, y, name, fontsize=FS_DIST, color=DGREY, style='italic', ha='center', va='center',
                rotation=rot, zorder=20, path_effects=HALO)
    for name, x, y, rot in NEIGH:
        ax.text(x, y, name, fontsize=FS_NEIGH, color='#7a7a7a', ha='center', va='center', rotation=rot,
                zorder=9)
    for i, k in enumerate(reversed(list(strata))):      # background strata first, jhum on top
        g = gdf[gdf.stratum == k]
        ax.scatter(g.geometry.x, g.geometry.y, s=DOT, c=strata[k][0], lw=0.35, edgecolors='white',
                   zorder=12 + i)
    ax.set_xlim(LON0, LON1); ax.set_ylim(LAT0, LAT1); ax.set_aspect(KLAT)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(''); ax.set_ylabel('')
    for sp in ax.spines.values(): sp.set_linewidth(0.6); sp.set_color(INK)

axA = ax_mm(XA, MAP_Y, MW, MH); axB = ax_mm(XB, MAP_Y, MW, MH)
map_panel(axA, p25, True, STRATA25); map_panel(axB, p10, False, STRATA10)
heading(XA, 'a', f'Validation sample, 2025 (n = {len(p25)})')
heading(XB, 'b', f'Validation sample, 2010 (n = {len(p10)})')

# one plain scale bar (panel b, lower right) and a small north arrow above it
KM_PER_DEG = 111.32 * np.cos(np.radians(22.1))
def scalebar(ax, lon0, lat0, km=50, h=0.022):
    dl = km / KM_PER_DEG
    for i, c in enumerate((INK, 'white')):
        ax.add_patch(Rectangle((lon0 + i * dl / 2, lat0), dl / 2, h, facecolor=c, edgecolor=INK, lw=0.5, zorder=20))
    for v, xx in ((0, lon0), (km, lon0 + dl)):
        ax.text(xx, lat0 + h + 0.018, f'{v}' + (' km' if v == km else ''), ha='center', va='bottom',
                fontsize=FS_TXT - 0.5, zorder=20, path_effects=HALO)
scalebar(axB, LON1 - 0.10 - 30 / KM_PER_DEG, LAT0 + 0.05, km=30)
nx = LON1 - 0.07
axB.annotate('', xy=(nx, LAT0 + 0.40), xytext=(nx, LAT0 + 0.24),
             arrowprops=dict(arrowstyle='-|>,head_length=0.55,head_width=0.25', lw=0.8, color=INK), zorder=20)
axB.text(nx, LAT0 + 0.42, 'N', ha='center', va='bottom', fontsize=FS_TXT, fontweight='bold', zorder=20,
         path_effects=HALO)

# ------------------------------------------------------------------ c: India locator
heading(COL_X, 'c', 'Location')
IND_CRS = 'EPSG:7755'                                   # India NSF LCC (WGS 84)
ind = india.to_crs(IND_CRS); mz = miz.to_crs(IND_CRS)
ib = ind.total_bounds; pad = 50000
iw, ih = ib[2] - ib[0] + 2 * pad, ib[3] - ib[1] + 2 * pad
LOC_W = 30.0; LOC_H = LOC_W * ih / iw
axC = ax_mm(COL_X + (COL_W - LOC_W) / 2, TOP_Y - LOC_H, LOC_W, LOC_H)
ind.plot(ax=axC, facecolor='white', edgecolor='#7a7a7a', lw=0.5)
mz.plot(ax=axC, facecolor='#8B1A1A', edgecolor='#8B1A1A', lw=0.4, zorder=5)
mc = mz.geometry.iloc[0].representative_point()
axC.annotate('Mizoram', xy=(mc.x, mc.y - 60000), xytext=(mc.x - 230000, mc.y - 640000),
             fontsize=FS_TXT, fontweight='bold', ha='center', va='top',
             arrowprops=dict(arrowstyle='-', lw=0.5, color=INK, shrinkA=1, shrinkB=0), zorder=7)
axC.text(0.40, 0.50, 'INDIA', transform=axC.transAxes, fontsize=FS_LEG_H + 1, color='#7a7a7a',
         ha='center', va='center', fontweight='bold')
axC.set_xlim(ib[0] - pad, ib[2] + pad); axC.set_ylim(ib[1] - pad, ib[3] + pad); axC.set_aspect('equal')
axC.axis('off')

# ------------------------------------------------------------------ legend (three blocks)
y = TOP_Y - LOC_H - 5.0
lg = ax_mm(COL_X, MAP_Y, COL_W, y - MAP_Y); lg.set_xlim(0, COL_W); lg.set_ylim(MAP_Y, y); lg.axis('off')
def head(yy, s): lg.text(0, yy, s, fontsize=FS_LEG_H, fontweight='bold', va='top')

cN = COL_W - 1.0                                       # right-aligned point counts
def strata_block(y, title, strata, counts, letters):
    head(y, title); lg.text(cN, y, 'n', fontsize=FS_TXT, ha='right', va='top', fontweight='bold')
    y -= 6.0
    for k, (col, lab) in strata.items():
        lg.scatter([1.3], [y], s=DOT, c=col, edgecolors='white', lw=0.35, clip_on=False)
        lg.text(3.4, y, (f'{k}  ' if letters else '') + lab, fontsize=FS_TXT, va='center')
        lg.text(cN, y, f'{counts.get(k, 0)}', fontsize=FS_TXT, va='center', ha='right')
        y -= 3.9
    lg.plot([3.4, COL_W], [y + 1.95, y + 1.95], lw=0.4, color=INK)
    lg.text(3.4, y - 0.2, 'Total', fontsize=FS_TXT, va='center', fontweight='bold')
    lg.text(cN, y - 0.2, f'{sum(counts.values())}', fontsize=FS_TXT, va='center', ha='right', fontweight='bold')
    return y - 7.5
y = strata_block(y, 'Sampling strata, 2025', STRATA25, N25, True)
y = strata_block(y, 'Sampling strata, 2010', STRATA10, N10, False)

head(y, 'Elevation (m)'); y -= 5.4
cb_h = 2.6; cb_w = COL_W - 3.0
lg.imshow(TINT(np.linspace(0, 1, 256)[None, :])[..., :3], extent=(0, cb_w, y - cb_h, y), aspect='auto')
lg.add_patch(Rectangle((0, y - cb_h), cb_w, cb_h, fill=False, lw=0.4, edgecolor=INK))
for v in (40, 500, 1000, 1500, 2000):
    xx = (v - Z_LO) / (Z_HI - Z_LO) * cb_w
    lg.plot([xx, xx], [y - cb_h, y - cb_h - 0.8], lw=0.4, color=INK)
    lg.text(xx, y - cb_h - 1.1, f'{v:,}', fontsize=FS_TXT, ha='center', va='top')
y -= cb_h + 7.5

head(y, 'Boundaries'); y -= 5.6
for lab, kw in [('State', dict(lw=LW_STATE, color=INK, solid_capstyle='butt')),
                ('District', dict(lw=LW_DIST, color=GREY, solid_capstyle='butt')),
                ('Neighboring state', dict(lw=LW_THIN, color='#8c8c8c', ls='--')),
                ('International', dict(lw=LW_THIN, color='#8c8c8c'))]:
    lg.add_line(Line2D([0.3, 6.8], [y, y], **kw))
    lg.text(8.8, y, lab, fontsize=FS_TXT, va='center')
    y -= 4.3

print(f'{FIG_W} x {FIG_H:.1f} mm; maps {MW:.1f} x {MH:.1f} mm; locator {LOC_H:.1f} mm; '
      f'legend bottom y = {y:.1f} mm (map bottom {MAP_Y}); font {SERIF}')
fig.savefig(f'{OUT}.pdf', dpi=600)
fig.savefig(f'{OUT}_500dpi.png', dpi=500)
fig.savefig(f'{OUT}_500dpi.tif', dpi=500, pil_kwargs={'compression': 'tiff_lzw'})
