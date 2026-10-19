"""Figure 9 (v4): mapped jhum (2025) and fallow dynamics (1989-2025) at 30 m.

Layout: two state maps side by side and a right-hand column with the detail windows and legends.
Text in Times New Roman (Liberation Serif, metric-identical, where Times is not installed).
Nothing is printed on top of the state: titles sit above the frames, statistics and legends in
the column, and window labels outside the state with leader lines.

Usage: python plot_fig9.py <fig9_layers> <relief_layers_30m> <gee_exports_dir> <out_prefix>
"""
import sys
import numpy as np, geopandas as gpd, rasterio
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe
from rasterio.features import geometry_mask, shapes
from shapely.geometry import shape, box as sbox
from pyproj import Transformer

L9, LR, GEE, OUT = sys.argv[1:5]
CRS = 'EPSG:32646'; MM = 1 / 25.4

installed = {f.name for f in font_manager.fontManager.ttflist}
SERIF = 'Times New Roman' if 'Times New Roman' in installed else 'Liberation Serif'
mpl.rcParams.update({'font.family': SERIF, 'font.size': 9, 'axes.linewidth': 0.6,
                     'pdf.fonttype': 42, 'lines.dashed_pattern': (4, 1.5)})
FS_LETTER, FS_TITLE, FS_TEXT, FS_SMALL, FS_AXIS = 13, 10, 9, 8, 8
INK, GREY = '#1A1A1A', '#555555'
HALO_S = [pe.withStroke(linewidth=1.6, foreground='white', alpha=0.95)]
OUTLINE = '#00D8F0'          # cleared in the map year (2025)
OUTLINE_PREV = '#FFD21F'     # cleared the year before (fallow age 1): bare-looking, but not current-year jhum
JHUM_COL = np.array([0.66, 0.06, 0.09])
F_CMAP = LinearSegmentedColormap.from_list('fallow', ['#67000d', '#cb181d', '#fb6a4a', '#fdd0a2',
                                                      '#c7e9c0', '#41ab5d', '#00441b'])
F_MAX = 25

# ------------------------------------------------------------------ vectors
miz = gpd.read_file(f'{LR}/mizoram.gpkg'); STATE = miz.geometry.iloc[0]
dist = gpd.read_file(f'{LR}/districts.gpkg')
bgd = gpd.read_file(f'{LR}/bgd.gpkg').to_crs(CRS); mmr = gpd.read_file(f'{LR}/mmr.gpkg').to_crs(CRS)
ne = gpd.read_file(f'{LR}/ind_ne_states.gpkg').to_crs(CRS)
xmin, ymin, xmax, ymax = STATE.bounds; PAD = 5000
WIN = (xmin - PAD, xmax + PAD, ymin - PAD, ymax + PAD)
asp = (WIN[3] - WIN[2]) / (WIN[1] - WIN[0])

# ------------------------------------------------------------------ relief base (30 m)
with rasterio.open(f'{LR}/hillshade_30m.tif') as r:
    w = rasterio.windows.from_bounds(WIN[0], WIN[2], WIN[1], WIN[3], r.transform).round_offsets().round_lengths()
    shade = np.nan_to_num(r.read(1, window=w), nan=0.8); wtr = r.window_transform(w)
H, W = shade.shape
BASE_EXT = (wtr.c, wtr.c + W * 30, wtr.f - H * 30, wtr.f)
inside = ~geometry_mask([STATE], shade.shape, wtr)
base_out = (0.81 + 0.16 * shade)[..., None] * np.array([0.962, 0.960, 0.953])
base_in = (0.70 + 0.31 * shade)[..., None] * np.array([0.992, 0.987, 0.975])
BASE = (np.clip(np.where(inside[..., None], base_in, base_out), 0, 1) * 255).astype('uint8')
del shade

# ------------------------------------------------------------------ 30 m data layers
with rasterio.open(f'{GEE}/FIG9_jhum_2025.tif') as r:
    j30 = r.read(1); jb = r.bounds
J_RGBA = np.zeros(j30.shape + (4,), 'uint8'); J_RGBA[j30 == 1] = [*(JHUM_COL * 255).astype('uint8'), 255]
D_EXT = (jb.left, jb.right, jb.bottom, jb.top); del j30
with rasterio.open(f'{L9}/fallow_age_30m_clean.tif') as r:
    fa = r.read(1)
F_RGBA = (F_CMAP(np.clip(np.nan_to_num(fa), 0, F_MAX) / F_MAX) * 255).astype('uint8')
F_RGBA[..., 3] = np.where(np.isfinite(fa), 255, 0).astype('uint8'); del fa

# ------------------------------------------------------------------ detail windows
INSETS = []
up = lambda a: np.repeat(np.repeat(a, 8, axis=0), 8, axis=1)      # crisp pixels in any viewer
for key, lab in (('A', 'c1'), ('B', 'c2')):
    with rasterio.open(f'{GEE}/FIG9_inset_{key}_jhum.tif') as r:
        j = r.read(1); itr = r.transform; ib = r.bounds
    f = rasterio.open(f'{GEE}/FIG9_inset_{key}_fallow_age.tif').read(1).astype('float32')
    f[(f == 0) & (j == 0)] = np.nan                     # v1 nodata fix
    rgb = np.moveaxis(rasterio.open(f'{GEE}/FIG9_inset_{key}_landsat_rgb.tif').read().astype('float32'), 0, -1)
    lo, hi = np.percentile(rgb, 1, axis=(0, 1)), np.percentile(rgb, 99.5, axis=(0, 1))
    rgb = np.clip((rgb - lo) / (hi - lo), 0, 1) ** 0.85
    polys = [shape(g) for g, v in shapes(j.astype('uint8'), mask=j == 1, transform=itr) if v == 1]
    prev = (f == 1)                                      # last cleared in 2024
    polys_prev = [shape(g) for g, v in shapes(prev.astype('uint8'), mask=prev, transform=itr) if v == 1]
    INSETS.append(dict(lab=lab, rgb=up(rgb), fa=up(f), polys=polys, polys_prev=polys_prev,
                       ext=(ib.left, ib.right, ib.bottom, ib.top)))

# ------------------------------------------------------------------ layout (mm)
FIG_W = 190
M_L, M_R, M_T, M_B = 6.0, 5.0, 4.0, 4.0          # outer margins
LAT_LAB = 6.0                                     # room for latitude labels
MW = 57.5; MH = MW * asp
GAP_AB, GAP_C = 5.0, 8.0
TITLE_H, LON_LAB = 10.0, 6.0
FIG_H = M_T + TITLE_H + MH + LON_LAB + M_B
fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))
def ax_mm(x, y, w, h): return fig.add_axes([x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H])
def ftxt(x, y, s, **kw): return fig.text(x / FIG_W, y / FIG_H, s, **kw)

XA = M_L + LAT_LAB; XB = XA + MW + GAP_AB; MAP_Y = M_B + LON_LAB
axA = ax_mm(XA, MAP_Y, MW, MH); axB = ax_mm(XB, MAP_Y, MW, MH)
COL_X = XB + MW + GAP_C; COL_W = FIG_W - M_R - COL_X
TOP_Y = MAP_Y + MH                                # top of map frames

to_utm = Transformer.from_crs(4326, CRS, always_xy=True)
def graticule(ax, left_labels):
    for lon in np.arange(92.0, 94.01, 0.5):
        la = np.linspace(21.5, 25.0, 200); x, y = to_utm.transform(np.full_like(la, lon), la)
        ax.plot(x, y, lw=0.3, color='#8f8f8f', alpha=0.45, zorder=4, ls='--')
        xb = np.interp(WIN[2], y, x)
        if WIN[0] + 3000 < xb < WIN[1] - 3000:
            ax.annotate(f'{lon:.1f}°E', (xb, WIN[2]), xytext=(0, -3), textcoords='offset points',
                        ha='center', va='top', fontsize=FS_AXIS, color=GREY, annotation_clip=False)
            ax.plot([xb, xb], [WIN[2], WIN[2] + 1800], lw=0.6, color=INK, zorder=6)
    for lat in np.arange(22.0, 25.01, 0.5):
        lo = np.linspace(91.5, 94.0, 200); x, y = to_utm.transform(lo, np.full_like(lo, lat))
        ax.plot(x, y, lw=0.3, color='#8f8f8f', alpha=0.45, zorder=4, ls='--')
        yl = np.interp(WIN[0], x, y)
        if WIN[2] + 3000 < yl < WIN[3] - 3000:
            ax.plot([WIN[0], WIN[0] + 1800], [yl, yl], lw=0.6, color=INK, zorder=6)
            if left_labels:
                ax.annotate(f'{lat:.1f}°N', (WIN[0], yl), xytext=(-3, 0), textcoords='offset points',
                            ha='right', va='center', fontsize=FS_AXIS, color=GREY, rotation=90,
                            annotation_clip=False)

STATE_BUF = STATE.buffer(1500)
def safe_label(ax, name, lon, lat, rot):
    """Neighbour name, drawn only if its whole footprint stays clear of the state."""
    x, y = to_utm.transform(lon, lat)
    t = ax.text(x, y, name, fontsize=FS_SMALL, color='#7d7d7d', ha='center', va='center', rotation=rot,
                zorder=8, fontweight='bold', path_effects=[pe.withStroke(linewidth=1.5, foreground='#f4f4f2')])
    fig.canvas.draw()
    bb = t.get_window_extent().transformed(ax.transData.inverted())
    fp = sbox(bb.x0, bb.y0, bb.x1, bb.y1)
    if fp.intersects(STATE_BUF) or not sbox(*[WIN[0], WIN[2], WIN[1], WIN[3]]).contains(fp):
        t.remove(); print('  dropped label (would overlap):', name)

overlays = []
def map_panel(ax, x_mm, data_rgba, letter, title):
    ax.imshow(BASE, extent=BASE_EXT, interpolation='antialiased', zorder=1)
    overlays.append(ax.imshow(data_rgba, extent=D_EXT, interpolation='none', zorder=2))
    for g in (bgd, mmr): g.boundary.plot(ax=ax, lw=0.45, color='#8a8a8a', zorder=3)
    ne.boundary.plot(ax=ax, lw=0.4, color='#a0a0a0', ls='--', zorder=3)
    dist.boundary.plot(ax=ax, lw=0.45, color='#3d3d3d', alpha=0.45, zorder=5)
    miz.boundary.plot(ax=ax, lw=1.0, color=INK, zorder=7)
    graticule(ax, left_labels=(letter == 'a'))
    ax.set_xlim(WIN[0], WIN[1]); ax.set_ylim(WIN[2], WIN[3]); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(''); ax.set_ylabel('')
    for sp in ax.spines.values(): sp.set_linewidth(0.7); sp.set_color(INK)
    # title row above the frame: letter + title, no overlap with the map
    ftxt(x_mm, TOP_Y + 2.0, letter, fontsize=FS_LETTER, fontweight='bold', va='bottom', ha='left', color=INK)
    ftxt(x_mm + 5.0, TOP_Y + 2.3, title, fontsize=FS_TITLE, fontweight='bold', va='bottom', ha='left',
         color=INK)

map_panel(axA, XA, J_RGBA, 'a', 'Jhum clearing, 2025')
map_panel(axB, XB, F_RGBA, 'b', 'Years since last clearing')
for ax in (axA, axB):
    for name, lon, lat, rot in [('BANGLADESH', 92.27, 23.05, 90), ('MYANMAR', 93.38, 22.45, 90),
                                ('MANIPUR', 93.30, 24.45, 0), ('ASSAM', 92.55, 24.50, 0)]:
        safe_label(ax, name, lon, lat, rot)

# detail-window frames and leader-line labels placed outside the state (east side)
LABEL_X = WIN[1] - 2500
LABEL_LAT = {'c1': 22.93, 'c2': 22.67}          # both labels sit east of the state, on the Myanmar side
for d in INSETS:
    l, r_, b_, t = d['ext']
    axA.add_patch(Rectangle((l, b_), r_ - l, t - b_, fill=False, lw=0.9, edgecolor=INK, zorder=20))
    _, ly = to_utm.transform(93.40, LABEL_LAT[d['lab']])
    cy = (b_ + t) / 2
    elbow = LABEL_X - 9000
    xs, ys = ([r_, elbow, LABEL_X - 4200], [cy, ly, ly]) if abs(cy - ly) > 2000 else \
             ([r_, LABEL_X - 4200], [cy, cy])
    axA.plot(xs, ys, lw=0.5, color=INK, zorder=19)
    axA.text(LABEL_X, ly if abs(cy - ly) > 2000 else cy, d['lab'], fontsize=FS_TEXT, fontweight='bold',
             ha='right', va='center', color=INK, zorder=21, path_effects=HALO_S)

# scale bar and north arrow: panel b, south-east corner (outside the state)
def scalebar(ax, x0, y0, ticks, h, fs, lab_idx):
    for i in range(len(ticks) - 1):
        ax.add_patch(Rectangle((x0 + ticks[i] * 1000, y0), (ticks[i + 1] - ticks[i]) * 1000, h,
                     facecolor=INK if i % 2 == 0 else 'white', edgecolor=INK, lw=0.45, zorder=25))
    for i in lab_idx:
        last = i == len(ticks) - 1
        ax.text(x0 + ticks[i] * 1000, y0 + h * 1.8, f'{ticks[i]:g} km' if last else f'{ticks[i]:g}',
                fontsize=fs, ha='right' if last else 'center', va='bottom', color=INK, zorder=25,
                path_effects=HALO_S)
scalebar(axB, WIN[1] - 38000, WIN[2] + 4500, (0, 10, 20, 30), 1300, FS_SMALL, (0, 3))
axB.annotate('', xy=(WIN[1] - 6000, WIN[2] + 30000), xytext=(WIN[1] - 6000, WIN[2] + 19000),
             arrowprops=dict(arrowstyle='-|>,head_length=0.6,head_width=0.28', lw=0.9, color=INK), zorder=25)
axB.text(WIN[1] - 6000, WIN[2] + 31500, 'N', fontsize=FS_TEXT + 1, ha='center', va='bottom',
         fontweight='bold', color=INK, zorder=25, path_effects=HALO_S)

# ------------------------------------------------------------------ right column
ftxt(COL_X, TOP_Y + 2.0, 'c', fontsize=FS_LETTER, fontweight='bold', va='bottom', ha='left', color=INK)
ftxt(COL_X + 5.0, TOP_Y + 2.3, 'Detail windows (5 × 5 km)', fontsize=FS_TITLE, fontweight='bold',
     va='bottom', ha='left', color=INK)
TG = 2.5; TILE = (COL_W - TG) / 2
y = TOP_Y - 1.0
ftxt(COL_X + TILE / 2, y, 'Landsat 8/9', fontsize=FS_SMALL, ha='center', va='top', color=GREY)
ftxt(COL_X + TILE + TG + TILE / 2, y, 'Years since clearing', fontsize=FS_SMALL, ha='center', va='top', color=GREY)
y -= 4.5
for k, d in enumerate(INSETS):
    y -= TILE
    tA = ax_mm(COL_X, y, TILE, TILE); tB = ax_mm(COL_X + TILE + TG, y, TILE, TILE)
    tA.imshow(d['rgb'], extent=d['ext'], interpolation='none')
    if d['polys_prev']:
        gpd.GeoSeries(d['polys_prev'], crs=CRS).boundary.plot(ax=tA, lw=0.5, color=OUTLINE_PREV)
    gpd.GeoSeries(d['polys'], crs=CRS).boundary.plot(ax=tA, lw=0.6, color=OUTLINE)
    tB.set_facecolor('#efede6')
    tB.imshow(np.ma.masked_invalid(d['fa']), extent=d['ext'], cmap=F_CMAP, vmin=0, vmax=F_MAX,
              interpolation='none')
    for t in (tA, tB):
        t.set_xlim(d['ext'][:2]); t.set_ylim(d['ext'][2:]); t.set_xticks([]); t.set_yticks([])
        t.set_xlabel(''); t.set_ylabel('')
        for sp in t.spines.values(): sp.set_linewidth(0.7); sp.set_color(INK)
    # row label to the left of the tiles, in the gutter (not on the image)
    ftxt(COL_X - 1.2, y + TILE / 2, d['lab'], fontsize=FS_TEXT, fontweight='bold', ha='right',
         va='center', rotation=90, color=INK)
    y -= 3.0
# tile scale bar, below the tiles in its own row
sb = ax_mm(COL_X + TILE + TG, y - 2.5, TILE, 4.0); sb.set_xlim(0, 2.5); sb.set_ylim(0, 1); sb.axis('off')
for i, (a0, a1) in enumerate([(0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2)]):
    sb.add_patch(Rectangle((a0 + 0.05, 0.55), a1 - a0, 0.22, facecolor=INK if i % 2 == 0 else 'white',
                           edgecolor=INK, lw=0.45))
sb.text(0.05, 0.40, '0', fontsize=FS_SMALL, ha='center', va='top'); sb.text(2.05, 0.40, '2 km', fontsize=FS_SMALL,
        ha='center', va='top')
y -= 9.0

# legends (all text in the column, none on the maps)
def head(yy, s):
    ftxt(COL_X, yy, s, fontsize=FS_TEXT, fontweight='bold', va='top', color=INK)
lg = ax_mm(COL_X, M_B, COL_W, y - M_B); lg.set_xlim(0, COL_W); lg.set_ylim(M_B, y); lg.axis('off')

head(y, 'a  Jhum clearing, 2025'); y -= 5.6
lg.add_patch(Rectangle((0, y - 2.6), 5.5, 2.6, facecolor=JHUM_COL))
lg.text(7.5, y - 1.3, 'Mapped jhum (30 m)', fontsize=FS_SMALL, va='center')
y -= 5.0
lg.text(0, y, 'Mapped area 528.2 km²\nError-adjusted 578 km² (461–724)', fontsize=FS_SMALL, va='top',
        color='#333', linespacing=1.3)
y -= 11.5

head(y, 'b, c  Years since last clearing'); y -= 5.6
cb_h = 2.8
lg.imshow(F_CMAP(np.linspace(0, 1, 256)[None, :]), extent=(0, COL_W, y - cb_h, y), aspect='auto')
lg.add_patch(Rectangle((0, y - cb_h), COL_W, cb_h, fill=False, lw=0.4, edgecolor='#777'))
for v in range(0, 26, 5):
    xx = v / F_MAX * COL_W
    lg.plot([xx, xx], [y - cb_h, y - cb_h - 0.8], lw=0.4, color=INK)
    lg.text(xx, y - cb_h - 1.1, f'{v}' if v < 25 else '≥25', fontsize=FS_SMALL,
            ha='center' if 0 < v < 25 else ('left' if v == 0 else 'right'), va='top')
y -= cb_h + 6.2
lg.text(0, y, 'Blank: no clearing detected,\n1989–2025', fontsize=FS_SMALL, va='top', color=GREY,
        linespacing=1.3)
y -= 10.5

lg.plot([0, 5.5], [y - 1.3, y - 1.3], lw=1.3, color=OUTLINE)
lg.text(7.5, y - 1.3, 'Cleared in 2025 (mapped jhum)', fontsize=FS_SMALL, va='center')
y -= 4.6
lg.plot([0, 5.5], [y - 1.3, y - 1.3], lw=1.3, color=OUTLINE_PREV)
lg.text(7.5, y - 1.3, 'Cleared in 2024 (c)', fontsize=FS_SMALL, va='center')
y -= 6.0
lg.text(0, y, 'Panels a and b at full 30 m\nresolution in the digital version',
        fontsize=FS_SMALL - 0.5, va='top', color=GREY, linespacing=1.3, style='italic')

# ------------------------------------------------------------------ save
fig.savefig(f'{OUT}.pdf', dpi=400)        # relief at 400 dpi; 30 m data layers embedded unresampled
print(f'saved {OUT}: {FIG_W} x {FIG_H:.1f} mm; maps {MW} x {MH:.1f} mm; tiles {TILE:.1f} mm; font {SERIF}')
print(f'legend column bottom at y = {y:.1f} mm (must stay above {M_B})')
