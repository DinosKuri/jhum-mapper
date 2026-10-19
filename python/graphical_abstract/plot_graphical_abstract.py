"""Graphical abstract (Elsevier: landscape, >= 1328 x 531 px, legible at 13 x 5 cm).

Four steps, left to right: composites -> rule -> 30 m map -> validation and 37-year record.
Usage: python plot_graphical_abstract.py <figure_data_dir> <fig9_layers> <relief_layers> <gee_exports> <out_prefix>
"""
import sys
import numpy as np, pandas as pd, geopandas as gpd, rasterio
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon as MPoly, Circle
import matplotlib.patheffects as pe
from rasterio.features import geometry_mask, shapes
from shapely.geometry import shape
from scipy.stats import theilslopes

FD, L9, LR, GEE, OUT = sys.argv[1:6]
MM = 1 / 25.4
installed = {f.name for f in font_manager.fontManager.ttflist}
SERIF = 'Times New Roman' if 'Times New Roman' in installed else 'Liberation Serif'
mpl.rcParams.update({'font.family': SERIF, 'font.size': 11, 'axes.linewidth': 0.8, 'pdf.fonttype': 42,
                     'mathtext.fontset': 'custom', 'mathtext.rm': SERIF, 'mathtext.it': f'{SERIF}:italic',
                     'mathtext.bf': f'{SERIF}:bold'})
INK, GREY, MUTE = '#1A1A1A', '#5a5a5a', '#8a8a8a'
RED, BLUE, BG, CARD = '#C1292E', '#2E5E8C', '#F3F1EC', '#FFFFFF'
F_CMAP = LinearSegmentedColormap.from_list('fallow', ['#67000d', '#cb181d', '#fb6a4a', '#fdd0a2',
                                                      '#c7e9c0', '#41ab5d', '#00441b'])

# ------------------------------------------------------------------ canvas (mm)
W, H = 260.0, 108.0
fig = plt.figure(figsize=(W * MM, H * MM)); fig.patch.set_facecolor(BG)
def ax_mm(x, y, w, h, **kw): return fig.add_axes([x / W, y / H, w / W, h / H], **kw)
def T(x, y, s, **kw): return fig.text(x / W, y / H, s, **kw)

# title strip
T(8, H - 7.5, 'Annual mapping of shifting cultivation (jhum) from harmonized Landsat',
  fontsize=16, fontweight='bold', color=INK, va='center')
T(W - 8, H - 7.5, 'Mizoram, India · 1989–2025 · 30 m', fontsize=12.5, color=GREY, va='center', ha='right',
  style='italic')

TOP, BOT = H - 15.0, 7.0
CH = TOP - BOT
CARDS = [(8, 52), (66, 50), (122, 58), (186, 66)]     # x, width
STEPS = ['Symmetric anniversary\ncomposites', 'Two-variable\nspectral rule', '30 m annual maps',
         'Blind validation &\n37-year record']

def card(i):
    x, w = CARDS[i]
    fig.patches.append(FancyBboxPatch((x / W, BOT / H), w / W, CH / H, transform=fig.transFigure,
                       boxstyle='round,pad=0,rounding_size=0.012', facecolor=CARD, edgecolor='#d8d4cb',
                       lw=0.8, zorder=-5))
    cx, cy = x + 6.0, TOP - 6.5
    fig.add_artist(Circle((cx * MM, cy * MM), 3.4 * MM, transform=fig.dpi_scale_trans, facecolor=RED,
                          edgecolor='none', zorder=2))
    T(cx, cy, str(i + 1), fontsize=12, fontweight='bold', color='white', ha='center', va='center', zorder=3)
    T(x + 11.5, cy, STEPS[i], fontsize=11.5, fontweight='bold', color=INK, va='center', linespacing=1.0)

def chevron(x):
    yc = BOT + CH / 2
    fig.patches.append(MPoly([(x / W, (yc + 4) / H), ((x + 3.2) / W, yc / H), (x / W, (yc - 4) / H)],
                             closed=True, transform=fig.transFigure, facecolor=RED, edgecolor='none'))
for i in range(4): card(i)
for x in (61.4, 117.4, 181.4): chevron(x)

# ------------------------------------------------------------------ 1  composites
x0, w0 = CARDS[0]
a1 = ax_mm(x0 + 3, BOT + 19, w0 - 6, CH - 34); a1.set_xlim(0, 10); a1.set_ylim(0, 10); a1.axis('off')
MONTHS = 'NDJFMAM'
for yy, col, lab in [(7.6, '#4575b4', 'year $y$−1'), (3.2, RED, 'year $y$')]:
    a1.add_patch(FancyBboxPatch((0.3, yy - 0.75), 7.4, 1.5, boxstyle='round,pad=0,rounding_size=0.5',
                                facecolor=col, edgecolor='none'))
    for k, m in enumerate(MONTHS):
        a1.text(0.95 + k * 1.03, yy, m, color='white', fontsize=10, ha='center', va='center', fontweight='bold')
    a1.add_patch(Circle((8.85, yy), 0.95, facecolor='white', edgecolor=col, lw=1.4))
    a1.text(8.85, yy, 'P25', color=col, fontsize=9.5, ha='center', va='center', fontweight='bold')
    a1.text(0.3, yy + 1.15, f'Nov–May, {lab}', color=INK, fontsize=10, ha='left', va='bottom')
a1.annotate('', xy=(8.85, 4.25), xytext=(8.85, 6.55), arrowprops=dict(arrowstyle='-|>', lw=1.2, color=INK))
a1.text(7.75, 5.4, 'compare', fontsize=9.5, color=INK, ha='right', va='center', style='italic')
a1.text(5.0, 0.6, 'same months, same percentile:\nleaf-off cancels, burns remain', fontsize=9.5,
        color=GREY, ha='center', va='center', style='italic', linespacing=1.15)
T(x0 + w0 / 2, BOT + 13.5, r'$\mathrm{dNBR}=\mathrm{NBR_{pre}}-\mathrm{NBR_{post}}$', fontsize=12,
  color=INK, ha='center', va='center')
T(x0 + w0 / 2, BOT + 6.0, 'Landsat 5 · 7 · 8 · 9 (harmonized)', fontsize=10, color=GREY,
  ha='center', va='center')

# ------------------------------------------------------------------ 2  rule
x0, w0 = CARDS[1]
d = pd.read_csv(f'{FD}/fig3_featurespace_points.csv')
d = d[d.set == 'calibration'] if 'set' in d else d
a2 = ax_mm(x0 + 12, BOT + 16, w0 - 16, CH - 32)
a2.add_patch(Rectangle((-0.3, 0.20), 0.56, 1.2, facecolor=RED, alpha=0.10, edgecolor='none', zorder=0))
a2.plot([-0.3, 0.26, 0.26], [0.20, 0.20, 1.2], color=RED, lw=1.2, zorder=1)
nj = d[d.ref_jhum == 0]; jj = d[d.ref_jhum == 1]
a2.scatter(nj.post_wide25, nj.dnbr_wide25, s=6, c='#9aa5b1', lw=0, alpha=0.75, zorder=2)
a2.scatter(jj.post_wide25, jj.dnbr_wide25, s=8, c=RED, lw=0.25, edgecolors='white', zorder=3)
a2.set_xlim(-0.25, 0.80); a2.set_ylim(-0.15, 1.0)
a2.set_xlabel(r'$\mathrm{NBR_{post}}$', fontsize=10.5, labelpad=1)
a2.set_ylabel('dNBR', fontsize=10.5, labelpad=1)
a2.tick_params(labelsize=8.5, length=2, width=0.6, pad=1)
a2.set_xticks([0, 0.25, 0.5]); a2.set_yticks([0, 0.5, 1.0])
for s in ('top', 'right'): a2.spines[s].set_visible(False)
a2.text(0.30, 0.86, 'jhum', color=RED, fontsize=10.5, fontweight='bold', va='center')
a2.text(0.79, 0.72, 'not jhum', color='#6f7b88', fontsize=10.5, va='center', ha='right', fontweight='bold')
T(x0 + w0 / 2, BOT + 5.5, r'$\mathrm{NBR_{post}}\leq0.26$  and  $\mathrm{dNBR}\geq0.20$', fontsize=10,
  color=INK, ha='center', va='center')

# ------------------------------------------------------------------ 3  map
x0, w0 = CARDS[2]
miz = gpd.read_file(f'{LR}/mizoram.gpkg').to_crs('EPSG:32646'); ST = miz.geometry.iloc[0]
b = ST.bounds; pad = 3000
WIN = (b[0] - pad, b[2] + pad, b[1] - pad, b[3] + pad)
with rasterio.open(f'{L9}/fallow_age_30m_clean.tif') as r:
    fa = r.read(1, out_shape=(r.height // 2, r.width // 2), resampling=rasterio.enums.Resampling.nearest)
    fb = r.bounds
FR = (F_CMAP(np.clip(np.nan_to_num(fa), 0, 25) / 25) * 255).astype('uint8')
FR[..., 3] = np.where(np.isfinite(fa), 255, 0); del fa
with rasterio.open(f'{LR}/hillshade_25m.tif') as r:
    sh = r.read(1, out_shape=(r.height // 4, r.width // 4)); sb = r.bounds
    shtr = r.transform * r.transform.scale(4, 4)
ins = ~geometry_mask([ST], sh.shape, shtr)
base = (0.80 + 0.18 * np.nan_to_num(sh, nan=0.8))[..., None] * np.array([0.97, 0.965, 0.955])
base = np.where(ins[..., None], np.clip(base * 1.03, 0, 1), np.clip(base * 0.985, 0, 1))
mh = CH - 25; mw = mh * (WIN[1] - WIN[0]) / (WIN[3] - WIN[2])
a3 = ax_mm(x0 + 3.0, BOT + 11.0, mw, mh)
a3.imshow(base, extent=(sb.left, sb.right, sb.bottom, sb.top), interpolation='antialiased')
a3.imshow(FR, extent=(fb.left, fb.right, fb.bottom, fb.top), interpolation='antialiased')
miz.boundary.plot(ax=a3, lw=0.8, color=INK)
a3.set_xlim(WIN[0], WIN[1]); a3.set_ylim(WIN[2], WIN[3]); a3.set_aspect('equal'); a3.axis('off')
# detail tile (Landsat with mapped jhum outlined), linked to its location
with rasterio.open(f'{GEE}/FIG9_inset_A_jhum.tif') as r:
    j = r.read(1); itr = r.transform; ib = r.bounds
rgb = np.moveaxis(rasterio.open(f'{GEE}/FIG9_inset_A_landsat_rgb.tif').read().astype('float32'), 0, -1)
lo, hi = np.percentile(rgb, 1, axis=(0, 1)), np.percentile(rgb, 99.5, axis=(0, 1))
rgb = np.clip((rgb - lo) / (hi - lo), 0, 1) ** 0.85
rgb = np.repeat(np.repeat(rgb, 6, 0), 6, 1)
polys = [shape(g) for g, v in shapes(j.astype('uint8'), mask=j == 1, transform=itr) if v == 1]
tw = w0 - mw - 9.0
tx, ty = x0 + 6.0 + mw, TOP - 16.0 - tw
a3t = ax_mm(tx, ty, tw, tw)
a3t.imshow(rgb, extent=(ib.left, ib.right, ib.bottom, ib.top), interpolation='none')
gpd.GeoSeries(polys, crs='EPSG:32646').boundary.plot(ax=a3t, lw=0.7, color='#00D8F0')
a3t.set_xlim(ib.left, ib.right); a3t.set_ylim(ib.bottom, ib.top); a3t.set_xticks([]); a3t.set_yticks([])
a3t.set_xlabel(''); a3t.set_ylabel('')
for s in a3t.spines.values(): s.set_color(INK); s.set_linewidth(0.8)
a3.add_patch(Rectangle((ib.left, ib.bottom), ib.right - ib.left, ib.top - ib.bottom, fill=False, lw=0.8,
                       edgecolor=INK, zorder=5))
T(tx + tw / 2, ty - 1.2, '5 km detail,\njhum outlined', fontsize=9.5, color=GREY, ha='center', va='top',
  linespacing=1.1)
# fallow-age key
kx, ky, kw_ = tx, BOT + 20.0, tw
ka = ax_mm(kx, ky, kw_, 2.6)
ka.imshow(F_CMAP(np.linspace(0, 1, 128)[None, :]), aspect='auto', extent=(0, 25, 0, 1))
ka.set_yticks([]); ka.set_xticks([0, 25]); ka.set_xticklabels(['0', '25+'])
ka.tick_params(labelsize=8.5, length=1.5, pad=1)
T(kx, ky + 3.6, 'Years since\nlast clearing', fontsize=9.5, color=INK, va='bottom', linespacing=1.05)
T(x0 + w0 / 2, BOT + 5.0, '56 % of the state cleared\nat least once since 1989', fontsize=9.5,
  color=GREY, ha='center', va='center', style='italic', linespacing=1.1)

# ------------------------------------------------------------------ 4  validation + series
x0, w0 = CARDS[3]
v = pd.read_csv(f'{FD}/fig10_two_year_validation.csv')
# the file also holds the superseded 2025 three-stratum round, kept there for Fig 5;
# the abstract shows only the two published assessments
v = v[~((v.year == 2025) & (v.design == 'three strata'))].reset_index(drop=True)
for k, row in enumerate(v.itertuples()):
    bx = x0 + 5 + k * 33.5
    by = TOP - 30.5
    fig.patches.append(FancyBboxPatch((bx / W, by / H), 31 / W, 15.5 / H, transform=fig.transFigure,
                       boxstyle='round,pad=0,rounding_size=0.01', facecolor='#F8F1F0' if k else '#EEF2F7',
                       edgecolor='none', zorder=-4))
    T(bx + 15.5, by + 12.0, f'{row.year}  ·  {row.sensors}', fontsize=9.5, color=GREY, ha='center', va='center')
    T(bx + 15.5, by + 6.6, f'F1 = {row.F1:.2f}', fontsize=15, fontweight='bold',
      color=BLUE if k == 0 else RED, ha='center', va='center')
    T(bx + 15.5, by + 2.3, f'{row.n_points} blind points', fontsize=9, color=GREY, ha='center', va='center')

s = pd.read_csv(f'{FD}/fig4_annual_series.csv')
full = s[s.year >= 1998]
a4 = ax_mm(x0 + 12.5, BOT + 8.0, w0 - 17, 33)
a4.fill_between(s.year, 0, s.primary_km2, color=RED, alpha=0.10, lw=0)
a4.plot(s.year, s.primary_km2, color=RED, lw=1.3)
a4.scatter(s.year, s.primary_km2, s=7, color=RED, zorder=3, lw=0)
sl, ic, _, _ = theilslopes(full.primary_km2, full.year)
xx = np.array([1998, 2025]); a4.plot(xx, ic + sl * xx, color=INK, lw=1.0, ls=(0, (4, 2)))
for (y0, y1), m in [((1998, 2009), 975), ((2020, 2025), 443)]:
    a4.plot([y0, y1], [m, m], color=INK, lw=2.2, solid_capstyle='butt', alpha=0.85)
a4.text(2003.5, 1730, '975 km²', fontsize=10.5, fontweight='bold', ha='center', va='center', color=INK)
a4.text(2022.5, 800, '443 km²', fontsize=10.5, fontweight='bold', ha='center', va='center', color=INK)
T(x0 + w0 / 2, BOT + 44.5, f'Trend 1998–2025: −{abs(sl):.1f} km$^2$ yr$^{{-1}}$ ($p$ < 0.001)', fontsize=9.5,
  color=GREY, ha='center', va='center')
a4.set_xlim(1988, 2026); a4.set_ylim(0, 1850)
a4.set_xticks([1990, 2000, 2010, 2020]); a4.set_yticks([0, 500, 1000, 1500])
a4.tick_params(labelsize=8.5, length=2, width=0.6, pad=1)
a4.set_ylabel('Mapped jhum (km$^2$)', fontsize=10, labelpad=1)
for sp in ('top', 'right'): a4.spines[sp].set_visible(False)

fig.savefig(f'{OUT}.pdf', dpi=600, facecolor=BG)
fig.savefig(f'{OUT}.png', dpi=600, facecolor=BG)
print(f'saved {OUT}: {W} x {H} mm; Sen slope 1998-2025 = {sl:.2f}')
