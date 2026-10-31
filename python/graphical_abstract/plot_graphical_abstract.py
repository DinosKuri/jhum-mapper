"""Graphical abstract, Remote Sensing of Environment.

Canvas 132.1 x 52.8 mm (2.5:1), rendered at 600 dpi -> 3120 x 1248 px RGB TIFF.
Every number is read from the archived figure data or computed from the 2025 map.

Panels, left to right (the order of the title):
  1  symmetric anniversary composites: separation ratio of three anniversary designs
     at the 400 calibration points
  2  dual-sensor validation strata: producer's accuracy of the same 2025 Landsat map
     under the two stratifications, with 95 % intervals
  3  37 years of Landsat: the 2025 map (share of each 2 km cell cleared, area-averaged)
     and the 1989-2025 mapped-area series with the two error-adjusted estimates

Usage: python plot_graphical_abstract.py <figure_data_dir> <jhum_2025.tif> <layers_dir> <out_prefix>
  <layers_dir> must hold mizoram.gpkg (EPSG:32646).
"""
import sys, os
import numpy as np, pandas as pd, geopandas as gpd, rasterio
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.stats import theilslopes
from PIL import Image

FD, JHUM, LAY, OUT = sys.argv[1:5]
MM = 1 / 25.4
installed = {f.name for f in font_manager.fontManager.ttflist}
SERIF = 'Times New Roman' if 'Times New Roman' in installed else 'Liberation Serif'   # metric-equivalent
mpl.rcParams.update({
    'font.family': SERIF, 'font.size': 8, 'axes.linewidth': 0.5, 'pdf.fonttype': 42,
    'mathtext.fontset': 'custom', 'mathtext.rm': SERIF, 'mathtext.it': f'{SERIF}:italic', 'mathtext.bf': f'{SERIF}:bold',
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5, 'axes.edgecolor': '#555555', 'xtick.color': '#555555',
    'ytick.color': '#555555', 'xtick.labelcolor': '#1A1A1A', 'ytick.labelcolor': '#1A1A1A'})
INK, SEC = '#1A1A1A', '#555555'
RED, BLUE, BLUE_T = '#C1272D', '#2C5F8E', '#9FB4CC'
LAND, BAND, CHEV = '#E3E3E3', '#EFEFEF', '#9AA5B1'
HEAD, BIG, BASE, NOTE = 10.5, 16, 8, 7

W, H = 132.1, 52.8
fig = plt.figure(figsize=(W * MM, H * MM)); fig.patch.set_facecolor('white')
def ax_mm(x, y, w, h): return fig.add_axes([x / W, y / H, w / W, h / H])
def T(x, y, s, **kw): return fig.text(x / W, y / H, s, **kw)
def clean(ax, keep=('bottom',)):
    for sp in ('top', 'right', 'left', 'bottom'):
        ax.spines[sp].set_visible(sp in keep)
    ax.set_facecolor('none')

# panel frame: x positions (mm)
P1, P2, P3 = (2.0, 29.5), (35.0, 36.0), (75.0, 55.1)
TOP = 49.6                                             # header top
heads = (('Compositing', 'symmetry'), ('Dual-sensor', 'stratification'), ('37 years of Landsat,', '1989–2025'))
for (x, w), (l1, l2) in zip((P1, P2, P3), heads):
    T(x, TOP, l1 + '\n' + l2, fontsize=HEAD, fontweight='bold', color=INK, va='top', linespacing=0.95)
for xc in (P1[0] + P1[1] + 1.8, P2[0] + P2[1] + 2.0):     # chevrons mean only "leads to"
    T(xc, 25, '›', fontsize=17, color=CHEV, ha='center', va='center')
def mark(x, y, **kw):                                  # legend glyph drawn as a marker, in mm
    fig.add_artist(Line2D([x / W], [y / H], transform=fig.transFigure, **kw))
def dash(x, y, w=4.0, **kw):
    fig.add_artist(Line2D([x / W, (x + w) / W], [y / H, y / H], transform=fig.transFigure, **kw))

# ------------------------------------------------------------------ panel 1
x0, w0 = P1
months = ['N', 'D', 'J', 'F', 'M', 'A', 'M']
sw, mx = 2.55, x0 + 6.6
for lab, col, yy in (('$t$−1', BLUE, 36.0), ('$t$', RED, 32.0)):
    for i, m in enumerate(months):
        fig.patches.append(Rectangle(((mx + i * sw) / W, yy / H), (sw - 0.3) / W, 3.0 / H,
                                     transform=fig.transFigure, fc=col, ec='none'))
        T(mx + i * sw + (sw - 0.3) / 2, yy + 1.5, m, fontsize=6.5, color='white', ha='center', va='center', fontweight='bold')
    T(x0, yy + 1.5, lab, fontsize=BASE + 1, color=col, va='center')
    T(mx + 7 * sw + 0.5, yy + 1.5, 'P25', fontsize=BASE, color=INK, va='center', fontweight='bold')
c = pd.read_csv(f'{FD}/fig3_compositing.csv').set_index('order')
designs = [(c.loc[3], '75/25, Mar–May', BLUE_T), (c.loc[4], '25/25, Mar–May', BLUE_T), (c.loc[5], '25/25, Nov–May', BLUE)]
a1 = ax_mm(x0, 7.0, w0, 22.0); clean(a1, keep=())
for i, (row, lab, col) in enumerate(designs):
    yy = 2 - i
    a1.barh(yy, row.separation_ratio, height=0.36, color=col)
    a1.text(0, yy + 0.25, lab, fontsize=BASE, color=INK, va='bottom', ha='left')
    big = i == 2
    a1.text(row.separation_ratio + 0.5, yy - 0.02, f'{row.separation_ratio:.1f}', va='center', ha='left',
            fontsize=BIG if big else BASE + 1, fontweight='bold', color=BLUE if big else INK)
a1.set_xlim(0, 25.5); a1.set_ylim(-0.3, 2.75); a1.set_xticks([]); a1.set_yticks([])
T(x0, 5.4, 'jhum ÷ background dNBR,', fontsize=BASE, color=SEC, va='center')
T(x0, 2.6, '400 calibration points', fontsize=BASE, color=SEC, va='center')

# ------------------------------------------------------------------ panel 2
x0, w0 = P2
v = pd.read_csv(f'{FD}/fig4_validation_designs.csv')
old = v[(v.year == 2025) & (v.design == 'three strata')].iloc[0]
new = v[(v.year == 2025) & (v.design.str.startswith('five'))].iloc[0]
a2 = ax_mm(x0, 22.0, w0, 18.4); clean(a2)
rows = [('Landsat-only strata, $n$ = 310', old, False), ('Landsat + Sentinel-2, $n$ = 520', new, True)]
for i, (lab, r, filled) in enumerate(rows):
    yy = 1 - i
    a2.plot([r.PA_lo, r.PA_hi], [yy, yy], color=BLUE, lw=1.0, solid_capstyle='butt')
    a2.plot(r.PA, yy, 'o', ms=5.8, mfc=BLUE if filled else 'white', mec=BLUE, mew=1.1, zorder=4)
    a2.text(0.4, yy + 0.22, lab, fontsize=BASE, color=INK, va='bottom', ha='left')
    if filled:
        a2.text(r.PA_hi + 0.02, yy, f'{r.PA:.3f}', va='center', ha='left', fontsize=BIG, fontweight='bold', color=BLUE)
    else:
        a2.text(r.PA_lo - 0.02, yy, f'{r.PA:.3f}', va='center', ha='right', fontsize=BASE + 1, fontweight='bold', color=INK)
a2.set_xlim(0.4, 1.0); a2.set_ylim(-0.35, 1.75); a2.set_yticks([])
a2.set_xticks([0.4, 0.6, 0.8, 1.0]); a2.set_xticklabels(['0.4', '0.6', '0.8', '1.0']); a2.tick_params(labelsize=BASE)
T(x0, 16.6, "Producer's accuracy (95 % CI)", fontsize=BASE, color=SEC, va='center')
T(x0, 11.6, 'Omissions seen: 1/240 vs 32/435', fontsize=BASE, color=INK, va='center')
T(x0, 8.4, f"User's accuracy {new.UA:.3f}", fontsize=BASE, color=INK, va='center')
T(x0, 5.3, f'(95 % CI {new.UA_lo:.3f}–{new.UA_hi:.3f})', fontsize=BASE, color=SEC, va='center')

# ------------------------------------------------------------------ panel 3: map
x0, w0 = P3
miz = gpd.read_file(f'{LAY}/mizoram.gpkg').to_crs('EPSG:32646')
AGG = 67                                               # 30 m -> 2.01 km cells, area-averaged share
with rasterio.open(JHUM) as r:
    a = r.read(1); hh, ww = (r.height // AGG) * AGG, (r.width // AGG) * AGG
    share = (a[:hh, :ww] == 1).reshape(hh // AGG, AGG, ww // AGG, AGG).mean(axis=(1, 3)) * 100
    bnd = r.bounds; ext = (bnd.left, bnd.left + ww * 30, bnd.top - hh * 30, bnd.top)
from rasterio.features import geometry_mask
from affine import Affine
tr = Affine(AGG * 30, 0, bnd.left, 0, -AGG * 30, bnd.top)
inside = ~geometry_mask(miz.geometry, out_shape=share.shape, transform=tr)
mx0, my0, mx1, my1 = miz.total_bounds
mh = 33.0; mw = mh * (mx1 - mx0) / (my1 - my0)
am = ax_mm(x0, 3.4, mw, mh); am.set_axis_off()
miz.plot(ax=am, color=LAND, lw=0, zorder=1)
cmap = mpl.colormaps['Reds']
am.imshow(np.ma.masked_where(~inside | (share <= 0), share), extent=ext, cmap=cmap, vmin=0, vmax=20,
          interpolation='nearest', zorder=2)
miz.boundary.plot(ax=am, lw=0.6, color=SEC, zorder=3)
am.set_xlim(mx0 - 1000, mx1 + 1000); am.set_ylim(my0 - 1000, my1 + 1000); am.set_aspect('equal')
sx, sy = mx0 - 800, my0 + 2500
am.plot([sx, sx + 50000], [sy, sy], color=INK, lw=1.3, solid_capstyle='butt', zorder=5)
am.text(sx + 25000, sy + 2500, '50 km', ha='center', va='bottom', fontsize=BASE - 1, zorder=5)
nx, ny = mx0 + 6000, my1 - 36000
am.annotate('', xy=(nx, ny + 17000), xytext=(nx, ny), arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.6, mutation_scale=6), zorder=5)
am.text(nx, ny + 18500, 'N', ha='center', va='bottom', fontsize=BASE - 1, fontweight='bold')
T(x0, 42.0, 'Cleared, 2025', fontsize=BASE, color=RED, va='center', fontweight='bold')
cb = ax_mm(x0, 38.6, mw, 1.3)
mpl.colorbar.ColorbarBase(cb, cmap=cmap, norm=mpl.colors.Normalize(0, 20), orientation='horizontal', ticks=[0, 10, 20])
cb.set_xticklabels(['0', '10', '20 %']); cb.tick_params(labelsize=NOTE, length=1.2, pad=0.6)
[sp.set_linewidth(0.4) for sp in cb.spines.values()]

# ------------------------------------------------------------------ panel 3: series
s = pd.read_csv(f'{FD}/fig6_annual_series.csv')
gx = x0 + mw + 7.0; gw = W - 1.5 - gx
ag = ax_mm(gx, 17.2, gw, 19.6); clean(ag, keep=('bottom', 'left'))
ag.axvspan(1988.5, 1997.5, color=BAND, lw=0, zorder=0)
ag.plot(s.year, s.primary_km2, color=RED, lw=1.1, zorder=3)
full = s[s.year >= 1998]
for y0, y1, lab in ((1998, 2009, '975'), (2020, 2025, '443')):
    m = full[(full.year >= y0) & (full.year <= y1)].primary_km2.mean()
    assert f'{m:.0f}' == lab, (lab, m)
    ag.plot([y0, y1], [m, m], color=BLUE, lw=2.4, solid_capstyle='butt', zorder=4)
import matplotlib.patheffects as pe
HALO = [pe.withStroke(linewidth=2.2, foreground='white')]
ag.text(2003.5, 1060, '975', ha='center', va='bottom', fontsize=BASE + 1, fontweight='bold', color=BLUE, path_effects=HALO, zorder=7)
ag.text(2025.5, 760, '443', ha='right', va='bottom', fontsize=13, fontweight='bold', color=BLUE, path_effects=HALO, zorder=7)
sen = theilslopes(full.primary_km2, full.year)
xx = np.array([1998, 2025]); ag.plot(xx, sen.intercept + sen.slope * xx, color=INK, lw=0.8, ls=(0, (3, 2)), zorder=5)
ea = v.set_index(['year', 'design'])
for yr, des in ((2010, 'three strata'), (2025, 'five strata; both sensors')):
    r = ea.loc[(yr, des)]
    ag.errorbar(yr, r.adjusted_km2, yerr=[[r.adjusted_km2 - r.adjusted_lo], [r.adjusted_hi - r.adjusted_km2]], fmt='D',
                ms=3.4, color=INK, mfc=INK, elinewidth=1.0, capsize=0, zorder=6)
ag.set_xlim(1988.5, 2026); ag.set_ylim(0, 1600); ag.set_yticks([0, 500, 1000, 1500]); ag.set_yticklabels(['0', '500', '1000', '1500'])
ag.set_xticks([1990, 2005, 2020]); ag.tick_params(labelsize=BASE)
T(W - 1.5, 42.0, 'Mapped jhum (km$^2$)', fontsize=BASE, color=RED, va='center', ha='right', fontweight='bold')
lx = x0 + mw + 1.2
dash(lx, 11.6, 3.6, color=INK, lw=0.8, ls=(0, (3, 2)))
T(lx + 4.6, 11.6, 'Sen ' + f'{sen.slope:.1f}'.replace('-', '−') + ' km$^2$ yr$^{-1}$', fontsize=BASE, color=INK, va='center')
mark(lx + 1.8, 8.4, marker='D', ms=3.4, color=INK, ls='none')
T(lx + 4.6, 8.4, 'error-adjusted area, 95 % CI', fontsize=BASE, color=INK, va='center')
T(lx + 4.6, 4.9, r'2025: $\mathbf{578}$ (461–724) km$^2$', fontsize=BASE + 1, color=INK, va='center')
T(W - 1.5, 0.9, 'Map lines delineate study areas and do not necessarily depict accepted national boundaries.',
  fontsize=NOTE, color=SEC, ha='right', va='center')

# ------------------------------------------------------------------ export
fig.savefig(f'{OUT}.pdf', facecolor='white', metadata={'Title': 'Graphical abstract', 'Author': 'L. Khiangte'})
tmp = f'{OUT}_tmp.png'
fig.savefig(tmp, dpi=600, facecolor='white')
im_ = Image.open(tmp).convert('RGB'); os.remove(tmp)
if im_.size != (3120, 1248): im_ = im_.resize((3120, 1248), Image.LANCZOS)
im_.save(f'{OUT}.tif', dpi=(600, 600), compression='tiff_lzw')
im_.resize((500, 200), Image.LANCZOS).save(f'{OUT}_500x200_check.png')
print('saved', OUT, im_.size, f'Sen {sen.slope:.2f}')
