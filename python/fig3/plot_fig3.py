"""Figure 3 (final): visual validation, PlanetScope / Sentinel-2 / Landsat 8/9, four sites side by side.

Every tile is cut from the georeferenced source to the same 1.24 x 1.24 km square (UTM 46N) centred
on the validation point, so all tiles of a site show identical ground. All text in Times New Roman
(Liberation Serif where Times is absent), added at final size. Scale: one scale bar per sensor row.
Usage: python plot_fig3_final.py <tiles_dir> <landsat_dir> <out_prefix>
"""
import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, rasterio
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe
from rasterio.warp import Resampling
from chips import window, chip_rgb

TD, LD, OUT = sys.argv[1:4]
installed = {f.name for f in font_manager.fontManager.ttflist}
SERIF = 'Times New Roman' if 'Times New Roman' in installed else 'Liberation Serif'
mpl.rcParams.update({'font.family': SERIF, 'pdf.fonttype': 42, 'axes.linewidth': 0.5,
                     'mathtext.fontset': 'custom', 'mathtext.rm': SERIF})
MM = 1 / 25.4; INK = 'black'
HALF = 620; PX = 520                                   # 1.24 km window; chip pixels (> 600 dpi at tile size)
HALO = [pe.withStroke(linewidth=1.6, foreground='black')]

# (label, file) per tile; Planet DEC_25 / NOV_25 / SEP_25 files are the 2024 basemaps (verified against imagery)
SITES = {
 'a': dict(pid='V2025_051', lon=93.277870, lat=23.323005, case='Detected, wide branch',
           info=(0.170, 0.207, '(narrow dNBR 0.082)'),
           P=[('Nov 2024', 'V2025_51_CONFIRMJHUM/NOV_2024'), ('Dec 2024', 'V2025_51_CONFIRMJHUM/DEC_2024'),
              ('Jan 2025', 'V2025_51_CONFIRMJHUM/JAN_2025'), ('Mar 2025', 'V2025_51_CONFIRMJHUM/MARCH_2025'),
              ('May 2025', 'V2025_51_CONFIRMJHUM/MAY_2025'), ('Sep 2025', 'V2025_51_CONFIRMJHUM/SEP_2025')],
           S=[('Nov 2024', 'V2025_51_CONFIRMJHUM_sentinnel/nov_24'), ('Dec 2024', 'V2025_51_CONFIRMJHUM_sentinnel/dec_24'),
              ('Jan 2025', 'V2025_51_CONFIRMJHUM_sentinnel/jan_25'), ('Mar 2025', 'V2025_51_CONFIRMJHUM_sentinnel/march_25'),
              ('May 2025', 'V2025_51_CONFIRMJHUM_sentinnel/may_25'), ('Sep 2025', 'V2025_51_CONFIRMJHUM_sentinnel/sept_25')],
           L=['NOV2024', 'DEC2024', 'JAN2025', 'MAR2025', 'MAY2025', 'SEP2025'], lfile='a_V2025_051'),
 'b': dict(pid='V2025_030', lon=92.909740, lat=22.668672, case='Detected, narrow branch only',
           info=(0.388, 0.147, '(narrow dNBR 0.577)'),
           P=[('Nov 2024', 'V2025_030/NOV_24'), ('Dec 2024', 'V2025_030/DEC_25'), ('Feb 2025', 'V2025_030/FEB_25'),
              ('Mar 2025', 'V2025_030/MARCH'), ('Apr 2025', 'V2025_030/APR_25'), ('Sep 2025', 'V2025_030/sept')],
           S=[('Nov 2024', 'V2025_030_sentinnel/nov_24'), ('Dec 2024', 'V2025_030_sentinnel/dec_24'),
              ('Jan 2025', 'V2025_030_sentinnel/jan_25'), ('Mar 2025', 'V2025_030_sentinnel/mar_25'),
              ('May 2025', 'V2025_030_sentinnel/may_25'), ('Aug 2025', 'V2025_030_sentinnel/aug_25')],
           L=['NOV2024', 'DEC2024', 'JAN2025', 'MAR2025', 'MAY2025', 'AUG2025'], lfile='b_V2025_030'),
 'c': dict(pid='V2025_015', lon=92.666387, lat=23.851484, case='Commission error, deciduous leaf-off',
           info=(0.242, 0.207, '(threshold 0.20)'),
           P=[('Nov 2024', 'V2025_015/NOV_24'), ('Dec 2024', 'V2025_015/DEC_25'), ('Jan 2025', 'V2025_015/JAN_25'),
              ('Mar 2025', 'V2025_015/MARCH_25'), ('May 2025', 'V2025_015/MAY_25'), ('Aug 2025', 'V2025_015/AUG_25')],
           S=[('Nov 2024', 'V2025_015_sentinnel/Nov_24'), ('Dec 2024', 'V2025_015_sentinnel/dec_24'),
              ('Jan 2025', 'V2025_015_sentinnel/jan_25'), ('Mar 2025', 'V2025_015_sentinnel/mar_25'),
              ('May 2025', 'V2025_015_sentinnel/may_25'), ('Aug 2025', 'V2025_015_sentinnel/aug_25')],
           L=['NOV2024', 'DEC2024', 'JAN2025', 'MAR2025', 'MAY2025', 'AUG2025'], lfile='c_V2025_015'),
 'd': dict(pid='V2025_100', lon=92.586077, lat=22.451729, case='Omission, young fallow',
           info=(0.313, 0.119, '(threshold 0.20)'),
           P=[('Sep 2024', 'V2025_100/SEP_25'), ('Nov 2024', 'V2025_100/NOV_25'), ('Jan 2025', 'V2025_100/JAN_25'),
              ('Mar 2025', 'V2025_100/MARCH_25'), ('May 2025', 'V2025_100/MAY_25'), ('Aug 2025', 'V2025_100/AUG_25')],
           S=[('Sep 2024', 'V2025_100_sentinnel/sept_24'), ('Nov 2024', 'V2025_100_sentinnel/nov_24'),
              ('Jan 2025', 'V2025_100_sentinnel/jan_25'), ('Mar 2025', 'V2025_100_sentinnel/march_24'),
              ('May 2025', 'V2025_100_sentinnel/may_25'), ('Aug 2025', 'V2025_100_sentinnel/aug_25')],
           L=['SEP2024', 'NOV2024', 'JAN2025', 'MAR2025', 'MAY2025', 'AUG2025'], lfile='d_V2025_100'),
}
LAB = {'SEP2024': 'Sep 2024', 'NOV2024': 'Nov 2024', 'DEC2024': 'Dec 2024', 'JAN2025': 'Jan 2025',
       'MAR2025': 'Mar 2025', 'MAY2025': 'May 2025', 'AUG2025': 'Aug 2025', 'SEP2025': 'Sep 2025'}

# ---------------------------------------------------------------- build all chips
import glob
chips = {}
for k, s in SITES.items():
    W = window(s['lon'], s['lat'], HALF); s['W'] = W
    for sens in ('P', 'S'):
        for i, (lab, f) in enumerate(s[sens]):
            rgb, tr = chip_rgb(f'{TD}/{k}/{f}.tif', W, PX)
            chips[k, sens, i] = (lab, rgb); s['tr'] = tr
    lf = glob.glob(f'{LD}/*FIG3_landsat_{s["lfile"]}*.tif')[0]
    r = rasterio.open(lf); idx = {d: j + 1 for j, d in enumerate(r.descriptions)}
    for i, m in enumerate(s['L']):
        rgb, _ = chip_rgb(lf, W, PX, bands=[idx[f'{m}_{b}'] for b in 'RGB'], resampling=Resampling.nearest,
                          scale=10000.0, src_nodata=-9999)
        chips[k, 'L', i] = (LAB[m], rgb)

# uniform stretch per sensor (same for all four sites): Planet basemaps are already display-ready
def pooled(sens, cloud_cut):
    v = np.concatenate([chips[k, sens, i][1][::4, ::4].reshape(-1, 3) for k in SITES for i in range(6)])
    v = v[np.isfinite(v).all(1) & (v.max(1) < cloud_cut)]             # leave clouds out of the statistics
    return np.percentile(v, 1, axis=0), np.percentile(v, 99.5, axis=0)
STRETCH = {'P': (np.zeros(3), np.ones(3))}                             # Planet basemaps: already display-ready
lo, hi = pooled('S', 0.6)
STRETCH['S'] = (lo, hi)                                                # per band, pooled over all sites: removes the blue cast
lo, hi = pooled('L', 0.25)
STRETCH['L'] = (np.zeros(3), np.full(3, 1.1 * hi[0]))                  # one value for all bands: natural colour
GAMMA = {'P': 1.0, 'S': 1.0, 'L': 1.3}
def show(sens, rgb):
    lo, hi = STRETCH[sens]
    a = np.clip((rgb - lo) / (hi - lo), 0, 1) ** (1 / GAMMA[sens])
    nod = ~np.isfinite(rgb).all(-1); a[nod] = 0.75
    return np.nan_to_num(a, nan=0.75)
def cloudy(rgb):
    return (np.nan_to_num(rgb).mean(-1) > 0.25).mean() > 0.6

# ---------------------------------------------------------------- layout (mm)
FW = 190.0
ML, MR, MT, MB = 10.0, 2.0, 2.0, 2.0
PGAP, TG, BG = 3.6, 0.7, 2.4          # gap between panels, between tiles, between sensor blocks
PW = (FW - ML - MR - 3 * PGAP) / 4
T = (PW - 2 * TG) / 3                 # tile size
HEAD = 21.0
BH = 2 * T + TG
FH = MT + HEAD + 3 * BH + 2 * BG + MB
fig = plt.figure(figsize=(FW * MM, FH * MM))
def ax_mm(x, y, w, h): return fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH])   # y from top
def ftxt(x, y, s, **kw): return fig.text(x / FW, 1 - y / FH, s, **kw)

SENS = [('P', 'PlanetScope', '4.8 m'), ('S', 'Sentinel-2', '10–20 m'), ('L', 'Landsat 8/9', '30 m')]
for pi, k in enumerate('abcd'):
    s = SITES[k]; x0 = ML + pi * (PW + PGAP)
    post, dnbr, extra = s['info']
    ftxt(x0, MT + 0.3, f'({k})', fontsize=11, fontweight='bold', va='top')
    ftxt(x0 + 6.2, MT + 0.7, s['pid'], fontsize=9, fontweight='bold', va='top')
    import textwrap
    ftxt(x0, MT + 4.6, textwrap.fill(s['case'], 27), fontsize=8, style='italic', va='top',
         linespacing=1.15)
    ftxt(x0, MT + 12.2, f"{s['lat']:.4f} N, {s['lon']:.4f} E", fontsize=7.5, va='top')
    ftxt(x0, MT + 15.3, f"NBR$_{{post}}$ {post:.3f}   dNBR {dnbr:.3f}", fontsize=7.5, va='top')
    ftxt(x0, MT + 18.4, extra, fontsize=7.5, va='top')
    y = MT + HEAD
    for si, (sens, name, res) in enumerate(SENS):
        for i in range(6):
            r_, c_ = divmod(i, 3)
            ax = ax_mm(x0 + c_ * (T + TG), y + r_ * (T + TG), T, T)
            lab, rgb = chips[k, sens, i]
            if sens == 'L' and cloudy(rgb):
                ax.imshow(np.full((PX, PX, 3), 0.82), interpolation='nearest')
                ax.text(0.5, 0.5, 'no cloud-free\nLandsat image', transform=ax.transAxes, ha='center',
                        va='center', fontsize=6.5, color='#333333', linespacing=1.1)
            else:
                ax.imshow(show(sens, rgb), interpolation='nearest' if sens == 'L' else 'antialiased')
            ax.text(0.04, 0.96, lab, transform=ax.transAxes, fontsize=6.8, fontweight='bold', color='white',
                    va='top', ha='left', path_effects=HALO)
            # validation point: open ring (the pixel underneath stays visible)
            cpx = cpy = PX / 2 if not (sens == 'L' and cloudy(rgb)) else -100
            ax.plot(cpx, cpy, 'o', ms=4.2, mfc='none', mec='black', mew=1.6)
            ax.plot(cpx, cpy, 'o', ms=4.2, mfc='none', mec='#FFD21F', mew=0.9)
            ax.set_xlim(0, PX); ax.set_ylim(PX, 0); ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values(): sp.set_linewidth(0.4); sp.set_color('black')
            # one scale bar per sensor row (panel a, lower-left tile) and a north arrow (first tile)
            if k == 'a' and i == 3:
                m = PX / (2 * HALF)                         # pixels per metre
                bx, by, bh = 0.06 * PX, 0.88 * PX, 0.035 * PX
                for j, (a0, a1) in enumerate([(0, 250), (250, 500)]):
                    ax.add_patch(Rectangle((bx + a0 * m, by), (a1 - a0) * m, bh,
                                 facecolor='black' if j == 0 else 'white', edgecolor='black', lw=0.5))
                ax.text(bx, by - 0.02 * PX, '0', fontsize=6.5, color='white', ha='center', va='bottom',
                        path_effects=HALO)
                ax.text(bx + 500 * m, by - 0.02 * PX, '500 m', fontsize=6.5, color='white', ha='center',
                        va='bottom', path_effects=HALO)
            if k == 'a' and i == 0 and si == 0:
                ax.annotate('', xy=(0.88 * PX, 0.10 * PX), xytext=(0.88 * PX, 0.30 * PX),
                            arrowprops=dict(arrowstyle='-|>,head_width=0.25,head_length=0.4', lw=1.0, color='white',
                                            path_effects=[pe.withStroke(linewidth=2.2, foreground='black')]))
                ax.text(0.88 * PX, 0.36 * PX, 'N', fontsize=6.8, fontweight='bold', color='white', ha='center',
                        va='top', path_effects=HALO)
        if pi == 0:                                    # sensor label in the left margin
            ftxt(ML - 2.0, y + BH / 2, f'{name}\n({res})', fontsize=8.5, fontweight='bold', rotation=90,
                 ha='right', va='center', ma='center', linespacing=1.15)
        y += BH + BG

fig.savefig(f'{OUT}.pdf', dpi=600)
fig.savefig(f'{OUT}.png', dpi=600)
print(f'saved {OUT}: {FW:.0f} x {FH:.1f} mm; tile {T:.1f} mm = {2*HALF} m -> scale 1:{2*HALF/(T/1000):,.0f}')
