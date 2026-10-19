"""Figures 2, 4, 5, 6, 7, 8 in the house style of Figs 1, 3 and 9:
Times New Roman (Liberation Serif fallback), all text black, 9-10 pt at print size, 600 dpi PNG + vector PDF.
Usage: python plot_line_art.py <figure_data_dir> <out_dir>
Manuscript numbering: Fig 4 = feature space, Fig 5 = accuracy 2010 vs 2025, Fig 8 = annual series.
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from scipy.stats import theilslopes

D, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
installed = {f.name for f in font_manager.fontManager.ttflist}
SERIF = 'Times New Roman' if 'Times New Roman' in installed else 'Liberation Serif'
INK = '#000000'
mpl.rcParams.update({
    'font.family': SERIF, 'mathtext.fontset': 'custom', 'mathtext.rm': SERIF, 'mathtext.it': SERIF + ':italic',
    'font.size': 9.5, 'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'axes.edgecolor': INK,
    'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.major.size': 3, 'ytick.major.size': 3, 'xtick.direction': 'out', 'ytick.direction': 'out',
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
    'pdf.fonttype': 42, 'savefig.dpi': 600})
MM = 1 / 25.4
BLUE, RED, GREEN, ORANGE = '#2c5f8e', '#c1272d', '#1f7a3a', '#e08e2b'
GRID = dict(color='#e3e3e3', lw=0.5, zorder=0)


PAD_IN = 0.04
def save(fig, name, width_mm):
    # adjust the figure width until the trimmed (tight) output is exactly width_mm wide
    target = width_mm * MM
    for _ in range(6):
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        out_w = bb.width + 2 * PAD_IN
        if abs(out_w - target) < 0.002: break
        w, h = fig.get_size_inches()
        fig.set_size_inches(w + (target - out_w), h)
    fig.savefig(f'{OUT}/{name}.pdf', bbox_inches='tight', pad_inches=PAD_IN)
    fig.savefig(f'{OUT}/{name}_600dpi.png', dpi=600, bbox_inches='tight', pad_inches=PAD_IN)
    plt.close(fig)
    print('saved', name)


def letter(ax, s, x=-0.10, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')


# ------------------------------------------------------------------ Fig 2: composite windows (even month axis)
# month index 0 = 1 Nov of year Y-2; a bar [a, b] covers the start of month a to the end of month b-1
fig, ax = plt.subplots(figsize=(140 * MM, 60 * MM))
rows = {'Wide pair': 2.0, 'Narrow pair': 1.0, 'Jhum calendar': 0.0}
H = 0.46
def bar(a, b, y, col, txt, tcol='white', bold=True):
    ax.add_patch(Rectangle((a, y - H / 2), b - a, H, fc=col, ec='none', zorder=3))
    ax.text((a + b) / 2, y, txt, ha='center', va='center', color=tcol, fontsize=8.5,
            fontweight='bold' if bold else 'normal', zorder=4)
bar(0, 7, 2, BLUE, 'Before (Nov–May)')
bar(12, 19, 2, RED, 'After (Nov–May)')
bar(4, 7, 1, '#8fb3d3', 'Before', INK)
bar(16, 19, 1, '#f2a07b', 'After', INK)
bar(12, 15, 0, '#bdbdbd', 'Slash', INK)
bar(15, 18, 0, '#5a5a5a', 'Burn')
bar(19, 23, 0, GREEN, 'Crop')
# 12-month offset arrows
for (a, b, y) in [(0, 12, 2.42), (4, 16, 1.42)]:
    ax.annotate('', xy=(b, y), xytext=(a, y), arrowprops=dict(arrowstyle='<->', lw=0.7, color=INK,
                                                             shrinkA=0, shrinkB=0), zorder=2)
ax.text(6, 2.47, '12 months', ha='center', va='bottom', fontsize=9, style='italic')
ax.text(10, 1.47, '12 months', ha='center', va='bottom', fontsize=9, style='italic')
ax.axvspan(12, 19, color='#f6e7e7', zorder=0, lw=0)
ax.text(15.5, -0.55, 'clearing season of mapping year Y', ha='center', va='center', fontsize=9, style='italic')
ax.set_yticks(list(rows.values()), list(rows.keys()))
ax.tick_params(axis='y', length=0, labelsize=10)
months = ['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O'] * 2
ax.set_xticks(np.arange(0, 24) + 0.5, months, fontsize=8)
ax.tick_params(axis='x', length=0, pad=2)
for m in range(0, 25):
    ax.plot([m, m], [-0.82, -0.74], color=INK, lw=0.6, clip_on=False)
for x, s in [(0, 'Y−2'), (2, 'Y−1'), (14, 'Y')]:
    ax.text(x + 0.05, -1.18, s, fontsize=9, ha='left', va='top', transform=ax.transData, fontweight='bold')
ax.set_xlim(-0.2, 24); ax.set_ylim(-0.78, 2.75)
ax.spines['left'].set_visible(False); ax.spines['bottom'].set_position(('data', -0.78))
ax.set_xlabel('Month', labelpad=16)
save(fig, 'Fig2_windows', 140)

# ------------------------------------------------------------------ Fig 4: calibration feature space
fs = pd.read_csv(f'{D}/fig3_featurespace_points.csv')
c = fs[fs.set == 'calibration']
fig, ax = plt.subplots(figsize=(90 * MM, 80 * MM))
ax.add_patch(Rectangle((-0.65, 0.20), 0.65 + 0.26, 1.2, fc='#f6e4e4', ec='none', zorder=0))
ax.plot([-0.65, 0.26, 0.26], [0.20, 0.20, 1.4], ls='--', lw=0.9, color=INK, zorder=2)
n0, n1 = c[c.ref_jhum == 0], c[c.ref_jhum == 1]
ax.scatter(n0.post_wide25, n0.dnbr_wide25, s=11, marker='o', fc=GREEN, ec='none', alpha=0.55, zorder=3,
           label=f'not jhum (n = {len(n0)})')
ax.scatter(n1.post_wide25, n1.dnbr_wide25, s=20, marker='^', fc=RED, ec='white', lw=0.3, zorder=4,
           label=f'jhum (n = {len(n1)})')
ax.text(-0.62, -0.06, 'Shaded: branch A\ndecision region', fontsize=9, fontweight='bold', va='top', linespacing=1.1)
ax.text(-0.62, -0.20, 'NBR$_{\\mathrm{post}}$ ≤ 0.26\ndNBR ≥ 0.20', fontsize=9, va='top', linespacing=1.3)
ax.set_xlim(-0.65, 0.8); ax.set_ylim(-0.4, 1.02)
ax.set_xlabel(r'NBR$_{\mathrm{post}}$ (wide window, 25th percentile)')
ax.set_ylabel('dNBR (wide window)')
ax.grid(**GRID); ax.set_axisbelow(True)
ax.legend(loc='upper right', handletextpad=0.3, borderaxespad=0.3)
save(fig, 'Fig4_featurespace', 90)

# ------------------------------------------------------------------ Fig 5: accuracy by sampling design
acc = pd.read_csv(f'{D}/fig10_two_year_validation.csv')
fig, ax = plt.subplots(figsize=(140 * MM, 76 * MM))
metrics = [('UA', "User's\naccuracy"), ('PA', "Producer's\naccuracy"), ('F1', 'F1')]
cols = [BLUE, ORANGE, RED]
n = len(acc); w = 0.78 / n
for j, (r, col) in enumerate(zip(acc.itertuples(), cols)):
    xs = np.arange(3) + (j - (n - 1) / 2) * w
    vals = [r.UA, r.PA, r.F1]
    lab = f'{r.year}  {r.sensors}\n{r.design}  (n = {r.n_points})'
    ax.bar(xs, vals, w * 0.88, color=col, zorder=3, label=lab)
    for k, (m, _) in enumerate(metrics):
        top = vals[k]
        if m != 'F1':
            lo, hi = getattr(r, f'{m}_lo'), getattr(r, f'{m}_hi')
            ax.plot([xs[k]] * 2, [lo, hi], color=INK, lw=0.8, zorder=4)
            for yy in (lo, hi):
                ax.plot([xs[k] - 0.035, xs[k] + 0.035], [yy, yy], color=INK, lw=0.8, zorder=4)
            top = hi
        ax.text(xs[k], top + 0.02, f'{vals[k]:.2f}', ha='center', va='bottom', fontsize=8.5)
ax.set_xticks(range(3), [m[1] for m in metrics], fontsize=10)
ax.set_ylim(0, 1.26); ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.set_ylabel('Accuracy')
ax.yaxis.grid(True, **GRID); ax.set_axisbelow(True)
ax.legend(loc='lower center', ncol=3, bbox_to_anchor=(0.5, 1.0), handlelength=1.1,
          columnspacing=1.4, fontsize=8)
save(fig, 'Fig5_accuracy_eras', 140)

# ------------------------------------------------------------------ Fig 6: compositing designs
cp = pd.read_csv(f'{D}/fig5_compositing.csv').sort_values('order')
labs = ['v1\nmedian\nadjacent', 'v2\nasym.\nOct–Jan', 'v3\nanniv.\nMar–May\n75/25',
        'sym25\nanniv.\nMar–May\n25/25', 'v6 (final)\nanniv.\nNov–May\n25/25']
fig, (a1, a2) = plt.subplots(1, 2, figsize=(170 * MM, 72 * MM), gridspec_kw=dict(width_ratios=[1.65, 1], wspace=0.30))
x = np.arange(len(cp)); w = 0.36
a1.bar(x - w / 2, cp.jhum_dnbr_median, w * 0.94, color=RED, label='jhum', zorder=3)
a1.bar(x + w / 2, cp.background_dnbr_median, w * 0.94, color=GREEN, label='non-jhum', zorder=3)
for i, (j, b) in enumerate(zip(cp.jhum_dnbr_median, cp.background_dnbr_median)):
    a1.text(i - w / 2, j + 0.012, f'{j:.3f}', ha='center', va='bottom', fontsize=8.5)
    a1.text(i + w / 2, b + 0.012, f'{b:.3f}', ha='center', va='bottom', fontsize=8.5)
a1.set_xticks(x, labs, fontsize=8.5); a1.set_ylim(0, 0.80); a1.set_ylabel('Median dNBR')
a1.yaxis.grid(True, **GRID); a1.set_axisbelow(True)
a1.legend(loc='upper left', ncol=2, handlelength=1.2)
letter(a1, 'a', -0.13)
yy = np.arange(len(cp))[::-1]
a2.barh(yy, cp.separation_ratio, 0.55, color=BLUE, zorder=3)
for y_, v in zip(yy, cp.separation_ratio):
    a2.text(v * 1.15, y_, f'{v:g}', va='center', fontsize=9)
a2.set_xscale('log'); a2.set_xlim(1, 100)
a2.set_yticks(yy, ['v1', 'v2', 'v3', 'sym25', 'v6'], fontsize=10)
a2.set_xlabel('Separation ratio (log scale)')
a2.xaxis.grid(True, **GRID); a2.set_axisbelow(True)
letter(a2, 'b', -0.16)
save(fig, 'Fig6_compositing', 170)

# ------------------------------------------------------------------ Fig 7: Landsat against Sentinel-2
a7 = pd.read_csv(f'{D}/figure7a_accuracy.csv')
b7 = pd.read_csv(f'{D}/figure7b_disagreement.csv')
CONN, BAND, PALE = '#9a9a9a', '#eef2f7', '#cfcfcf'
fig = plt.figure(figsize=(140 * MM, 86 * MM))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.18], left=0.205, right=0.885,
                      bottom=0.100, top=0.915, hspace=0.62)

ax = fig.add_subplot(gs[0, 0])
rows7 = [('Landsat 30 m', True), ('Sentinel-2 20 m', False)]
X0, X1 = 0.35, 1.00
brk = np.arange(X0, X1 + 1e-9, 0.05)
ax.set_xlim(X0, X1); ax.set_ylim(-1.15, 1.55)
for x in brk: ax.axvline(x, color='#e3e3e3', lw=0.5, zorder=0)
for i7, (key, hl) in enumerate(rows7):
    y = len(rows7) - 1 - i7
    if hl: ax.axhspan(y - 0.46, y + 0.46, color=BAND, zorder=0.5)
    g = a7[a7.sensor == key].set_index('metric')
    ua, pa = g.loc['users'], g.loc['producers']
    for rec, col, dy in ((pa, ORANGE, -0.30), (ua, BLUE, 0.30)):
        ax.plot([rec.lo, rec.hi], [y + dy] * 2, color=col, lw=0.8, zorder=2)
        for xc in (rec.lo, rec.hi):
            ax.plot([xc] * 2, [y + dy - 0.075, y + dy + 0.075], color=col, lw=0.8, zorder=2)
    ax.plot([pa.value, ua.value], [y, y], color=CONN, lw=2.6, solid_capstyle='round', zorder=3)
    ax.plot(pa.value, y, 'o', ms=6.2, color=ORANGE, mec='white', mew=0.6, zorder=4)
    ax.plot(ua.value, y, 'o', ms=6.2, color=BLUE, mec='white', mew=0.6, zorder=4)
    ax.text(1.055, y, '%.2f' % g.loc['f1'].value, transform=ax.get_yaxis_transform(),
            ha='center', va='center', fontsize=8.5, clip_on=False)
ax.text(1.055, 1.045, 'F1', transform=ax.transAxes, ha='center', va='bottom',
        fontsize=8.5, fontweight='bold', clip_on=False)
ax.text(0.0, 1.045, 'Exact McNemar $p$ = 0.00008 (90 discordant of 520 paired points)',
        transform=ax.transAxes, ha='left', va='bottom', fontsize=7.5, color='#4a4a4a')
ax.set_yticks(range(len(rows7))); ax.set_yticklabels([r[0] for r in reversed(rows7)])
ax.set_xticks(brk); ax.set_xticklabels(['%.2f' % v for v in brk])
ax.set_xlabel('Accuracy (520-point dual-map reference sample)', labelpad=3)
for sp in ('top', 'right', 'left'): ax.spines[sp].set_visible(False)
ax.tick_params(axis='y', length=0, pad=4); ax.tick_params(axis='x', length=2.5, pad=2)
ax.legend(handles=[
    plt.Line2D([], [], ls='none', marker='o', ms=6.2, mfc=ORANGE, mec='white', mew=0.6,
               label="Producer's accuracy"),
    plt.Line2D([], [], ls='none', marker='o', ms=6.2, mfc=BLUE, mec='white', mew=0.6,
               label="User's accuracy"),
    plt.Line2D([], [], color='#6f6f6f', lw=0.8, label='95 % bootstrap interval')],
    loc='lower center', bbox_to_anchor=(0.5, -0.045), ncol=3, handlelength=1.1,
    handletextpad=0.45, columnspacing=1.6, fontsize=8)
ax.text(-0.255, 1.12, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold',
        ha='left', va='bottom')

ax2 = fig.add_subplot(gs[1, 0])
lab7 = {'Both maps': 'Both maps', 'Landsat only': 'Landsat 30 m only',
        'Sentinel-2 only': 'Sentinel-2 20 m only'}
ax2.set_xlim(0, 500); ax2.set_ylim(-1.15, 2.55)
for x in range(0, 501, 50): ax2.axvline(x, color='#e3e3e3', lw=0.5, zorder=0)
for i7, cls in enumerate(['Both maps', 'Landsat only', 'Sentinel-2 only']):
    y = 2 - i7
    r = b7[b7['class'] == cls].iloc[0]
    ax2.barh(y, r.real_jhum_km2, height=0.44, color=BLUE, lw=0, zorder=2)
    ax2.barh(y, r.not_jhum_km2, left=r.real_jhum_km2, height=0.44, color=PALE, lw=0, zorder=2)
    ax2.text(1.055, y, '%.0f' % (100 * r.frac_real), transform=ax2.get_yaxis_transform(),
             ha='center', va='center', fontsize=8.5, clip_on=False)
ax2.text(1.055, 1.045, '% reference\njhum', transform=ax2.transAxes, ha='center', va='bottom',
         fontsize=8.5, fontweight='bold', clip_on=False, linespacing=1.2)
ax2.set_yticks(range(3))
ax2.set_yticklabels([lab7[c] for c in reversed(['Both maps', 'Landsat only', 'Sentinel-2 only'])])
ax2.set_xticks(range(0, 501, 50))
ax2.set_xlabel('Mapped area in 2025 (km$^2$)', labelpad=3)
for sp in ('top', 'right', 'left'): ax2.spines[sp].set_visible(False)
ax2.tick_params(axis='y', length=0, pad=4); ax2.tick_params(axis='x', length=2.5, pad=2)
ax2.legend(handles=[Rectangle((0, 0), 1, 1, fc=BLUE, ec='none', label='Reference jhum'),
                    Rectangle((0, 0), 1, 1, fc=PALE, ec='none', label='Commission')],
           loc='lower center', bbox_to_anchor=(0.5, -0.045), ncol=2, handlelength=1.1,
           handleheight=0.9, handletextpad=0.45, columnspacing=1.8, fontsize=8)
ax2.text(-0.255, 1.12, 'b', transform=ax2.transAxes, fontsize=10, fontweight='bold',
         ha='left', va='bottom')
save(fig, 'Fig7_sensors', 140)

# ------------------------------------------------------------------ Fig 8: annual series
a = pd.read_csv(f'{D}/fig4_annual_series.csv')
fig, (ax, bx) = plt.subplots(2, 1, figsize=(190 * MM, 112 * MM), sharex=True,
                             gridspec_kw=dict(height_ratios=[3.1, 1], hspace=0.08))
ax.axvspan(1988.5, 1997.5, color='#efefef', zorder=0, lw=0)
ax.text(1993, 1570, 'Incomplete coverage', ha='center', va='center', fontsize=9, style='italic')
ax.fill_between(a.year, a.primary_km2, a.coverage_corrected_km2, color='#f5c9a6', lw=0, zorder=1,
                label='Coverage-corrected area')
ax.plot(a.year, a.primary_km2, color=RED, lw=1.3, marker='o', ms=3.2, zorder=3, label='Mapped area')
s98 = a[a.year >= 1998]
sl, ic, lo, hi = theilslopes(s98.primary_km2, s98.year)
xx = np.array([1998, 2025])
ax.plot(xx, ic + sl * xx, ls='--', color=INK, lw=1.0, zorder=2,
        label=f"Sen's slope 1998–2025, {sl:.1f} km² yr$^{{-1}}$".replace('-', '−'))
from matplotlib.lines import Line2D
means = []
for (y0, y1) in [(1998, 2009), (2010, 2019), (2020, 2025)]:
    m = a[(a.year >= y0) & (a.year <= y1)].primary_km2.mean(); means.append((y0, y1, m))
    ax.plot([y0 - 0.3, y1 + 0.3], [m, m], color=BLUE, lw=2.4, zorder=2, solid_capstyle='butt')
ax.text(2015.6, 1170, 'Period means\n' + '\n'.join(f'{y0}–{y1}:  {m:,.0f} km²' for y0, y1, m in means),
        ha='left', va='top', fontsize=9, linespacing=1.35)
ax.set_ylabel('Jhum area (km²)')
ax.set_ylim(0, 1650); ax.yaxis.grid(True, **GRID); ax.set_axisbelow(True)
h, l = ax.get_legend_handles_labels()
order = [1, 0, 2]
hh = [h[i] for i in order] + [Line2D([], [], color=BLUE, lw=2.4)]
ll = [l[i] for i in order] + ['Period mean']
ax.legend(hh, ll, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=4, handlelength=1.6, columnspacing=1.2)
letter(ax, 'a', -0.115, 1.06)
bx.bar(a.year, a.coverage_pct, 0.72, color='#9fb4cc', zorder=2)
bx.axhline(99, color=BLUE, lw=0.8, ls=':', zorder=3)
bx.text(2025.9, 99, '99 %', va='center', ha='left', fontsize=8.5)
bx.set_ylim(50, 102); bx.set_yticks([60, 80, 100])
bx.set_ylabel('Coverage (%)'); bx.yaxis.grid(True, **GRID); bx.set_axisbelow(True)
bx.set_xlim(1988.3, 2025.7); bx.set_xticks(range(1990, 2026, 5)); bx.set_xlabel('Year')
letter(bx, 'b', -0.115, 0.95)
save(fig, 'Fig8_annual_series', 190)
print('font', SERIF, 'sen', round(sl, 2), round(lo, 1), round(hi, 1))
