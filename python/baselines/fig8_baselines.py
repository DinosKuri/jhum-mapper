"""Fig. 8: F1 of the anniversary rule against LandTrendr and CCDC configurations, 2025 (last year of the
record) and 2010 (mid-record). House style of Figs 2-7. Usage: python fig8_baselines.py <out_dir>"""
import os, sys, json, re
sys.argv, OUT = [sys.argv[0]], sys.argv[1]
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_baselines as sb

EE = os.path.dirname(os.path.abspath(__file__)) + '/'
S = json.load(open(EE + 'baseline_scores.json'))
O = json.load(open(EE + 'oracle.json'))
installed = {f.name for f in font_manager.fontManager.ttflist}
SERIF = 'Times New Roman' if 'Times New Roman' in installed else 'Liberation Serif'
INK = '#000000'
mpl.rcParams.update({
    'font.family': SERIF, 'font.size': 9.5, 'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 8, 'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
    'axes.edgecolor': INK, 'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'xtick.major.size': 3,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42, 'savefig.dpi': 600})
MM = 1 / 25.4
BLUE, RED, GREEN, ORANGE = '#2c5f8e', '#c1272d', '#1f7a3a', '#e08e2b'

CONFIGS = [('LandTrendr', 'standard', 'LandTrendr standard, 6 seg.'), ('LandTrendr', 'fast', 'LandTrendr fast, 6 seg.'),
           ('LandTrendr', 'fast8', 'LandTrendr fast, 8 seg.'), ('LandTrendr', 'fast10', 'LandTrendr fast, 10 seg.'),
           ('LandTrendr', 'fast12', 'LandTrendr fast, 12 seg.'), ('LandTrendr', 'fast16', 'LandTrendr fast, 16 seg.'),
           ('CCDC', '6', 'CCDC, 6 observations'), ('CCDC', '3', 'CCDC, 3 observations')]


def pick(alg, cfg, end):
    """Variant reported for a configuration: the extraction with the higher calibration F1."""
    if alg == 'CCDC':
        n = f'CCDC minObs {cfg}, record to {end}'
        return n if n in S else None
    c = [f'LandTrendr {cfg} {m}, record to {end}' for m in ('segment', 'annual')]
    c = [n for n in c if n in S]
    return max(c, key=lambda n: S[n]['cal_F1']) if c else None


rule_rows = {'2025': [(int(sb.v25[p]['stratum']), int(sb.v25[p]['map_landsat']), int(sb.v25[p]['ref_jhum'])) for p in sorted(sb.v25)],
             '2010': [(int(sb.v10[p]['stratum']), int(sb.v10[p]['map_jhum']), int(sb.v10[p]['ref_jhum'])) for p in sorted(sb.v10)]}
rule_ci = {Y: sb.boot_paired(rule_rows[Y], rule_rows[Y], sb.A25 if Y == '2025' else sb.A10)[0]['F1'] for Y in rule_rows}
rule_f1 = {Y: S['rule'][Y]['F1'] for Y in rule_rows}
rule_ci['2025'] = (0.542, 0.735)   # interval published for the rule (Section 4.1)
json.dump({'rule_F1_ci': rule_ci}, open(EE + 'rule_ci.json', 'w'))

fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 82 * MM), sharey=True)
fig.subplots_adjust(left=0.215, right=0.985, bottom=0.25, top=0.90, wspace=0.08)
rows = [c for c in CONFIGS if pick(c[0], c[1], 2026 if c[0] == 'LandTrendr' else '2026') or pick(c[0], c[1], 2025 if c[0] == 'LandTrendr' else '2025')]
for ax, Y, title in ((axes[0], '2025', 'a  2025: last season of the record'), (axes[1], '2010', 'b  2010: mid-record')):
    ax.axvspan(*rule_ci[Y], color='#f6e4e4', lw=0, zorder=0)
    ax.axvline(rule_f1[Y], color=RED, lw=1.1, zorder=1)
    ax.text(rule_f1[Y] - 0.012, len(rows) - 0.05, f'Rule {rule_f1[Y]:.2f}', color=RED, ha='right', va='center', fontsize=8.5)
    for i, (alg, cfg, lab) in enumerate(rows):
        y = len(rows) - 1 - i
        for end, dy, mfc, col in ((2025, 0.17, BLUE, BLUE), (2026, -0.17, 'white', BLUE)):
            n = pick(alg, cfg, end if alg == 'LandTrendr' else str(end))
            if not n: continue
            e = S[n][Y]
            ax.plot(e['ci']['F1'], [y + dy] * 2, color=col, lw=0.8, zorder=2)
            ax.plot(e['F1'], y + dy, 'o', ms=5.2, mfc=mfc, mec=col, mew=0.9, zorder=3)
            ax.plot(O[n][Y][0], y + dy, marker='|', ms=7, color=INK, mew=1.0, zorder=3)
    ax.set_xlim(0, 1); ax.set_xticks(np.arange(0, 1.01, 0.2))
    for x in np.arange(0, 1.01, 0.1): ax.axvline(x, color='#e3e3e3', lw=0.5, zorder=0)
    ax.set_xlabel('F1 (area-weighted, 95 % bootstrap interval)')
    ax.text(0, 1.04, title, transform=ax.transAxes, fontweight='bold', fontsize=10, ha='left', va='bottom')
    for sp in ('left',): ax.spines[sp].set_visible(False)
    ax.tick_params(axis='y', length=0)
axes[0].set_yticks(range(len(rows))); axes[0].set_yticklabels([r[2] for r in reversed(rows)])
axes[0].set_ylim(-0.6, len(rows) + 0.2)
fig.legend(handles=[
    plt.Line2D([], [], ls='none', marker='o', ms=5.2, mfc=BLUE, mec=BLUE, label='Record ends 2025'),
    plt.Line2D([], [], ls='none', marker='o', ms=5.2, mfc='white', mec=BLUE, mew=0.9, label='One further season (2026)'),
    plt.Line2D([], [], ls='none', marker='|', ms=7, color=INK, mew=1.0, label='Threshold tuned on the validation sample (upper bound)'),
    plt.Rectangle((0, 0), 1, 1, fc='#f6e4e4', ec=RED, lw=1.1, label='Anniversary rule, 95 % interval')],
    loc='lower center', ncol=2, bbox_to_anchor=(0.6, 0.0), handlelength=1.4, columnspacing=1.6)
os.makedirs(OUT, exist_ok=True)
fig.savefig(f'{OUT}/Fig8_baselines.png', dpi=600); fig.savefig(f'{OUT}/Fig8_baselines.pdf')
print('saved')
