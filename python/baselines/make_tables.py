"""Build Table 3 (main text) and Table S20 (all variants) from baseline_scores.json.

Selection rule, fixed before looking at validation results: for each algorithm and record end,
report the variant with the highest calibration F1 (400 calibration points, 2025).
"""
import json, re

S = json.load(open('baseline_scores.json'))
f = lambda x: f'{x:.3f}'
ci = lambda c: f'[{c[0]:.3f}, {c[1]:.3f}]'


def family(name):
    m = re.match(r'(LandTrendr|CCDC).*record to (\d{4})', name)
    return (m.group(1), m.group(2)) if m else None


best = {}
for name, r in S.items():
    if name == 'rule':
        continue
    fam = family(name)
    if fam not in best or r['cal_F1'] > S[best[fam]]['cal_F1']:
        best[fam] = name

rows_main = []
rule = S['rule']
rows_main.append(['Anniversary rule (this study)', '–', f(rule['2025']['UA']), f(rule['2025']['PA']), f(rule['2025']['F1']), '–',
                  f(rule['2010']['UA']), f(rule['2010']['PA']), f(rule['2010']['F1']), '–'])
for fam in [('LandTrendr', '2025'), ('LandTrendr', '2026'), ('CCDC', '2025'), ('CCDC', '2026')]:
    if fam not in best:
        continue
    n = best[fam]; r = S[n]
    lab = n.replace(', record to', '; record to')
    row = [lab, f"{r['cal_F1']:.3f}"]
    for Y in ('2025', '2010'):
        e = r[Y]; d = e['diff_rule_minus_baseline']['F1']
        row += [f"{f(e['UA'])}", f"{f(e['PA'])}", f"{f(e['F1'])} {ci(e['ci']['F1'])}", f"{d[0]:+.3f} [{d[1][0]:+.3f}, {d[1][1]:+.3f}]"]
    rows_main.append(row)

O = json.load(open('oracle.json'))
rows_all = []
for name, r in S.items():
    if name == 'rule':
        continue
    row = [name, f"{r['threshold']:.4g}", f"{r['cal_F1']:.3f}"]
    for Y in ('2025', '2010'):
        e = r[Y]; d = e['diff_rule_minus_baseline']['F1']; mc = e['mcnemar_rule_correct_baseline_wrong, baseline_correct_rule_wrong, p']
        # discordant counts are unweighted and the sample is stratified on the rule's map, so no test is reported
        row += [f(e['UA']), f(e['PA']), f"{f(e['F1'])} {ci(e['ci']['F1'])}", f"{d[0]:+.3f} [{d[1][0]:+.3f}, {d[1][1]:+.3f}]",
                f(O[name][Y][0]), f"{mc[0]}/{mc[1]}", f"{e['mapped']:.0f}"]
    rows_all.append(row)
json.dump({'best': {f'{k[0]} {k[1]}': v for k, v in best.items()}, 'main': rows_main, 'all': rows_all},
          open('tables.json', 'w'), indent=1)
for r in rows_main:
    print(' | '.join(r))
