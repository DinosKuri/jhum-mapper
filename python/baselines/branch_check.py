"""Accuracy of each branch of the rule alone, 2025 and 2010 (answers whether the narrow branch
explains the advantage over LandTrendr and CCDC). Needs rule_points_2025.json (rule_at_points.py)
and the data deposit. Usage: python branch_check.py"""
import json, csv
import score_baselines as sb

o = json.load(open('rule_points_2025.json'))
res = {}
rows = {'both': [], 'wide': [], 'narrow': []}
for p in sorted(sb.v25):
    r, q = sb.v25[p], o[p]
    s, y = int(r['stratum']), int(r['ref_jhum'])
    A = int(q['pw'] is not None and q['pw'] <= 0.26 and q['dw'] >= 0.20)
    B = int(q['pn'] is not None and q['pn'] <= 0.14 and q['dn'] >= 0.46)
    assert int(q['jh']) == int(r['map_landsat']), p
    rows['both'].append((s, int(r['map_landsat']), y)); rows['wide'].append((s, A, y)); rows['narrow'].append((s, B, y))
res['2025'] = rows
v10 = list(csv.DictReader(open(sb.DD + 'validation_points_160_2010.csv')))
res['2010'] = {'both': [(int(r['stratum']), int(r['map_jhum']), int(r['ref_jhum'])) for r in v10],
               'wide': [(int(r['stratum']), int(r['branch_wide']), int(r['ref_jhum'])) for r in v10],
               'narrow': [(int(r['stratum']), int(r['branch_narrow']), int(r['ref_jhum'])) for r in v10]}
out = {}
for Y, areas in (('2025', sb.A25), ('2010', sb.A10)):
    for b, rw in res[Y].items():
        e = sb.estimate(rw, areas)
        ci = sb.boot_paired(rw, rw, areas)[0]['F1']
        out[f'{Y} {b}'] = {k: e[k] for k in ('UA', 'PA', 'F1')} | {'F1_ci': ci}
        print(f"{Y} {b:6s} UA {e['UA']:.3f}  PA {e['PA']:.3f}  F1 {e['F1']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")
json.dump(out, open('branch_check.json', 'w'), indent=1, default=float)
