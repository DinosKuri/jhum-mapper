import json, re
S = json.load(open('baseline_scores.json'))
T = json.load(open('tables.json'))
RC = json.load(open('rule_ci.json'))['rule_F1_ci']
RC['2025'] = [0.542, 0.735]   # interval published for the rule (Section 4.1)
def lab(n):
    m = re.match(r'LandTrendr (\w+) (segment|annual), record to (\d+)', n)
    if m:
        seg = {'standard': 'standard, 6', 'fast': 'fast, 6'}.get(m.group(1), 'fast, ' + m.group(1)[4:])
        return f'LandTrendr, {seg.replace(", ", " set, ")} segments' + ('' if m.group(2) == 'annual' else ' (segment)'), m.group(3)
    m = re.match(r'CCDC minObs (\d), record to (\d+)', n)
    return f'CCDC, {m.group(1)} observations', m.group(2)
L = ['**Table 3.** The anniversary rule against the LandTrendr and CCDC configurations with the highest calibration F1 for each record end. ΔF1 is rule minus method, from a paired stratified bootstrap; intervals are 95 %. LandTrendr loss is read as the fitted change from the previous year. All 28 configurations are in Table S20.', '',
     '| Method | Year | UA | PA | F1 [95 % CI] | ΔF1 [95 % CI] |', '|:--|:-:|--:|--:|:--|:--|']
r = S['rule']
for Y in ('2025', '2010'):
    L.append(f"| {'Anniversary rule' if Y=='2025' else ''} | {Y} | {r[Y]['UA']:.3f} | {r[Y]['PA']:.3f} | {r[Y]['F1']:.3f} [{RC[Y][0]:.3f}, {RC[Y][1]:.3f}] | – |")
order = [k for k in ('LandTrendr 2025', 'LandTrendr 2026', 'CCDC 2025', 'CCDC 2026') if k in T['best']]
for k in order:
    n = T['best'][k]; e = S[n]; name, end = lab(n)
    for Y in ('2025', '2010'):
        d = e[Y]['diff_rule_minus_baseline']['F1']; c = e[Y]['ci']['F1']
        first = Y == '2025'
        L.append(f"| {name + '; record to ' + end + f' (calibration F1 {e['cal_F1']:.3f})' if first else ''} | {Y} | {e[Y]['UA']:.3f} | {e[Y]['PA']:.3f} | {e[Y]['F1']:.3f} [{c[0]:.3f}, {c[1]:.3f}] | {d[0]:+.3f} [{d[1][0]:+.3f}, {d[1][1]:+.3f}] |")
open('table3.md', 'w').write('\n'.join(L).replace('-0.', '−0.').replace('+-', '−') + '\n')
print('\n'.join(L))
