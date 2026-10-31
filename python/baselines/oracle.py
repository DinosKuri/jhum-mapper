"""Optimistic bound for the baselines: threshold chosen on the validation sample itself (area-weighted F1).
Not a fair estimate; it shows the best F1 each variant could reach with any single magnitude threshold."""
import os, json, sys
sys.argv = [sys.argv[0]]
import score_baselines as sb

def oracle(fn, data, ref, areas, Y):
    sc = {p: fn(data.get(p), Y) for p in ref}
    grid = sorted({0.0} | {v for v in sc.values() if v and v > 0})
    best = (-1, None)
    for t in grid:
        rows = [(int(ref[p]['stratum']), int(sc[p] is not None and sc[p] > t), int(ref[p]['ref_jhum'])) for p in sorted(ref)]
        f = sb.estimate(rows, areas)['F1']
        if f == f and f > best[0]: best = (f, t)
    return best

out = {}
for ps in ('standard', 'fast', 'fast8', 'fast10', 'fast12', 'fast16'):
    for end in (2025, 2026):
        f = f'lt_{ps}_{end}.json'
        if not os.path.exists(f): continue
        d = json.load(open(f))
        for mode in ('segment', 'annual'):
            fn = (lambda m: (lambda a, Y: sb.lt_loss(a, Y, m)))(mode)
            out[f'LandTrendr {ps} {mode}, record to {end}'] = {Y: oracle(fn, d, ref, A, int(Y)) for Y, ref, A in (('2025', sb.v25, sb.A25), ('2010', sb.v10, sb.A10))}
for mo in (6, 3):
    for end in ('2025', '2026'):
        f = f'ccdc_{mo}_{end}.json'
        if not os.path.exists(f): continue
        d = json.load(open(f))
        out[f'CCDC minObs {mo}, record to {end}'] = {Y: oracle(sb.ccdc_loss, d, ref, A, int(Y)) for Y, ref, A in (('2025', sb.v25, sb.A25), ('2010', sb.v10, sb.A10))}
json.dump(out, open('oracle.json', 'w'), indent=1)
for k, v in out.items(): print(f"{k:45s} 2025 {v['2025'][0]:.3f}  2010 {v['2010'][0]:.3f}")
