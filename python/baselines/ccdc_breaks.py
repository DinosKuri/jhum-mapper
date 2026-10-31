"""Share of 2025 validation points with an in-season CCDC break (NBR drop), by reference class (Section 5.1)."""
import json, statistics
import score_baselines as sb
for mo in (6, 3):
    d = json.load(open(f'ccdc_{mo}_2026.json'))
    for lab, name in (('1', 'jhum'), ('0', 'not jhum')):
        pts = [p for p in sb.v25 if sb.v25[p]['ref_jhum'] == lab]
        hit = sum(1 for p in pts if (sb.ccdc_loss(d.get(p), 2025) or 0) > 0)
        nb = statistics.mean(len([t for t in (d[p]['tBreak'] or []) if t]) for p in pts)
        print(f'minObs {mo}, {name}: {hit} of {len(pts)} points with an in-season break; {nb:.2f} breaks per pixel')
