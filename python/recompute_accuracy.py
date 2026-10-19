#!/usr/bin/env python3
"""
Recompute every 2025 accuracy figure from a pooled validation labels file.

Usage
    python recompute_accuracy.py <labels.csv> [stratum_areas.csv] [label_column]

    labels.csv        one row per point; needs the columns
                      stratum, map_landsat, map_s2 and a 0/1 label column
    stratum_areas.csv optional; columns stratum, area_km2
                      (defaults to the areas printed by the round-3 sampler)
    label_column      optional; defaults to ref_jhum

Prints: the stratum table, Landsat and Sentinel-2 accuracies with stratified
bootstrap intervals, the error-adjusted area, the four-product comparison
(Landsat / Sentinel-2 / union / intersection), where omission sits by stratum,
and the rule-of-three bound on any stratum that returned no jhum.

Estimators follow Olofsson et al. (2014); strata differ from map classes, so
Stehman (2014) applies. Intervals are percentile intervals from a stratified
bootstrap of 4000 resamples drawn within strata.
"""
import csv, io, sys, math, random, collections

DEFAULT_AREAS = {1: 425.6242859266101, 2: 102.52771525650927,
                 3: 208.47681360007508, 4: 8391.831029107088,
                 5: 12055.637613740282}
NAMES = {1: "A both maps", 2: "B Landsat only", 3: "C Sentinel-2 only",
         4: "D neither, within 300 m", 5: "E neither, beyond 300 m"}
NBOOT, SEED = 4000, 20261019


def load(path):
    with io.open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def estimate(rows, areas, lab, mapcol):
    """Area-weighted confusion cells, as proportions of the frame."""
    B = sum(areas.values())
    tp = fp = fn = 0.0
    counts = [0, 0, 0, 0]
    for s, area in areas.items():
        sub = [r for r in rows if int(r["stratum"]) == s]
        if not sub:
            continue
        W, n = area / B, len(sub)
        for r in sub:
            m, y = int(r[mapcol]), int(r[lab])
            if m and y:
                tp += W / n; counts[0] += 1
            elif m and not y:
                fp += W / n; counts[1] += 1
            elif y:
                fn += W / n; counts[2] += 1
            else:
                counts[3] += 1
    ua = tp / (tp + fp) if tp + fp else float("nan")
    pa = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * ua * pa / (ua + pa) if ua + pa else float("nan")
    return dict(OA=1 - (fp + fn), UA=ua, PA=pa, F1=f1,
                adj=B * (tp + fn), mapped=B * (tp + fp), counts=tuple(counts))


def boot(rows, areas, lab, mapcol, seed=SEED):
    random.seed(seed)
    B = sum(areas.values())
    cells = {s: [(int(r[mapcol]), int(r[lab])) for r in rows
                 if int(r["stratum"]) == s] for s in areas}
    keep = {s: v for s, v in cells.items() if v}
    out = []
    for _ in range(NBOOT):
        tp = fp = fn = 0.0
        for s, v in keep.items():
            W, n = areas[s] / B, len(v)
            for _ in range(n):
                m, y = v[random.randrange(n)]
                if m and y: tp += W / n
                elif m and not y: fp += W / n
                elif y: fn += W / n
        ua = tp / (tp + fp) if tp + fp else float("nan")
        pa = tp / (tp + fn) if tp + fn else float("nan")
        f1 = 2 * ua * pa / (ua + pa) if (ua + pa) and ua == ua and pa == pa else float("nan")
        out.append((1 - (fp + fn), ua, pa, f1, B * (tp + fn)))
    res = {}
    for i, k in enumerate(("OA", "UA", "PA", "F1", "adj")):
        v = sorted(x[i] for x in out if x[i] == x[i])
        res[k] = (v[int(0.025 * len(v))], v[int(0.975 * len(v))]) if v else (float("nan"),) * 2
    return res


def line(name, e, ci=None):
    s = "  %-22s OA %.3f  UA %.3f  PA %.3f  F1 %.3f   mapped %7.1f  adjusted %7.1f" % (
        name, e["OA"], e["UA"], e["PA"], e["F1"], e["mapped"], e["adj"])
    print(s)
    if ci:
        print("  %-22s    [%.3f, %.3f]  [%.3f, %.3f]  [%.3f, %.3f]  [%.3f, %.3f]"
              % ("", *ci["OA"], *ci["UA"], *ci["PA"], *ci["F1"]))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    rows = load(sys.argv[1])
    areas = dict(DEFAULT_AREAS)
    if len(sys.argv) > 2 and sys.argv[2] not in ("-", ""):
        areas = {int(r["stratum"]): float(r["area_km2"]) for r in load(sys.argv[2])}
    lab = sys.argv[3] if len(sys.argv) > 3 else "ref_jhum"
    if lab not in rows[0]:
        sys.exit("label column %r not in the file; columns are %s" % (lab, list(rows[0])))

    B = sum(areas.values())
    print("file            %s" % sys.argv[1])
    print("label column    %s" % lab)
    print("points          %d" % len(rows))
    print("frame           %.3f km2\n" % B)

    print("STRATUM TABLE")
    print("  %-26s %10s %8s %5s %8s %8s" % ("", "area km2", "weight", "n", "ref jhum", "km2/pt"))
    for s in sorted(areas):
        sub = [r for r in rows if int(r["stratum"]) == s]
        j = sum(1 for r in sub if r[lab] == "1")
        print("  %-26s %10.3f %8.5f %5d %8d %8.1f"
              % (NAMES.get(s, s), areas[s], areas[s] / B, len(sub), j,
                 areas[s] / len(sub) if sub else float("nan")))
    print("  %-26s %10.3f %8.5f %5d %8d"
          % ("total", B, 1.0, len(rows), sum(1 for r in rows if r[lab] == "1")))

    print("\nACCURACY  (95 % percentile intervals beneath each line)")
    for mapcol, nm in (("map_landsat", "Landsat 30 m"), ("map_s2", "Sentinel-2 20 m")):
        line(nm, estimate(rows, areas, lab, mapcol), boot(rows, areas, lab, mapcol))
    for mapcol, nm in (("map_landsat", "Landsat 30 m"), ("map_s2", "Sentinel-2 20 m")):
        e = estimate(rows, areas, lab, mapcol)
        print("  %-22s TP/FP/FN/TN counts %s" % (nm, e["counts"]))

    # union and intersection, derived from the two map columns
    for r in rows:
        r["_union"] = str(int(int(r["map_landsat"]) or int(r["map_s2"])))
        r["_inter"] = str(int(int(r["map_landsat"]) and int(r["map_s2"])))
    print("\nFOUR PRODUCTS")
    for mapcol, nm in (("map_landsat", "Landsat"), ("map_s2", "Sentinel-2"),
                       ("_union", "union, either"), ("_inter", "intersection, both")):
        line(nm, estimate(rows, areas, lab, mapcol))

    print("\nWHERE LANDSAT'S OMISSION SITS")
    tot = 0.0; parts = {}
    for s in sorted(areas):
        sub = [r for r in rows if int(r["stratum"]) == s]
        if not sub:
            continue
        W, n = areas[s] / B, len(sub)
        fn = sum(W / n for r in sub
                 if int(r["map_landsat"]) == 0 and r[lab] == "1")
        parts[s] = B * fn; tot += B * fn
    for s, v in parts.items():
        if v > 0:
            print("  %-26s %7.1f km2   %3.0f %% of all omission   (%.1f %% of the state)"
                  % (NAMES.get(s, s), v, 100 * v / tot, 100 * areas[s] / B))
    print("  %-26s %7.1f km2 omitted in total" % ("", tot))

    print("\nSTRATA THAT RETURNED NO JHUM")
    none = [s for s in sorted(areas)
            if [r for r in rows if int(r["stratum"]) == s]
            and not any(r[lab] == "1" for r in rows if int(r["stratum"]) == s)]
    if not none:
        print("  none - every stratum contains at least one reference jhum point")
    for s in none:
        n = len([r for r in rows if int(r["stratum"]) == s])
        print("  %-26s n=%d, zero positives; rule of three gives p <= %.4f, or <= %.0f km2 unseen"
              % (NAMES.get(s, s), n, 3 / n, 3 / n * areas[s]))

    # exact McNemar on the paired discordant points
    b = sum(1 for r in rows if (int(r["map_landsat"]) == int(r[lab])) and
            (int(r["map_s2"]) != int(r[lab])))
    c = sum(1 for r in rows if (int(r["map_s2"]) == int(r[lab])) and
            (int(r["map_landsat"]) != int(r[lab])))
    n = b + c
    if n:
        p = min(1.0, 2 * sum(math.comb(n, k) for k in range(min(b, c) + 1)) / 2 ** n)
        print("\nPAIRED TEST\n  %d discordant points; Landsat correct at %d, Sentinel-2 at %d; "
              "exact McNemar p = %.5f" % (n, b, c, p))


if __name__ == "__main__":
    main()
