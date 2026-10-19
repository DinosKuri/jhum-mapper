# Changelog

## v1.2.0 — 19 October 2026

The 2025 reference sample changed, and every 2025 figure changed with it.

### Changed
- **Round 2 of the 2025 validation sample was withdrawn in full and replaced.** Round 2 had been
  analyzed with round 1 for v1.1.0. It was then found that during its interpretation the
  PlanetScope basemap selector had made imagery from adjacent seasons available alongside the
  2025 season. The error could not be confined to identifiable points, so the whole 200-point
  round was set aside, after its estimates were known, and a third round of 200 was drawn in
  strata D and E with a new seed and a 100/100 allocation (`gee/20`; round 2 was 110/90) and
  interpreted under the same protocol. One of the withdrawn 200 was reference jhum. Pooling all
  three rounds (n = 720) gives 575 km² and Landsat producer's accuracy 0.617 (Table S17). The
  round is deposited unaltered in the dataset record under `withdrawn/`, and `gee/08` and
  `gee/20` carry a dated note.
- **All 520 reference labels were re-read once against the strictly point-based criterion of
  Text S1.** Eighteen labels changed, all in round 1, 17 of them from jhum to not jhum at points
  lying beside a clearing rather than within it; 9 of the 17 reversals fell in stratum B. The
  first-pass labels are deposited and reported as a sensitivity analysis (Tables S16, S17).
  Seventeen round-3 points carry notes beginning `reviewer_corrected`, recording a review step
  during interpretation (Section 2.3).
- 2025 results, Landsat: overall 0.981, user's 0.672, producer's 0.614, F1 0.642, error-adjusted
  area 578 km² (95 % CI 461–724). Sentinel-2: 0.987, 0.739, 0.811, 0.774. Exact McNemar on the 90
  discordant points, Landsat correct at 26 and Sentinel-2 at 64, p = 0.00008.
- Stratum areas for D and E are 8,391.8 and 12,055.6 km², from the round-3 sampler. The 3.6 km²
  shift in the D/E boundary against v1.1.0 changes no estimate to three decimal places.
- `figure_data/figure7a_accuracy.csv`, `figure7b_disagreement.csv` and the five-stratum row of
  `figure_data/fig10_two_year_validation.csv` rebuilt; Figures 5 and 7 regenerated.
- Supplementary workbook now runs to Table S17. Tables S2c, S3, S5, S9, S12 and S15 and the
  notes to Tables S1, S2, S13 and S14 rebuilt on the 520-point file with re-read labels; S17
  (sensitivity to sample composition and to R191) added.
- `figure_data/fig1_validation_points.csv` rebuilt from the 520-point file (it still held the
  withdrawn round), and `python/fig1/plot_fig1.py` now reads the stratum counts from the data.
- Headers of `gee/00`, `gee/06` and `gee/19` updated to the published values.
- `docs/reference_data_guide.docx` removed; it described the withdrawn round. The dataset
  README is the data dictionary.
- `python/fig9/plot_fig9.py` legend now reads 578 km² (461–724); it carried the superseded
  error-adjusted area.
- The README no longer attributes the suppression of commission error to calendar alignment
  alone. The operative constraint is applying the same percentile to both sides of the pair;
  Section 4.3 of the paper and `gee/00` already said so.

### Added
- `gee/20_dualmap_topup_round3_DE.js` — the round-3 sampler for strata D and E.
- `python/recompute_accuracy.py` — regenerates every 2025 figure from a labels file and a stratum
  area file, so the published numbers can be reproduced in one command.

## v1.1.0 — 16 October 2026

The 2025 validation design changed, and everything downstream of it changed with it.

### Added
- `gee/07_dualmap_sample_2025.js` — the published validation design: five strata defined
  jointly on the Landsat 30 m and Sentinel-2 20 m maps (A both, B Landsat only,
  C Sentinel-2 only, D neither but within 300 m, E remainder), with an agreement raster
  export for the figures.
- `gee/08_dualmap_topup_strata_DE.js` — the 200-point enlargement of strata D and E that
  brings the 2025 sample to 520 points.
- `gee/06_sentinel2_v6rule.js` — the Sentinel-2 implementation of the rule.
- `gee/14_historical_validation_2010_round2.js` and
  `gee/superseded/S3_validation_extension_far.js`, which were missing from v1.0.0.
- `python/` — the Matplotlib scripts that draw every figure, and
  `python/gee_exports/` — the Earth Engine exports they consume.
- `docs/reference_data_guide.docx` — data dictionary for the 520-point reference file.

### Changed
- **`gee/06` now carries the refit Sentinel-2 thresholds** (wide branch NBR_post ≤ 0.27,
  dNBR ≥ 0.25). v1.0.0 shipped a Sentinel-2 script running the Landsat wide branch
  (0.26 / 0.20), which reproduces none of the published Sentinel-2 numbers.
- GEE scripts renumbered so the published design comes before the superseded one. The four
  superseded scripts moved to `gee/superseded/` and carry a banner saying so.
- `README.md` rebuilt on the 520-point results. v1.0.0 published the superseded figures
  (310 points, user's 0.886, producer's 0.901).
- `figure_data/fig1_validation_points.csv` is now the 520 points with strata A–E.
- `figure_data/fig10_two_year_validation.csv` now holds three rows: 2010 three-stratum,
  2025 three-stratum and 2025 five-stratum, so Figure 5 shows the same map under two designs.
- `python/fig2_4-8/plot_line_art.py` now draws Figure 7 as the Landsat–Sentinel-2 comparison
  from `figure7a_accuracy.csv` and `figure7b_disagreement.csv`.
- `python/fig9/plot_fig9.py` legend now reads 528.2 km² mapped, 691 km² (566–854)
  error-adjusted.
- `python/fig1/plot_fig1.py` draws five strata in panel a and three in panel b.
- The strata preview layer in `gee/07` and `gee/08` is off by default; at state-wide zoom it
  overruns the interactive memory limit without affecting the exports.

### Removed
- `data/` — the reference tables now live in the separate CC BY 4.0 dataset record, because
  the licences differ. `figure_data/` keeps only what the figure scripts read.
- `figure_data/fig6_sensors.csv` and `fig7_accuracy_strata.csv` — inputs to the superseded
  Sentinel-2 chart and the three-stratum design. They remain in the dataset record.
- `docs/` manuscript drafts and `figures/` PNGs, which were superseded snapshots.

### Fixed
- **Table 4 gave the v3 design a user's accuracy of 0.886, the same value as v6.** The v3 map's
  own check points settle it: of the 100 mapped-jhum points drawn from that map, 50 are reference
  jhum, so its user's accuracy is 0.50. That is the over-calling named in the header of gee/04 as
  the reason the compositing experiment was run. The check-point key is now deposited as
  experiments/v3_preliminary_map_checkpoints_400_2025.csv.
- **The v3 row of the compositing experiment was wrong.** Table 4, Table S1 and Figure 6 gave the
  v3 non-jhum median dNBR as 0.004 and its separation ratio as 157. Recomputed from the deposited
  400 calibration points, `dnbr_v3` gives 0.567 / 0.122 and a ratio of 4.6. The corrected figure
  reverses the paper's conclusion in Section 4.3: pairing identical calendar months does not by
  itself suppress background change, and applying the same percentile to both sides is what does.
  That now agrees with Section 3.2, with the highlights, and with the header of `gee/00`.
  `figure_data/fig5_compositing.csv` carries the corrected values and gains the `sym25` design,
  which isolates percentile symmetry from window length; `plot_line_art.py` draws five designs at
  170 mm.

- `CITATION.cff` carried a placeholder ORCID and repository URL.
- The header of `gee/00` advertised the superseded 2025 accuracies as the validated ones, and
  `gee/09` still said "validated on 310".
- The header of `gee/06` contradicted its own code: it described the Landsat thresholds as the
  defaults and concluded from the 240-point test that Sentinel-2 did not outperform Landsat.
- `.gitignore` carried `gee_exports/`, which would have excluded `python/gee_exports/` from the
  repository. It now ignores `*.qgz` and `*.qgs` instead, because a saved QGIS project embeds
  any Planet API key in plain text.
- `validation_points_160_2010.csv` in the dataset record had no `stratum` or `map_jhum` column,
  so the published 2010 numbers could not be reproduced from it. Both are now merged in, and
  the three-stratum areas for the superseded 2025 design are deposited alongside it.
- Seven CSVs in `figure_data/` were read by no script and two of them carried superseded
  accuracies. Removed, and `figure_data/README.md` now maps each remaining file to its figure.
- The graphical abstract drew one card per row of `fig10_two_year_validation.csv`, which since
  v1.1.0 holds three rows; it now skips the superseded 2025 round.
- `python/fig3/plot_fig3.py` drew Figure 3 at 240 mm, wider than Elsevier's 190 mm maximum, so
  the journal would have scaled it down and taken its smallest labels under 5 pt. It now draws
  at 190 mm; the site description wraps, the index line is split in two, and no label is below
  6.5 pt. Verified by rendering.
- `python/fig1/prep_relief_30m.py` and `prep_boundaries_relief.py` never wrote `india.gpkg`,
  which `plot_fig1.py` reads for the locator inset, so Figure 1 could not be drawn from a clean
  checkout.
- `gee/18` and `python/gee_exports/fig1_fig9_export_layers.js` export an SRTM hillshade for a
  QGIS route that was not used for any published figure, and print a 2025 area of 527.4 km²
  against the paper's 528.2 km². Both now say so in their headers.

## v1.0.0 — 1 October 2026

First release: Earth Engine scripts for the rule, the compositing experiment, the
single-map validation design for 2025 and 2010, the sensor tests and the regional transfer
test, with the calibration and validation point tables.
