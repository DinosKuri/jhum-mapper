# Data dictionary

All coordinates are WGS 84 (EPSG:4326) pixel centres. Licence: CC BY 4.0 (see `LICENSE-data`).

## calibration_points_400_2025.csv
Stratified random sample from a preliminary (v3) map, interpreted on PlanetScope. Used only to
choose the compositing design and thresholds; never used for accuracy estimation.

| Column | Meaning |
|---|---|
| point_id | Identifier |
| lon, lat | Coordinates |
| stratum | Stratum of the preliminary map: 1 = mapped jhum (n = 100), 0 = background (n = 300) |
| ref_jhum | Reference label: 1 = cleared in the 2025 season (Nov 2024 – May 2025), 0 = not |
| ref_confidence | Interpreter confidence: h / m |
| ref_landcover | Reference land cover (controlled vocabulary) |
| notes | Interpretation tags (see `docs/interpretation_protocol.docx`) |
| pre_*, post_*, dnbr_* | NBR before, after and change under each compositing design: `v3` anniversary Mar–May 75/25; `sym25`, `sym50` anniversary Mar–May at the 25th / 50th percentile on both sides; `wide25`, `wide50` anniversary Nov–May (the final design uses `wide25`) |
| season_amplitude | Dry-season NBR amplitude |
| elevation, slope | SRTM terrain |

## validation_points_310_2025.csv and validation_points_160_2010.csv
Independent stratified random samples drawn from the final map, labelled blind (map class and
index values hidden), each point once.

| Column | Meaning |
|---|---|
| point_id | Identifier (V2025_/N2025_ = 2025, H2010_ = 2010) |
| year | 2010 file only |
| ref_jhum, ref_confidence, ref_landcover, notes | Reference interpretation |
| ref_source | 2010 file only: Google Earth Pro (140) or Landsat only (20) |
| stratum | 1 = mapped jhum; 2 = not mapped, within 300 m of mapped jhum; 3 = beyond 300 m |
| map_jhum | Map class at the point (1 = jhum) |
| branch_wide, branch_narrow | Which rule branch fired (1 = yes) |
| NBR_post_wide, dNBR_wide, NBR_post_narrow, dNBR_narrow | Rule inputs at the point |

Stratum areas (km²), the weights for area-weighted accuracy and error-adjusted area
(Olofsson et al., 2014; Stehman, 2014):

| Year | Stratum 1 | Stratum 2 | Stratum 3 |
|---|---|---|---|
| 2025 | 527.4 | 7,161.6 | 13,473.3 |
| 2010 | 375.1 | 5,118.3 | 15,679.6 |

Point counts: 2025 = 70 / 140 / 100; 2010 = 60 / 60 / 40.

## annual_jhum_area_1989_2025.csv

| Column | Meaning |
|---|---|
| year | Mapping year (burn season of that year) |
| jhum_km2 | Mapped area of the validated map, no minimum mapping unit (the series reported in the paper) |
| jhum_km2_mmu5 | The same map after a 5-pixel minimum mapping unit (optional filter) |
| valid_km2 | Area with a valid before/after image pair |
| coverage_pct | valid_km2 as a percentage of the state |
| jhum_km2_coverage_corrected | jhum_km2 scaled by coverage (assumes equal density in unobserved areas) |

Years 1989, 1991–92 and 1995–97 have 65–90 % coverage; use the corrected column with care.
The 2025 value (528.6 km²) is computed by summing pixel areas over the state; the 527.4 km² of
stratum 1 was computed on the sampling grid. The 0.2 % difference is a projection effect.

## experiments/

| File | Used for |
|---|---|
| sensor_consistency_L5_L7 / L7_L8 / L8_L9.csv | Cross-sensor NBR agreement over stable forest (2,000 pixels per pair) |
| sensor_transition_test.csv | Transition years rebuilt with and without each sensor |
| sentinel2_features_*.csv | Sentinel-2 features at the calibration and validation points |
| transfer_test_ne_india.csv | Exploratory application to seven north-east Indian states |
| control_recheck_2010.csv | Blind re-check of agreeing 2010 points (interpreter consistency) |
