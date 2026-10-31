# LandTrendr and CCDC baselines

Scripts for the comparison in Section 3.3, Section 4.4, Table 3, Fig. 8 and Table S20 of the paper.

Both algorithms run in Google Earth Engine at the 1,080 reference points (400 calibration, 520 validation for 2025, 160 validation for 2010) on the same harmonized Landsat collection as the anniversary rule. Preprocessing in `ee_common.py` is a line-for-line port of `gee/00_final_map_mizoram.js`. At the 520 points of 2025 the port reproduces the published map class of every point.

## Setup

```
pip install earthengine-api numpy matplotlib scipy
export EE_PROJECT=<your Earth Engine cloud project>        # or EE_SERVICE_ACCOUNT_KEY=<key file>
export JHUM_DATA=<folder of the data deposit, jhum-mapper-data-v1.2.0>
```

## Order

| Step | Script | Output |
|---|---|---|
| 1 | `python run_landtrendr.py standard fast fast8 fast10 fast12 fast16` | `lt_<set>_<end>.json` |
| 2 | `python run_ccdc.py` | `ccdc_<minObs>_<end>.json` |
| 3 | `python score_baselines.py` | `baseline_scores.json` |
| 4 | `python oracle.py` | `oracle.json` (validation-tuned upper bound) |
| 5 | `python fig8_baselines.py <out_dir>` | Fig. 8, `rule_ci.json` |
| 6 | `python make_tables.py` then `python table3_md.py` | `tables.json` (Table S20), `table3.md` (Table 3) |
| 7 | `python ccdc_breaks.py` | break counts quoted in Section 5.1 |
| 8 | `python rule_at_points.py` | `rule_points_2025.json`: both branches of the rule at the 520 points of 2025 (Earth Engine) |
| 9 | `python branch_check.py` | `branch_check.json`: accuracy of each branch alone, 2025 and 2010 |

Steps 3 to 7 and 9 need only the data deposit and the JSON files of steps 1, 2 and 8, not Earth Engine. The files used in the paper are in `outputs/`; copy them next to the scripts to rerun steps 3 to 7 and 9.

## Parameters

LandTrendr: one composite per year, the 25th percentile of per-scene NBR from 1 November (Y − 1) to 31 May (Y), multiplied by −1000. Parameter sets are listed in `run_landtrendr.py`.

CCDC: break bands GREEN, RED, NIR, SWIR1 and SWIR2; tmask bands GREEN and SWIR1; chi-square probability 0.99; minNumOfYearsScaler 1.33; lambda 0.002; minObservations 6 and 3.

Each threshold is tuned on the calibration points by area-weighted F1. Bootstrap seed 20261026 with 4,000 stratified resamples.
