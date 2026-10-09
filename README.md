# jhum-mapper

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23269371.svg)](https://doi.org/10.5281/zenodo.23269371)

Annual mapping of shifting-cultivation (*jhum*) clearing from harmonized Landsat 5/7/8/9,
calibrated and independently validated for Mizoram, north-east India, 1989–2025.

## The rule

A pixel is current-year jhum when either branch is true:

| Branch | Condition | Windows |
|---|---|---|
| A (wide) | `NBR_post ≤ 0.26` and `dNBR ≥ 0.20` | 1 Nov – 31 May, against the same months one year earlier |
| B (narrow) | `NBR_post ≤ 0.14` and `dNBR ≥ 0.46` | 1 Mar – 31 May, against the same months one year earlier |

`NBR = (NIR − SWIR2) / (NIR + SWIR2)`, `dNBR = NBR_pre − NBR_post`. Every composite is the 25th
percentile of per-scene NBR, applied on **both** sides of each comparison. TM and ETM+ are
harmonized to OLI (Roy et al., 2016).

## Accuracy

Blind stratified samples drawn from the final map (strata: mapped jhum, within 300 m of mapped
jhum, beyond 300 m); estimators of Olofsson et al. (2014) and Stehman (2014).

| | 2010 (TM / ETM+) | 2025 (OLI / OLI-2) |
|---|---|---|
| Reference points | 160 | 310 |
| Overall accuracy | 0.996 ± 0.008 | 0.995 ± 0.005 |
| User's accuracy (jhum) | 0.983 ± 0.033 | 0.886 ± 0.075 |
| Producer's accuracy (jhum) | 0.812 ± 0.299 | 0.901 ± 0.175 |
| F1 | 0.890 | 0.893 |
| Error-adjusted area | 454 ± 168 km² | 518 ± 108 km² |

Reference imagery: PlanetScope monthly basemaps (2025); Google Earth Pro and Landsat (2010).

## Quick start

1. Open the [Earth Engine Code Editor](https://code.earthengine.google.com).
2. Paste `gee/03_jhum_mapper_tool.js` into a new script and set the area and year:
   ```javascript
   var ADMIN_NAME = 'Mizoram';
   var COUNTRY    = 'India';
   var YEAR       = 2025;
   ```
3. Run. The Console prints an area report (Mizoram 2025 with the defaults: about 528 km²);
   the Tasks tab exports a GeoTIFF and the report.

**Elsewhere:** set the local jhum calendar (`CLEAR_START`, `BURN_START`, `BURN_END`) and check
known plots before trusting the output. Thresholds apply to Landsat surface reflectance.

## Contents

```
gee/                 Earth Engine scripts
data/                Reference points and the annual area series
data/experiments/    Inputs to the sensor, Sentinel-2 and transfer experiments
docs/                Reference interpretation protocol
```

| Script | Purpose |
|---|---|
| `01_final_map.js` | Final validated map for any year or the full series, area table, inspector |
| `02_annual_series_fallow_metrics.js` | 1989–2025 stack, clearing frequency, fallow age, cycle length |
| `03_jhum_mapper_tool.js` | Standalone tool: any area, year and local calendar |
| `04_validation_sample_2025.js` | Stratified blind validation sample, 2025 |
| `05_validation_extension_2025.js` | Additional points per stratum, pooled with the first round |
| `06_validation_sample_2010.js` | Stratified blind validation sample, 2010 |
| `07_blind_viewer.js` | Imagery only, for labelling (no map, no index values) |
| `08_compositing_experiment.js` | Compositing designs evaluated at the calibration points |
| `09_sentinel2_comparison.js` | Sentinel-2 replication under four cloud-mask settings |
| `10_sensor_consistency.js` | NBR agreement between sensor pairs over stable forest |
| `11_sensor_transition_test.js` | Transition years rebuilt with and without each sensor |
| `12_transfer_test_ne_india.js` | Exploratory application to north-east India (plausibility only) |

Scripts 04, 05, 08 and 09 read labelled points as Earth Engine table assets
(`projects/YOUR_PROJECT/assets/...`). Upload `data/calibration_points_400_2025.csv` or
`data/validation_points_310_2025.csv` and edit the path. Annual rasters are not stored here;
`01` and `02` regenerate them.

## Known limits

- Accuracy was measured for 2010 and 2025; other years rely on identical processing.
- The mapped class is dry-season clearing by fire on forested slopes; plots cleared without
  burning are missed and a few non-cultivation fires are included.
- Deciduous leaf-off is the main source of commission (four of eight errors in 2025).
- Six early years (1989, 1991–92, 1995–97) have 65–90 % image coverage.
- Mapped area at sensor transitions depends on observation density (up to 47 % in 2013); the
  1998–2025 trend is robust to excluding those years.

## Citation

Please cite the archived release (also available from GitHub's *Cite this repository*):

> Khiangte, L. (2026). *jhum-mapper: code and reference data for annual jhum mapping in Mizoram
> (1989–2025)* (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23269371

The accompanying paper is under review; its reference will be added on publication.

## Licence

Code: MIT (`LICENSE`). Data and documentation: CC BY 4.0 (`LICENSE-data`).
PlanetScope imagery is not redistributed (Planet Education and Research programme). Landsat data
courtesy of the U.S. Geological Survey; processing in Google Earth Engine.
