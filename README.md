# jhum-mapper

Code for mapping annual shifting-cultivation (jhum) clearing in Mizoram, north-east India,
1989–2025, from harmonized Landsat 5/7/8/9 in Google Earth Engine, with a 2025 Sentinel-2
comparison.

Accompanies: Khiangte, L.,Zoramkhuma, Sati, V.P.  (2026) *Anniversary-composite detection of annual shifting cultivation: compositing symmetry, dual-sensor stratification and
37 years of Landsat.*
Submitted.

Reference data are deposited separately as a dataset; see **Reference data** below.

## The rule

A pixel is mapped as current-year jhum when either branch fires.

| Branch | Window | Condition |
|---|---|---|
| Wide | 1 Nov – 31 May against the same months one year earlier | NBR_post ≤ 0.26 **and** dNBR ≥ 0.20 |
| Narrow | 1 Mar – 31 May against the same months one year earlier | NBR_post ≤ 0.14 **and** dNBR ≥ 0.46 |

Both sides of each pair are the 25th percentile of per-scene NBR. Applying the same
percentile to both sides is what suppresses seasonal browning of deciduous and bamboo
canopies; the paper shows that this symmetry, rather than the calendar alignment of the
windows on its own, is the principal control on commission error.

TM and ETM+ reflectance are put on the OLI scale with the coefficients of Roy et al. (2016):
near-infrared slope 0.8462, intercept 0.0412; shortwave-infrared 2 slope 0.9071, intercept
0.0172. Water is masked with JRC Global Surface Water occurrence ≥ 50 per cent. No minimum
mapping unit is applied to the published maps.

**Sentinel-2** uses B8A and B12 on the 20 m grid with Cloud Score+ `cs_cdf` ≥ 0.60, and a wide
branch refitted on the same 400 calibration points: NBR_post ≤ 0.27, dNBR ≥ 0.25. The narrow
branch is carried over from Landsat unchanged. Do not run the Landsat wide branch on
Sentinel-2; it reproduces none of the published Sentinel-2 numbers.

## Results this code reproduces (2025)

| | Landsat 30 m | Sentinel-2 20 m |
|---|---|---|
| Mapped area | 528.2 km² | 634.1 km² |
| Overall accuracy | 0.981 [0.974, 0.987] | 0.987 [0.980, 0.993] |
| User's accuracy | 0.672 [0.563, 0.774] | 0.739 [0.641, 0.835] |
| Producer's accuracy | 0.614 [0.491, 0.744] | 0.811 [0.653, 0.962] |
| F1 | 0.642 [0.542, 0.735] | 0.774 [0.672, 0.868] |

Error-adjusted jhum area, a property of the reference sample and therefore shared by both
maps: **578 km², 95 % CI 461–724**. Intervals are percentile intervals from a stratified
bootstrap of 4,000 resamples. Overall accuracy is uninformative at this prevalence: reference
jhum is about 2.7 per cent of the state, so a map that detected nothing would score 0.973.

Annual mapped area falls from a mean of 975 km² (1998–2009) to 443 km² (2020–2025);
Mann–Kendall τ = −0.513, p = 0.0001, Theil–Sen slope −20.8 km² yr⁻¹.

The 2025 figures above supersede the 310-point, three-stratum round published in v1.0.0 of
this repository (user's 0.672 against 0.886, producer's 0.614 against 0.901). The map is
unchanged; the sample, the strata and the reference labels are not (Section 2.3 of the paper). The earlier design stratified on the Landsat map
alone, so clearing that map missed could only be found in a thin buffer around it. The
scripts that produced it are kept under `gee/superseded/`, because the paper reports those
numbers as the comparison case.

## Layout

```
gee/                 Earth Engine scripts, in the order they were run
  00–05              mapping, calibration and the threshold, confuser and compositing tests
  06                 Sentinel-2 implementation of the rule (refit wide branch)
  07                 the 2025 validation design, round 1: five strata on both maps
  08                 round 2, an enlargement of strata D and E, later withdrawn
  09–18              annual series, portable tool, viewers, 2010 validation,
                     sensor consistency and transition, regional transfer, Fig 9 export
  19                 area and check points of the preliminary v3 map
  20                 round 3, the enlargement of strata D and E in the published sample
  superseded/        the earlier single-map validation design and Sentinel-2 test.
                     They produced numbers the paper still reports as the comparison
                     case (Table S2, Table S2c); they are not the published design.
python/              figure scripts (Matplotlib); see Running the figure scripts below
figure_data/         CSV inputs to the figures; figure_data/README.md maps each file
                     to its figure in the paper
colab/               the original threshold fitting on the 400 calibration points
docs/                interpretation guide for the 2025 dual-map sample
```

## Running the Earth Engine scripts

Open each in the Earth Engine Code Editor and press Run; exports appear in the Tasks tab.
Scripts use the FAO GAUL 2015 level-1 boundary as the sampling frame, not a local shapefile,
so nothing needs uploading. Two notes:

- `gee/07`, `gee/08` and `gee/20` print `STRATUM AREAS (km2)` before sampling. The final
  strata used in the paper are 425.6, 102.5, 208.5, 8391.8 and 12055.6 km²
  (`stratum_areas_2025.csv` in the dataset record). The script headers quote 8388.3 and
  12059.2 from an earlier run; D + E together agree to 0.1 km². If your printout differs
  materially, something upstream has changed and new points cannot be pooled with the
  published ones.
- The strata preview layer is off by default. Turning it on at state-wide zoom overruns the
  interactive memory limit, because `reproject` pins a 300 m distance kernel to the 30 m grid
  across 21,000 km². Those are tile errors in the map pane and do not affect the exports.

### Private Earth Engine assets

Six scripts read assets in the author's own Earth Engine project and will fail with "asset
not found" for anyone else. The paths are left in place because they record exactly what was
run; every such line is marked `// PRIVATE ASSET`. What to put there instead:

| Script | Asset | Replace with |
|---|---|---|
| `gee/01`, `gee/11` | `Jhum_1990toJhum_2025` | Your own digitised jhum polygons. Only used to draw training samples and to step through them visually. |
| `gee/02`, `gee/03` | `MIZ_Jhum_Samples_..._v3` | The table `gee/01` exports. Run `gee/01` first. |
| `gee/04` | `Label_2025` | `calibration_points_400_2025.csv` from the dataset record, uploaded as a table asset. |
| `gee/09` | `ASSET_DIR` | Any asset folder you can write to; it is only the export destination. |

Nothing in the published map depends on these. `gee/00`, `gee/05`, `gee/06`, `gee/07`,
`gee/08` and `gee/20` — the rule, the Sentinel-2 implementation and the validation samplers —
run with no asset of your own, on the FAO GAUL boundary.

## Running the figure scripts

Figures 2, 4, 5, 6, 7 and 8 need only `figure_data/` and run straight from a clone. Figures 1,
3 and 9 need rasters and vectors too large to ship here; `python/RUN_FIGURES_1_3_9.md` lists
every input, where to download it and the order to run.

```bash
pip install -r python/requirements.txt

# 1. Earth Engine: run python/gee_exports/*.js and download the GeoTIFFs
# 2. Boundaries and relief
python python/fig1/prep_districts_official.py <district_shapefile> layers_palsar/ layers30/
python python/fig1/prep_boundaries_relief.py  <geoboundaries_and_copernicus_dir> layers/
python python/fig1/prep_relief_30m.py         <geoboundaries_and_copernicus_dir> layers30/ 30
python python/fig1/prep_relief_palsar.py      <palsar_dem.tif> <copernicus_dem_30m.tif> \
                                              layers_palsar/mizoram.gpkg layers_palsar/
# 3. Figures
python python/fig1/plot_fig1.py  layers_palsar/ figure_data/ Fig1
python python/fig2_4-8/plot_line_art.py figure_data/ out/
python python/fig3/plot_fig3.py  <planet_and_sentinel_tiles> <landsat_chips> Fig3
python python/fig9/prep_fig9.py  <gee_exports> fig9_layers/ layers30/mizoram.gpkg
python python/fig9/plot_fig9.py  fig9_layers/ layers30/ <gee_exports> Fig9
pdftoppm -r 500 -png -singlefile Fig9.pdf Fig9_500dpi
python python/graphical_abstract/plot_graphical_abstract.py figure_data/ fig9_layers/ \
                                              layers_palsar/ <gee_exports> Graphical_abstract
```

Relief: Figure 1 uses the ALOS PALSAR 12.5 m DEM inside Mizoram, converted to orthometric
height with a +48.5 m geoid offset, and Copernicus GLO-30 beyond the state. Figure 9 and the
graphical abstract use GLO-30 at 30 m. In the published figures, SRTM GL1 v003 is used for
one thing only: the ≥ 10° slope constraint in the regional transfer test (`gee/17`). It also
appears as a hillshade in `gee/18` and `python/gee_exports/fig1_fig9_export_layers.js`, which
export layers for a QGIS route that was built first and then not used; both carry a header
saying so.

## Reference data

The reference points, the annual area series and the experiment inputs are deposited as a
separate dataset under CC BY 4.0, because the licences differ:

> Khiangte, L., 2026. *Reference points, annual area series and experiment inputs for annual
> jhum mapping in Mizoram (1989–2025)*, version 1.2.0 [dataset]. Zenodo.
> DOI: see the dataset record.

`figure_data/` duplicates the parts of that dataset the figure scripts need as inputs, so the
figures run without downloading the dataset. The dataset record is authoritative.

PlanetScope imagery used for reference interpretation cannot be redistributed under its
licence terms.

## How the figures were produced

All figures are data visualizations rendered by the Python scripts in `python/` directly from
the data described in the paper. No generative image model was used to create or alter any
image, figure or map; satellite image chips were contrast-stretched only. The Python code was
drafted with the assistance of a generative AI tool, as declared in the paper, and was
reviewed, tested and run by the author, who takes responsibility for it.

## Licence

Code: MIT, see `LICENSE`. Reference data: CC BY 4.0, in the dataset record.

## Citation

See `CITATION.cff`. Please cite the paper, the software record and the dataset record.
