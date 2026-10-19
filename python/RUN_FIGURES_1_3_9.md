# Running Figures 1, 3 and 9

Figures 2, 4, 5, 6, 7 and 8 need only `figure_data/` and run anywhere:

```bash
python python/fig2_4-8/plot_line_art.py figure_data/ out/
```

Figures 1, 3 and 9 need raster and vector data that is too large to ship with the code. This
file lists exactly what each one needs, where to get it, and the order to run.

```bash
pip install -r python/requirements.txt
```

A serif font must be installed. Times New Roman is used if present; Liberation Serif is
metric-identical and is the fallback (`apt install fonts-liberation`).

---

## Figure 1 — study area, relief and the reference points

### What you need

| Input | Source | Notes |
|---|---|---|
| `geoBoundaries-IND-ADM1.geojson` | geoboundaries.org | India states |
| `geoBoundaries-IND-ADM2.geojson` | geoboundaries.org | districts, clipped to Mizoram |
| `geoBoundaries-BGD-ADM0.geojson` | geoboundaries.org | Bangladesh outline |
| `geoBoundaries-MMR-ADM1.geojson` | geoboundaries.org | Myanmar, dissolved |
| `Copernicus_DSM_*.tif` | Copernicus GLO-30, OpenTopography or the ESA portal | all tiles covering 92.0–93.7 E, 21.8–24.7 N plus a 14 km pad |
| ALOS PALSAR RTC DEM, 12.5 m | ASF DAAC Vertex, search "ALOS PALSAR RTC high-res" over Mizoram | mosaic the tiles into one GeoTIFF first |
| Official district shapefile | your own copy, with LGD codes | optional; `prep_districts_official.py` uses it instead of geoBoundaries ADM2 |
| `figure_data/fig1_validation_points.csv` | in this repository | the 520 points |
| `figure_data/fig11_validation_2010_points.csv` | in this repository | the 160 points |

Put all four GeoJSONs and all the Copernicus tiles in one folder; call it `src/`.

### Run

```bash
mkdir -p layers30 layers_palsar

# 1. vectors + a 30 m Copernicus relief (also used by Figure 9)
python python/fig1/prep_relief_30m.py src/ layers30/ 30

# 2. official districts, if you have the shapefile; writes into both folders
python python/fig1/prep_districts_official.py <district_shapefile> layers_palsar/ layers30/

# 3. the vector layers Figure 1 draws must sit next to the PALSAR relief
cp layers30/bgd.gpkg layers30/mmr.gpkg layers30/ind_ne_states.gpkg \
   layers30/india.gpkg layers_palsar/
#    if you skipped step 2, copy these two as well:
cp layers30/mizoram.gpkg layers30/districts.gpkg layers_palsar/

# 4. PALSAR relief inside the state, Copernicus beyond, +48.5 m geoid offset
python python/fig1/prep_relief_palsar.py <palsar_mosaic.tif> <copernicus_30m.tif> \
                                         layers_palsar/mizoram.gpkg layers_palsar/

# 5. draw it
python python/fig1/plot_fig1.py layers_palsar/ figure_data/ Fig1
```

`layers_palsar/` must end up holding: `mizoram.gpkg`, `districts.gpkg`, `bgd.gpkg`, `mmr.gpkg`,
`ind_ne_states.gpkg`, `india.gpkg`, `dem_25m.tif`, `hillshade_25m.tif`. If `plot_fig1.py` stops
with a file-not-found, that list is the thing to check.

Outputs `Fig1.pdf` and `Fig1_500dpi.tif`.

**What changed, and why this has to be re-run.** Panel a now draws the five strata of the
published design, panel b the three strata of 2010. The version from 11 October drew three
classes in both panels and no longer matches Figure 1's caption or Section 2.3. Stratum counts
and n are read from `figure_data/fig1_validation_points.csv` (the 520 points of rounds 1 and 3:
A 40, B 45, C 45, D 200, E 190).

The geoid offset is not cosmetic. PALSAR RTC heights are above the WGS84 ellipsoid; Copernicus
is above mean sea level. Over Mizoram the difference is +48.5 m, and without it the two
sources step against each other at the state line.

---

## Figure 3 — visual validation chips, four sites

### What you need

| Input | Source |
|---|---|
| PlanetScope monthly basemap chips, four sites × six months | Planet Explorer, under your Education and Research licence |
| Sentinel-2 chips, same sites and months | Copernicus Browser, or export them from Earth Engine |
| Landsat 8/9 chips | `python/gee_exports/fig3_landsat_chips.js` in this repository |

The Planet and Sentinel-2 chips are read by folder and file name, and the script has those
names hard-coded in its `SITES` dictionary, for example
`V2025_51_CONFIRMJHUM/NOV_2024.tif` and `V2025_51_CONFIRMJHUM_sentinnel/nov_24.tif`. Keep the
folder names exactly as they are in the script, or edit `SITES` to match your own.

The four sites are V2025_051 (detected, wide branch), V2025_030 (detected, narrow branch only),
V2025_015 (commission error, deciduous leaf-off) and the fourth listed in the script.

### Run

```bash
# Landsat chips: open python/gee_exports/fig3_landsat_chips.js in the Earth Engine Code
# Editor, run it, then download the FIG3_landsat_* GeoTIFFs from Drive into landsat_chips/

python python/fig3/plot_fig3.py <planet_and_sentinel_dir> landsat_chips/ Fig3
```

**Width.** `plot_fig3.py` lays the figure out at 190 mm, Elsevier's full-page maximum, with
every label at 6.5 pt or larger: long site labels wrap, the index line splits over two rows.
Nothing in the data changes with the width; only the layout constants at the top of the script.

---

## Figure 9 — the map of jhum, fallow age and two detail insets

### What you need

| Input | Source |
|---|---|
| `FIG9_jhum_2025.tif` | `gee/18_fig9_export_for_qgis.js` |
| `FIG9_fallow_age_1989_2025.tif` | same |
| `FIG9_jhum_count_1989_2025.tif` | same |
| `FIG9_cycle_mean_1998_2025.tif` | same |
| `FIG9_inset_A_jhum.tif`, `FIG9_inset_A_fallow_age.tif`, `FIG9_inset_A_landsat_rgb.tif` | same |
| `FIG9_inset_B_jhum.tif`, `FIG9_inset_B_fallow_age.tif`, `FIG9_inset_B_landsat_rgb.tif` | same |
| `layers30/` with `mizoram.gpkg`, `districts.gpkg`, `bgd.gpkg`, `mmr.gpkg`, `ind_ne_states.gpkg`, `hillshade_30m.tif` | step 1 of Figure 1 above |

Use **`gee/18`**, not `python/gee_exports/fig1_fig9_export_layers.js`. The two scripts export
different things: `gee/18` writes the individual named rasters the Python reads, while the
other writes multi-band stacks for a QGIS route that produced no published figure.

The hillshade comes from `layers30/hillshade_30m.tif`, built from Copernicus GLO-30 — not from
the `FIG9_hillshade` raster that `gee/18` also exports, which is SRTM. Table 1 of the paper
credits Copernicus for this figure, so take the hillshade from `layers30/`.

### Run

```bash
# 1. Earth Engine: run gee/18, then download every FIG9_* GeoTIFF from Drive into gee_exports/
# 2. clean the fallow-age raster and build the 1 km state summary
python python/fig9/prep_fig9.py gee_exports/ fig9_layers/ layers30/mizoram.gpkg

# 3. draw it
python python/fig9/plot_fig9.py fig9_layers/ layers30/ gee_exports/ Fig9

# 4. rasterise for the journal
pdftoppm -r 500 -png -singlefile Fig9.pdf Fig9_500dpi
```

`prep_fig9.py` prints the 2025 mapped area and the ever-cleared area. The mapped area should
read close to 528 km². It also asserts that fallow age is 0 wherever 2025 jhum is 1; if that
assertion fails, the rasters came from different runs and should be re-exported together.

**What changed, and why this has to be re-run.** The legend text in `plot_fig9.py` now reads
`Mapped area 528.2 km² / Error-adjusted 578 km² (461–724)`. Any PNG rendered before
19 October carries a superseded area and contradicts the abstract, Table 2 and Section 4.1.

---

## Graphical abstract

Needs everything above, and runs last:

```bash
python python/graphical_abstract/plot_graphical_abstract.py figure_data/ fig9_layers/ \
                                   layers_palsar/ gee_exports/ Graphical_abstract
```

It now shows the 2010 assessment and the published 2025 assessment, and skips the superseded
2025 round that `fig10_two_year_validation.csv` also carries.
