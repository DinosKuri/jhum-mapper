# figure_data

CSV inputs to the Matplotlib scripts in `python/`. Since v1.3.0 each file name carries the number
of the figure it feeds in the current manuscript.

| File | Read by | Figure in the paper |
|---|---|---|
| `fig1a_validation_points_2025.csv` | `fig1/plot_fig1.py` | Figure 1a, the 520 reference points |
| `fig1b_validation_points_2010.csv` | `fig1/plot_fig1.py` | Figure 1b, the 160 reference points for 2010 |
| `fig3_compositing.csv` | `line_art/plot_line_art.py`, `graphical_abstract/` | Figure 3, compositing designs |
| `fig4_validation_designs.csv` | `line_art/plot_line_art.py`, `graphical_abstract/` | Figure 4, accuracy under three sampling designs |
| `fig5a_sensor_accuracy.csv`, `fig5b_sensor_disagreement.csv` | `line_art/plot_line_art.py` | Figure 5, Landsat against Sentinel-2 |
| `fig6_annual_series.csv` | `line_art/plot_line_art.py`, `graphical_abstract/` | Figure 6, annual series |
| `figS2_featurespace_points.csv` | `line_art/plot_line_art.py` | Figure S2, calibration feature space |

Figure 8 (LandTrendr and CCDC against the rule) reads the baseline outputs deposited in the dataset
record under `baselines/`; see `python/baselines/README.md`.

Two notes on `fig4_validation_designs.csv`. It holds three rows, not two: 2010 under the
three-stratum design, 2025 under the same three-stratum design, and 2025 under the published
five-stratum design. Figure 4 draws all three, because the point of that figure is that the
second and third assess the same map and differ only by design. The graphical abstract drops
the middle row and shows the two published assessments.

`fig6_annual_series.csv` gives 2025 as 528.57 km² in the scene projection; the paper reports
528.2 km², the same map measured on the declared EPSG:32646 30 m grid.

The full reference tables, including the point-level labels these summaries come from, are in
the separate dataset record under CC BY 4.0. The files here are the subset the figure scripts
read, so the figures run without downloading it.
