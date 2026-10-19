# figure_data

CSV inputs to the Matplotlib scripts in `python/`. The file names come from an earlier
figure order and no longer match the figure numbers in the paper. This table is the
mapping that matters.

| File | Read by | Figure in the paper |
|---|---|---|
| `fig3_featurespace_points.csv` | `fig2_4-8/plot_line_art.py`, `graphical_abstract/` | Figure 4, calibration feature space |
| `fig10_two_year_validation.csv` | `fig2_4-8/plot_line_art.py`, `graphical_abstract/` | Figure 5, accuracy under three sampling designs |
| `fig5_compositing.csv` | `fig2_4-8/plot_line_art.py` | Figure 6, compositing designs |
| `figure7a_accuracy.csv`, `figure7b_disagreement.csv` | `fig2_4-8/plot_line_art.py` | Figure 7, Landsat against Sentinel-2 |
| `fig4_annual_series.csv` | `fig2_4-8/plot_line_art.py`, `graphical_abstract/` | Figure 8, annual series |
| `fig1_validation_points.csv` | `fig1/plot_fig1.py` | Figure 1a, the 520 reference points |
| `fig11_validation_2010_points.csv` | `fig1/plot_fig1.py` | Figure 1b, the 160 reference points for 2010 |

Two notes on `fig10_two_year_validation.csv`. It holds three rows, not two: 2010 under the
three-stratum design, 2025 under the same three-stratum design, and 2025 under the published
five-stratum design. Figure 5 draws all three, because the point of that figure is that the
second and third assess the same map and differ only by design. The graphical abstract drops
the middle row and shows the two published assessments.

`fig4_annual_series.csv` gives 2025 as 528.57 km² in the scene projection; the paper reports
528.2 km², the same map measured on the declared EPSG:32646 30 m grid.

The full reference tables, including the point-level labels these summaries come from, are in
the separate dataset record under CC BY 4.0. The files here are the subset the figure scripts
read, so the figures run without downloading it.
