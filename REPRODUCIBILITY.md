# Reproducibility guide

Run each command from the repository root in the order shown in `README.md`.

`city_aggregation.py` creates `data/city_panel_10km.csv` from the public grid panel, city–grid weights and city reference CSV. `dea_performance.py` then estimates both spatial scales and writes one result CSV per resource system and scale.

The remaining scripts have deliberately narrow responsibilities:

- `appendix_c.py`: Table C1 and Figs. C1–C3.
- `analysis_4_1.py`: combined grid/city maps for all systems and benchmark years.
- `analysis_4_2.py`: regional Cropland summary and Fig. 6.
- `analysis_4_3.py`: grid/city SDE parameters, displacement and Fig. 7.
- `sensitivity.py`: independently re-estimated 50% and 70% scenarios, Tables D1–D2 and Figs. D1–D6.
- `robustness.py`: fixed-seed stratified subsampling and Fig. D7 using the same DEA functions.

All paths are relative. Generated CSVs encode missing values as `NA` and numbers with 12 significant digits. Package versions are pinned in `requirements.txt` and `environment.yml`. Solver-success and residual columns provide row-level numerical checks for generated outputs.

Map layers are supplied in `data/cartography_layers.csv` as projected WKT geometries. This single CSV contains the YRB boundary, upper/middle/lower reaches, the Yellow River and all 76 clipped city boundaries, so the plotting workflow does not depend on TIFF or shapefile inputs.
