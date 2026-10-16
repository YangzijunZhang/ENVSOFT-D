# A spatially explicit framework for evaluating resource performance across heterogeneous land-use grid cells

This repository contains the processed CSV inputs and concise Python code corresponding to the revised manuscript and Supplementary Information. No TIFF file or machine-specific path is required.

## Files

```text
data/grid_panel_10km.csv          Processed 10 km analysis panel
data/city_grid_weights_10km.csv   City–grid intersection weights
data/city_reference.csv           The 76 city identifiers and centroids
data/cartography_layers.csv       Basin, reaches, river and city boundaries (WKT)
code/city_aggregation.py          Creates the initial 76-city CSV panel
code/dea_performance.py           One DEA file for grid and city performance
code/appendix_c.py                Table C1 and Figs. C1-C3
code/analysis_4_1.py              Section 4.1 grid/city maps
code/analysis_4_2.py              Section 4.2 regional cropland analysis
code/analysis_4_3.py              Section 4.3 grid/city SDE analysis
code/sensitivity.py               50% and 70% land-share sensitivity
code/robustness.py                100-replicate frontier robustness
code/common.py                    Shared definitions; no analysis is performed here
code/map_utils.py                 Shared CSV-based cartography helpers
```

## Execution order

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python code/city_aggregation.py
python code/dea_performance.py --scale both
python code/appendix_c.py
python code/analysis_4_1.py
python code/analysis_4_2.py
python code/analysis_4_3.py
python code/sensitivity.py
python code/robustness.py
```

The scripts write generated performance tables and supplementary statistics to `results/`, and the manuscript/SI figures to `figures/`.

## Implemented methods

- Three pooled benchmark years: 2010, 2015 and 2020.
- Observed direction `g = (x_o, y_o, b_o)`, variable returns to scale, and weak disposability of the undesirable output.
- Score `TE = 1 / (1 + beta)`.
- Water group frontiers for Cropland, Forest and Grassland plus their meta-frontier; Energy for Impervious land; Food for Cropland.
- Variables enter the linear programmes directly in their physical units, with no standardisation or rescaling.
- The city benchmark uses all land-use classes and is estimated independently from the grid frontiers for all 76 cities.
- All maps use the processed WKT geometries in `data/cartography_layers.csv`; no TIFF, shapefile or machine-specific path is required.
- Section 4.1 maps group-frontier grid TE and independently estimated city TE for all three systems and years.
- Section 4.2 reports 2020 Cropland group-frontier scores for Regions A–D.
- Section 4.3 uses every evaluated DMU, weighted by its score, with one-standard-deviation SDE semi-axes and reports the city/reach containing each centre.
- Sensitivity re-estimates the same model independently under the 50% and 70% rules and produces Table D1 and Figs. D1–D6.
- Robustness uses 100 replicates, 80% sampling without replacement within year × land-type strata, seed `2602569`, and produces Fig. D7.

Repository URL: <https://github.com/YangzijunZhang/ENVSOFT-D>. See [RELEASE.md](RELEASE.md) for the GitHub–Zenodo release steps.
