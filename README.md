# A spatially explicit framework for evaluating resource performance across heterogeneous land-use grid cells

This repository contains the processed CSV inputs and concise Python code corresponding to the revised manuscript and Supplementary Information. 

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

