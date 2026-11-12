# A spatially explicit framework for evaluating resource performance across heterogeneous land-use grid cells

Data and code for the manuscript ENVSOFT-D-26-02569 (*Environmental Modelling & Software*) by Yangzijun Zhang and Wenyan Wu.

The repository contains the 10 km source rasters, the basin and city boundaries, the CSV inputs derived from them, and the Python code that reproduces the numbers, tables and figures in the manuscript and the Supplementary Information (SI). All paths are relative to the repository root.

- Source code: <https://github.com/YangzijunZhang/ENVSOFT-D>
- Archived release: v1.0.0 on Zenodo (the DOI is added here after archiving; see [RELEASE.md](RELEASE.md))
- Licence: MIT for the code. The data keep the licences of their sources (see [Data and licences](#data-and-licences)).

## Repository structure

```text
data/
  raster data/                    10 km source rasters for 2010, 2015 and 2020
    LULC/                         CLCD class fractions, dominant class and share, 50/70/80% assignments
    water/  energy/  food/        the nine input and output variables as 10 km cell totals
    adBoundary/YRB/               basin, reach, river and city boundaries (shapefiles)
  grid_panel_10km.csv             10 km analysis panel (8,622 cells x 3 years)
  city_grid_weights_10km.csv      city-grid intersection weights
  city_reference.csv              76 cities: identifiers, in-basin areas and centroids
  city_panel_10km.csv             76-city panel (created by city_aggregation.py)
  cartography_layers.csv          simplified map layers (WKT) used as figure backgrounds
code/
  raster_to_panel.py              rasters and boundaries -> grid panel, city weights, city reference
  city_aggregation.py             grid panel -> 76-city panel
  dea_performance.py              grid- and city-scale DEA scores
  analysis_4_1.py                 Section 4.1 maps
  analysis_4_2.py                 Section 4.2 and Fig. 6
  analysis_4_3.py                 Section 4.3 and Fig. 7
  appendix_c.py                   SI Table C1 and Figs. C1-C3
  appendix_d.py                   SI Appendix D: Tables D1-D2 and Figs. D1-D7
  specification_checks.py         SI Appendix B.3 and Appendix E
  city_efficiency_csv_to_shapefile_2010_2015_2020.py
                                  optional: city scores as shapefiles for GIS mapping
  common.py                       shared definitions and CSV helpers (no analysis)
  map_utils.py                    shared plotting helpers
results/                          generated tables (CSV)
figures/                          generated figures (PNG and PDF)
```

[MANIFEST.md](MANIFEST.md) maps each manuscript and SI item to its script, inputs and outputs. [METHODS.md](METHODS.md) summarises the computations, and [REPRODUCIBILITY.md](REPRODUCIBILITY.md) describes the checks and the expected numerical differences.

## Installation

The analyses were run with Python 3.14.3 on macOS (Apple silicon).

```bash
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
```

Alternatively, use conda: `conda env create -f environment.yml`.

## Reproducing the results

Run the scripts from the repository root in this order:

```bash
python code/raster_to_panel.py                 # 1. rebuild the CSV inputs from the rasters; compare with data/
python code/city_aggregation.py                # 2. 76-city panel
python code/dea_performance.py --scale both    # 3. grid and city scores
python code/appendix_c.py                      # 4. Table C1 and Figs. C1-C3
python code/analysis_4_1.py                    # 5. Section 4.1 maps
python code/analysis_4_2.py                    # 6. Section 4.2 and Fig. 6
python code/analysis_4_3.py                    # 7. Section 4.3 and Fig. 7
python code/specification_checks.py            # 8. Appendix B.3 and Appendix E
python code/appendix_d.py --workers 4          # 9. Appendix D
```

- Step 1 writes the rebuilt CSVs to `results/raster_to_panel/` together with a column-by-column comparison with the files in `data/`. It does not modify `data/` unless `--write-data` is given.
- Steps 2-9 read CSV files only.
- Step 9 re-estimates the frontiers under two alternative land-use rules and for 100 subsamples. It takes about 2.5 hours on one core. `--workers` runs the replicates in parallel and gives identical results. `python code/appendix_d.py --stages thresholds` runs Sections D.2-D.3 only (a few minutes).
- Optional: `python code/city_efficiency_csv_to_shapefile_2010_2015_2020.py` writes the city scores as shapefiles to `results/performance/city_shapefiles/`.

## Model summary

- Directional distance function with the observed direction g = (x, y, b), variable returns to scale and weak disposability of the undesirable output. One frontier is pooled over 2010, 2015 and 2020. The score is TE = 1 / (1 + β).
- Water: group frontiers for cropland, forest and grassland and a meta-frontier across them; the meta-frontier scores are reported. Energy: impervious land. Food: cropland. Energy and food use their group-frontier scores.
- Variables enter the linear programmes in their physical units, without rescaling, winsorisation or imputation.
- The city benchmark aggregates grid totals of all land-use types to 76 cities and is estimated independently of the grid frontiers. Cross-year city results use the same cities in all three years.

See [METHODS.md](METHODS.md) for details.

## Data and licences

The 10 km rasters were derived from the source datasets listed in Table 1 of the manuscript. The conversion of each dataset to 10 km cell totals is described in SI Appendix C.1 and recorded in the metadata tags of each raster. Reuse of the rasters and CSVs must follow the licence and citation requirements of the source datasets. The boundary files in `data/raster data/adBoundary/` keep the licences of their sources. The MIT License of this repository applies to the code only.

## Citation

Please cite the article and the archived release. Citation metadata are given in [CITATION.cff](CITATION.cff).
