# Manuscript-to-code map

Each row lists the script that produces a manuscript or Supplementary Information (SI) item, the inputs it reads and the files it writes. Run the scripts in the order given in [README.md](README.md).

## Inputs

| Item | Code | Input | Output |
| --- | --- | --- | --- |
| 10 km analysis panel (Section 3.3; SI Appendix C.1) | `code/raster_to_panel.py` | `data/raster data/` (LULC, water, energy and food rasters; `adBoundary/YRB/YRB.shp`) | `grid_panel_10km.csv` |
| City-grid weights and city reference (SI Appendix A.4, Eq. A5) | `code/raster_to_panel.py` | `data/raster data/adBoundary/YRB/YRB.shp`, `YRB_City.shp` | `city_grid_weights_10km.csv`, `city_reference.csv` |
| Checks of the rebuilt inputs | `code/raster_to_panel.py` | As above and the CSVs in `data/` | `results/raster_to_panel/comparison_with_reference.csv`, `qc_summary.csv`, `input_manifest.csv` |
| City panel (SI Appendix A.4, Eq. A6) | `code/city_aggregation.py` | `data/grid_panel_10km.csv`, `data/city_grid_weights_10km.csv`, `data/city_reference.csv` | `data/city_panel_10km.csv` |

`raster_to_panel.py` writes the three rebuilt CSVs to `results/raster_to_panel/` by default and to `data/` with `--write-data`.

## Manuscript

| Item | Code | Input | Output |
| --- | --- | --- | --- |
| Grid and city scores (Sections 2.1-2.2; SI Appendices A.1-A.4) | `code/dea_performance.py --scale both` | Grid and city panels | `results/performance/grid_{water,energy,food}.csv`, `city_{water,energy,food}.csv`, `grid_references.csv`, `city_references.csv` |
| Section 4.1, Figs. 3-5 (grid and city scores, 2010-2020) | `code/analysis_4_1.py` | `results/performance/`, `data/cartography_layers.csv` | `figures/wef_grid_city_performance_3years_combined.*` |
| Section 4.2, Fig. 6 | `code/analysis_4_2.py` | `results/performance/grid_*.csv` | `results/analysis_4_2/regional_scores.csv`, `figures/Fig_6.*` |
| Section 4.3, Fig. 7 (SI Appendices B.1-B.2) | `code/analysis_4_3.py` | `results/performance/` | `results/analysis_4_3/grid_sde_parameters.csv`, `grid_sde_displacement_2010_2020.csv`, `city_sde_parameters.csv`, `city_sde_displacement_2010_2020.csv`, `city_common_sample.csv`; `figures/Fig_7.*` |

## Supplementary Information

| Item | Code | Input | Output |
| --- | --- | --- | --- |
| Appendix B.3, Table B1 | `code/specification_checks.py` | `results/performance/`, grid and city panels | `results/appendix_e/B1_sde_decomposition.csv`, `sde_city_baseline.csv` |
| Appendix C.1, ET coverage | `code/specification_checks.py` | Grid panel | `results/appendix_e/C1_et_coverage.csv` |
| Appendix C.2, Table C1 | `code/appendix_c.py` | `results/performance/grid_*.csv` | `results/appendix_c/Table_C1.csv` |
| Appendix C.3-C.5, Figs. C1-C3 | `code/appendix_c.py` | `results/performance/grid_*.csv`, `data/cartography_layers.csv` | `results/appendix_c/Fig_C1_statistics.csv`, `figures/Fig_C1.*`-`Fig_C3.*` |
| Appendix D.1, design | `code/appendix_d.py` | - | `results/appendix_d/method_decisions.csv`, `run_config.json` |
| Appendix D.2, Table D1 | `code/appendix_d.py --stages thresholds` | Grid panel | `results/appendix_d/Table_D1.csv`, `Table_D1_long.csv`, `frontier_references.csv`, `scores/<rule>_<system>.csv` |
| Appendix D.3, Table D2 and Figs. D1-D6 | `code/appendix_d.py --stages thresholds` | Grid panel, `results/performance/` | `results/appendix_d/Table_D2.csv`, `grid_rank_correlations.csv`, `regional_medians_2020.csv`, `regional_order_retention.csv`, `sde_parameters.csv`, `sde_displacements.csv`, `sde_vs_dominant.csv`, `city_sde_parameters.csv`, `city_sde_displacements.csv`; `figures/Fig_D1.*`-`Fig_D6.*` |
| Appendix D.4, Fig. D7 | `code/appendix_d.py --stages subsample` | Grid panel, `results/performance/` | `results/appendix_d/subsample_replicate_metrics.csv`, `subsample_summary.csv`, `figures/Fig_D7.*` |
| Appendix D, values quoted in the text | `code/appendix_d.py --stages summary` | `results/appendix_d/` | `results/appendix_d/appendix_d_summary.md` |
| Appendix E.1, Table E1 | `code/specification_checks.py` | Grid panel, `results/performance/` | `results/appendix_e/E1_energy_coverage.csv` |
| Appendix E.2, Tables E2-E3 and Fig. E1 | `code/specification_checks.py` | Grid and city panels, `results/performance/` | `results/appendix_e/E2_ratio_correlations.csv`, `E2_water_erosion.csv`, `figures/Fig_E1.*` |
| Appendix E.3, returns to scale and disposability | `code/specification_checks.py` | As above | `results/appendix_e/E3_technology.csv`, `E3_zhengzhou.csv` |
| Appendix E.4, small grid cells and cities | `code/specification_checks.py` | As above | `results/appendix_e/E4_small_units.csv` |
| Appendix E.5, partly non-vegetated water DMUs | `code/specification_checks.py` | As above | `results/appendix_e/E5_nonvegetated.csv` |
| Appendix A.4, city DMUs that needed the fallback solver | `code/specification_checks.py` | `results/performance/city_*.csv` | `results/appendix_e/A4_solver.csv` |

## Optional

| Item | Code | Input | Output |
| --- | --- | --- | --- |
| City scores as shapefiles for GIS mapping | `code/city_efficiency_csv_to_shapefile_2010_2015_2020.py` | `results/performance/city_*.csv`, `data/raster data/adBoundary/YRB/YRB_City.shp`, `YRB.shp` | `results/performance/city_shapefiles/` |

`code/common.py` holds the shared definitions (years, systems, land-use scopes, assignment rules, CRS and CSV helpers), and `code/map_utils.py` reads the map layers in `data/cartography_layers.csv`. Every analysis script imports these definitions.
