# Manuscript-to-code map

| Manuscript/SI item | Code | CSV input | Generated result |
| --- | --- | --- | --- |
| City aggregation (Appendix A.4) | `code/city_aggregation.py` | Grid panel, city weights, city reference | `data/city_panel_10km.csv` |
| Grid and city VRS-DDF performance | `code/dea_performance.py` | Grid and city panels | `results/performance/` |
| Appendix C, Table C1 and Figs. C1–C3 | `code/appendix_c.py` | Primary grid performance CSVs | `results/appendix_c/`, `figures/Fig_C1.*`–`Fig_C3.*` |
| Section 4.1 grid/city maps | `code/analysis_4_1.py` | Primary grid/city performance CSVs and cartography CSV | `figures/wef_grid_city_performance_3years_combined.*` |
| Section 4.2 and Fig. 6 | `code/analysis_4_2.py` | Primary grid performance CSVs | `results/analysis_4_2/`, `figures/Fig_6.*` |
| Section 4.3 and Fig. 7 | `code/analysis_4_3.py` | Primary grid/city performance CSVs | `results/analysis_4_3/`, `figures/Fig_7.*` |
| Appendix D land assignment sensitivity | `code/sensitivity.py` | Grid panel and primary scores | `results/sensitivity/`, `figures/Fig_D1.*`–`Fig_D6.*` |
| Appendix D frontier robustness | `code/robustness.py` | Grid panel and primary scores | `results/robustness/`, `figures/Fig_D7.*` |

`code/common.py` contains the shared manuscript definitions and CSV helpers. `code/map_utils.py` reads the shared processed cartography CSV. Every analytical script imports these same definitions.
