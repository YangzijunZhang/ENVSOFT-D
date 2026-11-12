# Public input data

All public inputs are processed CSV files. No TIFF is required. Headers contain units and missing values are encoded as `NA`.

## `grid_panel_10km.csv`

Balanced panel of 8,622 in-basin grid cells for 2010, 2015 and 2020 (25,866 rows). Important fields are:

| Field | Meaning |
| --- | --- |
| `dmu_id`, `grid_id`, `year` | Decision-making unit, stable cell identifier and benchmark year |
| `grid_row`, `grid_col` | Position in the fixed 10 km grid |
| `x_center_m`, `y_center_m` | Albers equal-area cell centre in metres |
| `longitude_wgs84`, `latitude_wgs84` | Cell centre in WGS 84 |
| `boundary_fraction` | Fraction of the grid cell inside the basin |
| `area_inside_yrb_km2` | Cell area inside the basin |
| `lulc_dominant_*` | Primary dominant-class assignment |
| `lulc_threshold_50_*`, `lulc_threshold_70_*` | Sensitivity assignments at the stated minimum class share |
| `lulc_fraction_1_*` … `lulc_fraction_9_*` | Fractions for the nine CLCD classes |
| `water_input_et_volume_m3_yr` | Annual evapotranspiration volume, m3 yr-1 |
| `water_desirable_npp_t_c_yr` | Annual net primary production, t C yr-1 |
| `water_undesirable_soil_erosion_t_yr` | Annual soil erosion, t yr-1 |
| `energy_input_electricity_kwh_yr` | Annual electricity consumption, kWh yr-1 |
| `energy_desirable_gdp_2017_int_usd_yr` | Annual GDP, constant 2017 international USD yr-1 |
| `energy_undesirable_ghg_t_co2e_yr` | Annual greenhouse-gas emissions, t CO2e yr-1 |
| `food_input_nitrogen_kg_n_yr` | Annual nitrogen application, kg N yr-1 |
| `food_desirable_crop_production_kt_yr` | Annual production of four crops, kt yr-1 |
| `food_undesirable_n2o_t_n2o_yr` | Annual agricultural N2O emissions, t N2O yr-1 |

`water_input_et_volume_m3_yr` contains the annual MOD16A2 evapotranspiration volumes used in the Water analysis.

## `city_grid_weights_10km.csv`

Exact intersections between 76 clipped city geometries and in-basin 10 km grid cells. `allocation_weight` is the overlap area divided by the in-basin area of the grid cell. The same cell may intersect more than one city, but city weights for a cell do not exceed one in total.

## `city_reference.csv`

City identifiers, names, in-basin areas and centroids in the analysis Albers CRS. These centroids are used only for the city-scale standard-deviational ellipse.

## `cartography_layers.csv`

Processed vector geometries used only for figure backgrounds. Each row stores a layer name, feature identifier/name and WKT geometry in the same Albers equal-area CRS as the grid coordinates. The file contains the basin boundary, upper/middle/lower reaches, the Yellow River and the 76 city boundaries clipped to the basin. It replaces external TIFF and shapefile dependencies in the public plotting workflow.

## `city_panel_10km.csv` (generated)

`code/city_aggregation.py` creates this 76-city × 3-year panel from the three public CSVs above. It is the direct city-scale input to `code/dea_performance.py`. Coverage fractions are retained for every indicator so the aggregation remains auditable.

## Terms

The CSVs are processed research data derived from the sources listed in the manuscript and repository README. Redistribution and reuse must follow the licence and citation requirements of each source dataset. The repository's MIT License applies to code, not to third-party source rights.
