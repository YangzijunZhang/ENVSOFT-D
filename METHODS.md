# Computational methods

This file summarises what the code computes. The equations and their justification are given in Section 2 of the manuscript and in Appendices A-E of the Supplementary Information (SI).

## 1. From rasters to the analysis panel (`raster_to_panel.py`)

All rasters share one grid of 112 rows x 204 columns of 10 km cells in an Albers equal-area projection (standard parallels 25° N and 47° N, central meridian 105° E, WGS 84). The resource rasters hold annual totals for each cell, obtained from the source datasets by the conversions described in SI Appendix C.1. In cells crossing the basin boundary, the totals are already scaled to the in-basin part of the cell.

1. **In-basin fraction.** The basin polygon (`YRB.shp`) is projected to the grid CRS and intersected exactly with every cell. The in-basin fraction is stored in single precision, as in the original preprocessing. The 8,622 cells with a fraction greater than zero form the panel; 1,080 of them cross the basin boundary.
2. **Grid panel.** For each benchmark year (2010, 2015 and 2020), the script writes the cell identifiers, the cell-centre coordinates (Albers and WGS 84), the nine CLCD class fractions, the dominant class and its share, the 50%, 70% and 80% assignments and the nine resource variables. Raster values are copied unchanged. Missing values are written as `NA`.
3. **City tables.** Each of the 76 city polygons (`YRB_City.shp`) is clipped to the basin and intersected exactly with the grid. The allocation weight of cell *i* to city *c* is the overlap area divided by the in-basin area of the cell (SI Eq. A5). The city reference table holds the in-basin area and the centroid of each clipped city.

The MOD16A2 ET product provides no values for non-vegetated pixels. In each benchmark year, 69 cells have no valid ET and are therefore not evaluated in the water system (SI Appendix C.1).

## 2. Land-use assignment and evaluation scope

The primary analysis assigns each cell to its dominant CLCD class. The sensitivity analyses assign a cell only if its dominant class covers at least 50% or 70% of the valid LULC area of the cell; other cells are left unassigned.

The land-use scope of each resource system is fixed before estimation:

- Water: cropland, forest and grassland.
- Energy: impervious land.
- Food: cropland.

Only complete cases are evaluated. Inputs and desirable outputs must be greater than zero; undesirable outputs may be zero but not negative. A frontier is estimated only when at least nine eligible decision-making units (DMUs) are available.

## 3. Directional-distance DEA (`dea_performance.py`)

For DMU *o*, the pooled VRS directional-distance model maximises β:

```text
max beta

subject to
sum_j lambda_j x_j <= (1 - beta) x_o
sum_j lambda_j y_j >= (1 + beta) y_o
sum_j lambda_j b_j  = (1 - beta) b_o
sum_j lambda_j = 1
lambda_j >= 0
0 <= beta <= 1
```

The direction is the observed vector g = (x_o, y_o, b_o). The equality constraint on the undesirable output implements weak disposability, and the convexity constraint implements variable returns to scale. Observations from all three benchmark years form one intertemporal reference technology.

The reported score is

```text
TE = 1 / (1 + beta)
```

TE ranges from 0.5 to 1; higher values indicate stronger resource performance relative to the observed best practice. In the water system, group frontiers are estimated for cropland, forest and grassland, and the meta-frontier pools all three; the meta-frontier scores are used in the analysis. Energy and food each cover one land-use type and use their group-frontier scores.

The variables enter the linear programmes in the physical units given in the CSV headers, without standardisation, rescaling, winsorisation or imputation.

**Reference set and solver.** Exact duplicate input-output vectors and points inside the three-dimensional convex hull of the reference observations are omitted from the reference matrix, because they cannot change the VRS technology. The linear programmes are solved with HiGHS through `scipy.optimize.linprog`. If a solve fails, it is repeated with the full set of unique reference vectors, first with the default HiGHS method and then with the dual simplex (`highs-ds`) and interior-point (`highs-ipm`) methods. Each output row records whether the solve succeeded and the largest relative constraint residual.

## 4. City benchmark (`city_aggregation.py`, `dea_performance.py`)

City totals are allocated from grid totals with the exact intersection weights:

```text
w_ci         = overlap_area_ci / in_basin_area_i
city_total_c = sum_i (w_ci * grid_total_i)
```

Grid cells of all land-use classes contribute; no dominant-class filter is applied. For each resource system, a cell contributes to a city's totals only if all three variables of that system are valid in the cell, so the input and the outputs of a city cover the same area. Missing cells are not filled or rescaled; the covered share of each city's in-basin area is written as a diagnostic.

City inputs must be positive; desirable and undesirable city totals may be zero but not negative. For each resource system, the pooled city frontier is built from the cities that are eligible in all three benchmark years. Cross-year maps and SDEs use only cities with a successful, finite score in every benchmark year, so the same cities are compared in 2010, 2015 and 2020. Each resource system has its own city frontier, estimated independently of the grid-scale frontiers.

## 5. Standard-deviational ellipses (`analysis_4_3.py`)

For locations p_i and score weights TE_i, the weighted centre is `sum(TE_i p_i) / sum(TE_i)`, and the weighted covariance matrix uses the divisor `sum(TE_i)`. Its eigenvectors give the orientation of the ellipse, and the square roots of its ordered eigenvalues give the major and minor semi-axes. These are one-standard-deviation semi-axes; under a bivariate normal distribution the ellipse would contain about 39% of the total weight.

Grid SDEs use the water meta-frontier scores and the energy and food group-frontier scores of every evaluated DMU. City SDEs use the city scores of the common city sample. The displacement from 2010 to 2020 is reported as projected east and north components, geodesic distance on the WGS 84 ellipsoid, azimuth clockwise from north and an eight-point compass direction.

## 6. Sensitivity and robustness (`appendix_d.py`)

**Land-use assignment.** The frontiers are re-estimated independently under the 50% and 70% assignment rules with the same model and land-use scopes. Changes in DMU counts, score levels, rank correlations for cells evaluated under both rules, regional medians and SDE centres are compared with the primary analysis. City totals do not depend on the assignment rule, so the city results are not re-estimated.

**Subsampled frontiers.** With seed `2602569`, 100 replicates are drawn. Each replicate samples 80% of the eligible primary-analysis DMUs (rounded down) without replacement within each year x land-use stratum. The group frontiers and the water meta-frontier are re-estimated from each sample, and only the sampled DMUs are rescored. For each replicate, the Spearman rank correlation and the change in median score are computed against the scores of the same DMUs under the full frontier.

## 7. Specification diagnostics (`specification_checks.py`)

- **SI Appendix B.3.** The grid-scale displacement of each SDE centre is recalculated with only the cells evaluated in both 2010 and 2020, to separate score changes from changes in the set of evaluated cells.
- **SI Appendix E.** Coverage of the basin totals by the evaluated energy cells; rank correlations between the scores and the simple output-input ratios; water scores re-estimated without soil erosion; scores under constant returns to scale and under strong disposability of the undesirable output; exclusion of cells with less than 10 km² inside the basin and of cities with less than 100 km² inside the basin; and exclusion of water DMUs whose non-vegetated share exceeds 20%.
