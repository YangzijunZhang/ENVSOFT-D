# Computational methods

## Analysis panel and land assignment

The analysis uses one record per 10 km grid cell and benchmark year (2010, 2015 and 2020). Grid cells intersecting the Yellow River Basin are retained. The primary analysis assigns each cell to its dominant CLCD class. Sensitivity analyses require the assigned class to cover at least 50% or 70% of valid LULC area.

System scopes are fixed before DEA:

- Water: Cropland, Forest and Grassland.
- Energy: Impervious.
- Food: Cropland.

Rows are complete cases. Inputs and desirable outputs must be greater than zero; undesirable outputs may be zero but not negative. A frontier is estimated only when at least nine eligible decision-making units are available.

## Directional-distance DEA

For decision-making unit `o`, the pooled VRS directional-distance model maximises `beta`:

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

The direction is the observed vector `g = (x_o, y_o, b_o)`. Equality for the undesirable output implements weak disposability. The convexity constraint implements variable returns to scale. Observations from all three years form one intertemporal reference technology.

The bounded reporting score is:

```text
TE = 1 / (1 + beta)
```

`TE` ranges from 0.5 to 1, with higher values indicating stronger resource performance. The water meta-frontier pools all eligible Cropland, Forest and Grassland observations. Energy and Food each contain one land-use type and therefore use their group-frontier scores.

All variables are supplied to the linear programmes directly in the physical units shown in the CSV headers. There is no standardisation or rescaling. Exact duplicate input-output vectors and points lying inside the three-dimensional convex hull are omitted from the reference matrix because they cannot alter the VRS technology. If a reduced-reference solve fails, the code retries with every unique reference vector.

## City benchmark

Grid totals are allocated to cities using exact city–grid intersections:

```text
w_ci = overlap_area_ci / in_basin_grid_area_i
city_total_c = sum_i (w_ci * grid_total_i)
```

All land-cover classes are included and no dominant-class filter is applied. Finite grid-cell totals are summed using the intersection weights; missing cells are not filled or rescaled, and the covered-area fraction is written as a diagnostic. The city panel contains all 76 cities in each benchmark year. City inputs must be positive, while desirable and undesirable city totals may be zero but not negative. Each resource system uses its own pooled 228-observation city frontier, independent of the grid-scale frontiers.

## Standard-deviational ellipse

For locations `p_i` and performance weights `TE_i`, the weighted centre is `sum(TE_i p_i) / sum(TE_i)`. The weighted covariance matrix uses divisor `sum(TE_i)`. Its eigenvectors define orientation, and the square roots of its ordered eigenvalues define the major and minor semi-axes. These are one-weighted-standard-deviation semi-axes; under a bivariate normal distribution the ellipse would contain about 39% of total weight. Grid analyses use the water meta-frontier TE and the energy and food group-frontier TE. City analyses use the independently estimated city TE. Every evaluated decision-making unit is included.

Displacement is calculated from the 2010 and 2020 weighted centres. The published table includes projected east/north components, geodesic distance on WGS 84, azimuth and an eight-point compass direction.

## Robustness

The robustness analysis uses seed `2602569`. For each of 100 replicates, 80% of eligible primary-scenario observations (rounded down) are sampled without replacement within each year × land-type stratum. Pooled group and water meta-frontiers are re-estimated from each sample, and only sampled decision-making units are rescored. Spearman rank correlation and the change in median score are compared with the scores of the same units under the full frontier.
