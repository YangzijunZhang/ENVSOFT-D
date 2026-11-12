"""Build the 76-city input panel from the public 10 km CSV panel.

For each resource system, a grid cell contributes to a city's totals only if all
three variables of that system (input, desirable output, undesirable output) are
finite in the cell. The input and outputs of each city therefore cover the same
area. The coverage fraction reported for each variable is the share of the city's
in-basin area covered by these jointly valid grid cells.
"""

import numpy as np
import pandas as pd

from common import ROOT, SYSTEMS, YEARS, load_grid, read_csv, write_csv


def aggregate_city_panel() -> pd.DataFrame:
    grid = load_grid()
    weights = read_csv(
        ROOT / "data/city_grid_weights_10km.csv",
        dtype={"city_id": str, "province_id": str},
    )
    cities = read_csv(
        ROOT / "data/city_reference.csv",
        dtype={"city_id": str, "province_id": str},
    )
    required = {"city_id", "grid_row", "grid_col", "allocation_weight", "overlap_area_m2"}
    if not required.issubset(weights.columns):
        raise ValueError(f"City weights are missing: {sorted(required - set(weights.columns))}")

    weights = weights[
        ["city_id", "grid_row", "grid_col", "allocation_weight", "overlap_area_m2"]
    ]
    indicators = [column for system in SYSTEMS.values() for column in system.variables]
    output = []
    for year in YEARS:
        joined = weights.merge(
            grid.loc[grid["year"].eq(year), ["grid_row", "grid_col", *indicators]],
            on=["grid_row", "grid_col"],
            how="left",
            validate="many_to_one",
        )
        city_area = joined.groupby("city_id", observed=True)["overlap_area_m2"].sum()
        totals = pd.DataFrame(index=city_area.index)
        for system in SYSTEMS.values():
            variables = list(system.variables)
            valid = np.isfinite(joined[variables].to_numpy(float)).all(axis=1)
            valid_area = (
                joined["overlap_area_m2"].where(valid, 0).groupby(joined["city_id"], observed=True).sum()
            )
            for indicator in variables:
                contribution = (joined[indicator] * joined["allocation_weight"]).where(valid)
                totals[indicator] = contribution.groupby(joined["city_id"], observed=True).sum(
                    min_count=1
                )
                totals[f"{indicator}__coverage_fraction"] = valid_area / city_area
        totals["year"] = year
        output.append(totals.reset_index())

    panel = pd.concat(output, ignore_index=True).merge(
        cities, on="city_id", how="left", validate="many_to_one"
    )
    panel["scenario_dmu_id"] = panel["city_id"] + "__" + panel["year"].astype(str)
    panel = panel.sort_values(["year", "city_index"], kind="stable")
    if len(panel) != 76 * len(YEARS) or panel.duplicated(["city_id", "year"]).any():
        raise ValueError("City aggregation must produce 76 unique cities in each benchmark year")
    return panel


if __name__ == "__main__":
    write_csv(aggregate_city_panel(), ROOT / "data/city_panel_10km.csv")
