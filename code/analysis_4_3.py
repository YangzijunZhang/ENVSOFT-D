"""Section 4.3: score-weighted SDEs for all evaluated grid and city DMUs."""

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
import numpy as np
import pandas as pd
from pyproj import Geod, Transformer
from shapely.geometry import Point

from common import ALBERS_CRS, PRIMARY_SCORE, ROOT, SYSTEMS, YEARS, read_csv, write_csv
from map_utils import draw_base, layers


GEOD = Geod(ellps="WGS84")
YEAR_COLORS = {2010: "#1F77B4", 2015: "#2CA02C", 2020: "#D62728"}


def add_center_context(parameters: pd.DataFrame) -> pd.DataFrame:
    """Locate each SDE centre in the clipped city and reach polygons."""
    mapping = layers()
    cities = mapping[mapping["layer"].eq("city")].reset_index(drop=True)
    reaches = mapping[mapping["layer"].eq("reach")].reset_index(drop=True)
    rows = []
    for row in parameters.itertuples(index=False):
        point = Point(row.center_x_m, row.center_y_m)
        city_distance = cities["geometry"].map(lambda geometry: geometry.distance(point))
        containing_city = cities["geometry"].map(lambda geometry: geometry.covers(point))
        city_index = containing_city.idxmax() if containing_city.any() else city_distance.idxmin()
        city = cities.loc[city_index]
        reach_distance = reaches["geometry"].map(lambda geometry: geometry.distance(point))
        containing_reach = reaches["geometry"].map(lambda geometry: geometry.covers(point))
        reach_index = (
            containing_reach.idxmax() if containing_reach.any() else reach_distance.idxmin()
        )
        rows.append(
            {
                "center_city_id": city["feature_id"],
                "center_city_name": city["name"],
                "center_city_province": city["parent_name"],
                "center_city_distance_km": city_distance.loc[city_index] / 1000,
                "center_reach": reaches.loc[reach_index, "name"],
                "center_reach_distance_km": reach_distance.loc[reach_index] / 1000,
            }
        )
    return pd.concat(
        [parameters.reset_index(drop=True), pd.DataFrame(rows)], axis=1
    )


def sde(frame: pd.DataFrame, score: str) -> dict[str, float | int]:
    values = frame[["x_center_m", "y_center_m", score]].apply(pd.to_numeric, errors="coerce")
    values = values[np.isfinite(values).all(axis=1) & values[score].gt(0)]
    weights = values[score].to_numpy(float)
    coordinates = values[["x_center_m", "y_center_m"]].to_numpy(float)
    centre = np.average(coordinates, axis=0, weights=weights)
    centred = coordinates - centre
    covariance = (centred * weights[:, None]).T @ centred / weights.sum()
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 0)
    major_vector = eigenvectors[:, order[0]]
    major, minor = np.sqrt(eigenvalues)
    return {
        "n": len(values),
        "weight_sum": weights.sum(),
        "center_x_m": centre[0],
        "center_y_m": centre[1],
        "major_semi_axis_m": major,
        "minor_semi_axis_m": minor,
        "orientation_deg_ccw_from_east": math.degrees(
            math.atan2(major_vector[1], major_vector[0])
        )
        % 180,
        "ellipse_area_km2": math.pi * major * minor / 1_000_000,
    }


def sde_tables(
    frames: dict[str, pd.DataFrame], score_columns: dict[str, str], scenario: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for system, frame in frames.items():
        for year in YEARS:
            rows.append(
                {
                    "scenario": scenario,
                    "system": system,
                    "year": year,
                    "score_column": score_columns[system],
                    **sde(frame[frame["year"].eq(year)], score_columns[system]),
                }
            )
    parameters = pd.DataFrame(rows)
    transformer = Transformer.from_crs(ALBERS_CRS, "EPSG:4326", always_xy=True)
    lon, lat = transformer.transform(parameters["center_x_m"], parameters["center_y_m"])
    parameters["center_longitude_wgs84"] = lon
    parameters["center_latitude_wgs84"] = lat
    parameters = add_center_context(parameters)

    directions = ("N", "NE", "E", "SE", "S", "SW", "W", "NW")
    moves = []
    for system, group in parameters.groupby("system", observed=True):
        indexed = group.set_index("year")
        start, end = indexed.loc[2010], indexed.loc[2020]
        dx = end["center_x_m"] - start["center_x_m"]
        dy = end["center_y_m"] - start["center_y_m"]
        azimuth, _, distance = GEOD.inv(
            start["center_longitude_wgs84"],
            start["center_latitude_wgs84"],
            end["center_longitude_wgs84"],
            end["center_latitude_wgs84"],
        )
        azimuth %= 360
        moves.append(
            {
                "scenario": scenario,
                "system": system,
                "from_year": 2010,
                "to_year": 2020,
                "dx_east_km": dx / 1000,
                "dy_north_km": dy / 1000,
                "geodesic_displacement_km": distance / 1000,
                "azimuth_deg_clockwise_from_north": azimuth,
                "movement_direction_8point": directions[int((azimuth + 22.5) // 45) % 8],
                "from_center_city": start["center_city_name"],
                "to_center_city": end["center_city_name"],
                "from_reach": start["center_reach"],
                "to_reach": end["center_reach"],
            }
        )
    return parameters, pd.DataFrame(moves)


def _ellipse_panel(axis, parameters: pd.DataFrame, title: str) -> None:
    draw_base(axis, reaches=True, cities=True, river=True)
    parameters = parameters.sort_values("year")
    for row in parameters.itertuples(index=False):
        axis.add_patch(
            Ellipse(
                (row.center_x_m, row.center_y_m),
                2 * row.major_semi_axis_m,
                2 * row.minor_semi_axis_m,
                angle=row.orientation_deg_ccw_from_east,
                facecolor=YEAR_COLORS[row.year],
                edgecolor=YEAR_COLORS[row.year],
                alpha=0.16,
                linewidth=1.7,
                label=str(row.year),
                zorder=7,
            )
        )
        axis.add_patch(
            Ellipse(
                (row.center_x_m, row.center_y_m),
                2 * row.major_semi_axis_m,
                2 * row.minor_semi_axis_m,
                angle=row.orientation_deg_ccw_from_east,
                fill=False,
                edgecolor=YEAR_COLORS[row.year],
                linewidth=1.7,
                zorder=8,
            )
        )
    axis.plot(
        parameters["center_x_m"],
        parameters["center_y_m"],
        "x-",
        color="#111111",
        linewidth=0.8,
        markersize=5,
        markeredgewidth=1.2,
        zorder=9,
    )
    axis.set_title(title, fontsize=10, fontweight="bold", loc="left")


def sde_legend_handles():
    return [
        Line2D([0], [0], color=color, linewidth=1.8, label=str(year))
        for year, color in YEAR_COLORS.items()
    ] + [
        Line2D([0], [0], marker="x", color="#111111", linestyle="none", label="Ellipse centre"),
        Line2D([0], [0], color="#00BFFF", linewidth=1.2, label="Yellow River"),
        Line2D([0], [0], color="#A9A9A9", linewidth=0.8, label="City boundary"),
        Patch(facecolor="#C6C6C7", edgecolor="#69696A", label="Upper reach"),
        Patch(facecolor="#E4E4E4", edgecolor="#EFAFAF", label="Middle reach"),
        Patch(facecolor="#ECF7EB", edgecolor="#C2C7C2", label="Lower reach"),
    ]


def run_section_4_3() -> None:
    grid = {
        system: read_csv(ROOT / f"results/performance/grid_{system}.csv")
        for system in SYSTEMS
    }
    city = {
        system: read_csv(ROOT / f"results/performance/city_{system}.csv")
        for system in SYSTEMS
    }
    grid_parameters, grid_moves = sde_tables(grid, PRIMARY_SCORE, "dominant_grid")
    city_parameters, city_moves = sde_tables(
        city, {system: "city_score" for system in SYSTEMS}, "city"
    )
    output = ROOT / "results/analysis_4_3"
    write_csv(grid_parameters, output / "grid_sde_parameters.csv")
    write_csv(grid_moves, output / "grid_sde_displacement_2010_2020.csv")
    write_csv(city_parameters, output / "city_sde_parameters.csv")
    write_csv(city_moves, output / "city_sde_displacement_2010_2020.csv")

    fig = plt.figure(figsize=(11, 10), layout="constrained")
    layout = fig.add_gridspec(4, 2, height_ratios=(1, 1, 1, 0.20))
    panel = 0
    for row, system in enumerate(SYSTEMS):
        _ellipse_panel(
            fig.add_subplot(layout[row, 0]),
            grid_parameters[grid_parameters["system"].eq(system)],
            f"({chr(97 + panel)}) {system.capitalize()}: grid scale",
        )
        panel += 1
        _ellipse_panel(
            fig.add_subplot(layout[row, 1]),
            city_parameters[city_parameters["system"].eq(system)],
            f"({chr(97 + panel)}) {system.capitalize()}: city scale",
        )
        panel += 1
    summary = fig.add_subplot(layout[3, :])
    summary.axis("off")
    text = ["(g) 2010–2020 centre displacement"]
    for scale, frame in (("Grid", grid_moves), ("City", city_moves)):
        scale_text = []
        for row in frame.itertuples(index=False):
            scale_text.append(
                f"{row.system.capitalize()}: {row.geodesic_displacement_km:.3f} km "
                f"{row.movement_direction_8point}; "
                f"{row.from_center_city} to {row.to_center_city}"
            )
        text.append(f"{scale} — " + " | ".join(scale_text))
    summary.text(0, 1, "\n".join(text), va="top", fontsize=8.5)
    summary.legend(
        handles=sde_legend_handles(),
        loc="lower center",
        ncol=5,
        frameon=False,
    )
    fig.suptitle("Score-weighted standard deviational ellipses", fontweight="bold")
    figure_dir = ROOT / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig.savefig(figure_dir / "Fig_7.png", dpi=300, bbox_inches="tight")
    fig.savefig(figure_dir / "Fig_7.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    run_section_4_3()
