"""Shared map drawing from the repository's processed WKT-CSV layers."""

from functools import lru_cache

from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
from shapely import wkt
from shapely.geometry.polygon import orient

from common import ROOT, read_csv


REACH_STYLE = {
    "Upper": ("#C6C6C7", "#69696A"),
    "Middle": ("#E4E4E4", "#EFAFAF"),
    "Lower": ("#ECF7EB", "#C2C7C2"),
}
SYSTEM_COLORS = {"water": "#2E86C1", "energy": "#F1C40F", "food": "#27AE60"}


@lru_cache(maxsize=1)
def layers() -> pd.DataFrame:
    frame = read_csv(ROOT / "data/cartography_layers.csv", dtype=str)
    frame["geometry"] = frame["geometry_wkt"].map(wkt.loads)
    return frame


def _parts(geometry, kind: str):
    if geometry is None or geometry.is_empty:
        return
    if geometry.geom_type == kind:
        yield geometry
    elif geometry.geom_type in {f"Multi{kind}", "GeometryCollection"}:
        for part in geometry.geoms:
            yield from _parts(part, kind)


def _polygon_patch(polygon) -> PathPatch:
    polygon = orient(polygon, sign=1.0)
    vertices, codes = [], []
    for ring in (polygon.exterior, *polygon.interiors):
        xy = np.asarray(ring.coords, dtype=float)
        code = np.full(len(xy), MplPath.LINETO, dtype=np.uint8)
        code[0], code[-1] = MplPath.MOVETO, MplPath.CLOSEPOLY
        vertices.append(xy)
        codes.append(code)
    return PathPatch(MplPath(np.vstack(vertices), np.concatenate(codes)))


def add_polygons(
    axis,
    frame: pd.DataFrame,
    *,
    values=None,
    cmap=None,
    norm=None,
    facecolor="none",
    edgecolor="none",
    linewidth=0,
    zorder=1,
):
    patches, repeated = [], []
    for position, geometry in enumerate(frame["geometry"]):
        for polygon in _parts(geometry, "Polygon"):
            patches.append(_polygon_patch(polygon))
            if values is not None:
                repeated.append(values[position])
    collection = PatchCollection(
        patches,
        cmap=cmap,
        norm=norm,
        facecolor=facecolor if values is None else None,
        edgecolor=edgecolor,
        linewidth=linewidth,
        rasterized=values is not None,
        zorder=zorder,
    )
    if values is not None:
        collection.set_array(np.asarray(repeated, dtype=float))
    axis.add_collection(collection)
    return collection


def add_lines(axis, frame: pd.DataFrame, color, linewidth, zorder=5) -> None:
    segments = []
    for geometry in frame["geometry"]:
        for line in _parts(geometry, "LineString"):
            segments.append(np.asarray(line.coords, dtype=float))
    axis.add_collection(
        LineCollection(segments, colors=color, linewidths=linewidth, zorder=zorder)
    )


def basin_bounds() -> tuple[float, float, float, float]:
    return tuple(layers().loc[layers()["layer"].eq("basin"), "geometry"].iloc[0].bounds)


def set_extent(axis, pad=0.015) -> None:
    xmin, ymin, xmax, ymax = basin_bounds()
    axis.set_xlim(xmin - pad * (xmax - xmin), xmax + pad * (xmax - xmin))
    axis.set_ylim(ymin - pad * (ymax - ymin), ymax + pad * (ymax - ymin))
    axis.set_aspect("equal")
    axis.set_axis_off()


def draw_base(axis, *, reaches=False, cities=False, river=False) -> None:
    frame = layers()
    basin = frame[frame["layer"].eq("basin")]
    if reaches:
        for name, (face, edge) in REACH_STYLE.items():
            reach = frame[frame["layer"].eq("reach") & frame["name"].eq(name)]
            add_polygons(axis, reach, facecolor=face, edgecolor=edge, linewidth=0.45)
    else:
        add_polygons(axis, basin, facecolor="#ECECE9", edgecolor="none")
    if cities:
        city = frame[frame["layer"].eq("city")]
        add_polygons(axis, city, facecolor="none", edgecolor="#A9A9A9", linewidth=0.22, zorder=3)
    if river:
        add_lines(
            axis,
            frame[frame["layer"].eq("river")],
            color="#00BFFF",
            linewidth=1.1,
            zorder=5,
        )
    add_polygons(axis, basin, facecolor="none", edgecolor="#454545", linewidth=0.6, zorder=6)
    set_extent(axis)


def performance_cmap(system: str):
    return LinearSegmentedColormap.from_list(
        system, ["#FFFFFF", SYSTEM_COLORS[system]], N=256
    )


def add_grid_scores(axis, frame: pd.DataFrame, score: str, cmap, norm=None):
    norm = norm or Normalize(0.5, 1)
    valid = frame.dropna(subset=[score, "x_center_m", "y_center_m"])
    patches = [
        Rectangle((row.x_center_m - 5_000, row.y_center_m - 5_000), 10_000, 10_000)
        for row in valid.itertuples(index=False)
    ]
    collection = PatchCollection(
        patches, cmap=cmap, norm=norm, edgecolor="none", rasterized=True, zorder=2
    )
    collection.set_array(valid[score].to_numpy(float))
    basin = layers().loc[layers()["layer"].eq("basin"), "geometry"].iloc[0]
    clip_path = _polygon_patch(next(_parts(basin, "Polygon")))
    clip_path.set_transform(axis.transData)
    collection.set_clip_path(clip_path)
    axis.add_collection(collection)
    return collection


def add_city_scores(axis, scores: pd.DataFrame, score: str, cmap, norm=None):
    norm = norm or Normalize(0.5, 1)
    draw_base(axis)
    cities = layers()[layers()["layer"].eq("city")].copy()
    lookup = scores.assign(city_id=scores["city_id"].astype(str)).set_index("city_id")[score]
    cities["score"] = cities["feature_id"].map(lookup)
    missing = cities[cities["score"].isna()]
    valid = cities[cities["score"].notna()]
    add_polygons(axis, missing, facecolor="#D9D9D9", edgecolor="#444444", linewidth=0.35)
    collection = add_polygons(
        axis,
        valid,
        values=valid["score"].to_numpy(float),
        cmap=cmap,
        norm=norm,
        edgecolor="#333333",
        linewidth=0.35,
        zorder=2,
    )
    return collection
