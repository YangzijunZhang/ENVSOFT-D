"""Appendix D: sensitivity to the land-use assignment rule and robustness of the frontiers.

Reproduces every table, figure and number reported in Appendix D of the
Supplementary Information from the public CSV panel. No raster file is needed.

D.1  Design                      results/appendix_d/method_decisions.csv
                                 results/appendix_d/run_config.json
D.2  Table D1                    results/appendix_d/Table_D1.csv, Table_D1_long.csv
D.3  Table D2 and Figs. D1-D6    results/appendix_d/Table_D2.csv
                                 results/appendix_d/grid_rank_correlations.csv
                                 results/appendix_d/regional_medians_2020.csv
                                 results/appendix_d/regional_order_retention.csv
                                 results/appendix_d/sde_parameters.csv
                                 results/appendix_d/sde_displacements.csv
                                 results/appendix_d/sde_vs_dominant.csv
                                 results/appendix_d/scores/<rule>_<system>.csv
                                 figures/Fig_D1 ... Fig_D6
D.4  Fig. D7                     results/appendix_d/subsample_replicate_metrics.csv
                                 results/appendix_d/subsample_summary.csv
                                 figures/Fig_D7
Text values                      results/appendix_d/appendix_d_summary.md

Specification (identical to the main analysis): observed-value direction
g = (x_o, y_o, b_o), variable returns to scale, weak disposability of the
undesirable output, one frontier pooled over 2010, 2015 and 2020, TE = 1/(1+beta).
TE^M is reported for water and TE^g for energy and food. City totals do not use
a land-use assignment rule, so the city-scale results are not re-estimated.

Usage (from the repository root)
    python code/appendix_d.py                         # all stages
    python code/appendix_d.py --stages thresholds     # D.2-D.3, about 4 minutes
    python code/appendix_d.py --stages subsample --workers 8
    python code/appendix_d.py --stages summary        # rebuild the text values only

The subsampling stage solves the water frontiers 100 times and takes about
2.5 hours on one core; --workers runs replicates in parallel and gives
identical results, because every (system, replicate) task keeps its own seed.
"""

import argparse
import hashlib
import json
import math
import platform
from concurrent.futures import ProcessPoolExecutor
from importlib import metadata

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from scipy.stats import spearmanr

from analysis_4_1 import plot_grid_city_figure
from analysis_4_2 import plot_section_4_2, region_code
from analysis_4_3 import YEAR_COLORS, sde_tables
from common import (
    LAND_TYPES,
    PRIMARY_SCORE,
    REGIONS,
    ROOT,
    SYSTEMS,
    YEARS,
    common_city_score_sample,
    describe,
    grid_candidates,
    load_grid,
    read_csv,
    write_csv,
)
from dea_performance import score_grid_frontiers
from map_utils import draw_base


OUT = ROOT / "results/appendix_d"
ALTERNATIVE_RULES = ("threshold_50", "threshold_70")
RULE_LABEL = {"dominant": "Dominant-class", "threshold_50": "≥50%", "threshold_70": "≥70%"}
FIGURE_NUMBER = {"threshold_50": ("D1", "D3", "D5"), "threshold_70": ("D2", "D4", "D6")}
MINIMUM_DMUS = 9
SEED = 2602569
REPLICATES = 100
SAMPLE_FRACTION = 0.80
POINTS8 = ("N", "NE", "E", "SE", "S", "SW", "W", "NW")
POINTS16 = "N NNE NE ENE E ESE SE SSE S SSW SW WSW W WNW NW NNW".split()
DIRECTION_WORDS = {"N": "North", "NE": "North-east", "E": "East", "SE": "South-east",
                   "S": "South", "SW": "South-west", "W": "West", "NW": "North-west"}
SHOW_AZIMUTH = True  # report "North-east (53.2°)" in panel (g); set False for words only
SERIF = ["Times New Roman", "Liberation Serif", "TeX Gyre Termes", "DejaVu Serif"]
REACH_LEGEND = (("Upper reaches", "#C6C6C7", "#69696A"), ("Middle reaches", "#E4E4E4", "#EFAFAF"),
                ("Lower reaches", "#ECF7EB", "#C2C7C2"))
FRONTIERS = (
    ("water", "group", "Cropland", "Water–Cropland"),
    ("water", "group", "Forest", "Water–Forest"),
    ("water", "group", "Grassland", "Water–Grassland"),
    ("water", "meta", "Cropland,Forest,Grassland", "Water meta-frontier"),
    ("energy", "group", "Impervious", "Energy–Impervious"),
    ("food", "group", "Cropland", "Food–Cropland"),
)


# ----------------------------------------------------------------------------- helpers
def _spearman(first: pd.Series, second: pd.Series) -> float:
    values = pd.concat([first, second], axis=1).dropna()
    if len(values) < 2 or values.iloc[:, 0].nunique() < 2 or values.iloc[:, 1].nunique() < 2:
        return np.nan
    return float(spearmanr(values.iloc[:, 0], values.iloc[:, 1]).statistic)


def _save(fig, name: str) -> None:
    output = ROOT / "figures"
    output.mkdir(exist_ok=True)
    fig.savefig(output / f"Fig_{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output / f"Fig_{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def _direction(azimuth: float, points: tuple | list) -> str:
    width = 360 / len(points)
    return points[int((azimuth + width / 2) // width) % len(points)]


def _angle_difference(first: float, second: float) -> float:
    return abs((first - second + 180.0) % 360.0 - 180.0)


def _order(medians: pd.DataFrame) -> str:
    ranked = medians.sort_values(["median", "region_code"], ascending=[False, True])
    return " > ".join(ranked["region_code"])


def _sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _baseline_scores() -> dict[str, pd.DataFrame]:
    return {name: read_csv(ROOT / f"results/performance/grid_{name}.csv") for name in SYSTEMS}


# ----------------------------------------------------------------------------- D.1
def write_design() -> None:
    decisions = [
        ("land_use_assignment", "Dominant-class rule (primary), dominant share >= 50% and "
         "dominant share >= 70%; cells failing a threshold are excluded, not reclassified.", True),
        ("grid_size", "10 km in every run; the thresholds test within-cell land-use mixing directly.", False),
        ("direction_vector", "Observed values g = (x_o, y_o, b_o) in every run; beta is invariant to "
         "positive changes in measurement units, so no scaling constant is needed.", False),
        ("frontier", "One frontier pooled over 2010, 2015 and 2020, VRS, weak disposability of the "
         "undesirable output, as in the main analysis.", True),
        ("frontier_stability", f"{REPLICATES} replicates retaining {SAMPLE_FRACTION:.0%} of the DMUs "
         "(rounded down) within every year x land-use stratum; seed 2602569.", True),
        ("city_scale", "Not re-estimated: city totals include grid cells of all land-use types and do "
         "not depend on the land-use assignment rule.", False),
    ]
    write_csv(pd.DataFrame(decisions, columns=["topic", "decision", "requires_new_estimation"]),
              OUT / "method_decisions.csv")


def write_run_config(args: argparse.Namespace) -> None:
    files = ["data/grid_panel_10km.csv", "data/city_panel_10km.csv", "code/dea_performance.py",
             "code/common.py", "code/appendix_d.py"]
    packages = {}
    for name in ("numpy", "pandas", "scipy", "matplotlib", "pyproj", "shapely"):
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            packages[name] = None
    config = {
        "rules": ["dominant", *ALTERNATIVE_RULES],
        "years": list(YEARS),
        "minimum_dmus_reported_per_year": args.minimum_dmus,
        "subsampling": {"replicates": args.replicates, "sample_fraction": args.sample_fraction,
                        "strata": ["year", "lulc_code"], "seed": args.seed,
                        "scores_compared_for": "drawn DMUs only"},
        "python": platform.python_version(),
        "packages": packages,
        "sha256": {name: _sha256(ROOT / name) for name in files if (ROOT / name).exists()},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "run_config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")


# ----------------------------------------------------------------------------- D.2-D.3
def _d1_rows(frame: pd.DataFrame, rule: str, system: str, minimum: int) -> list[dict]:
    score = PRIMARY_SCORE[system]
    rows = []
    for (year, land_type), group in frame.groupby(["year", "land_type"], observed=True):
        summary = describe(group[score])
        rows.append({"lulc_rule": rule, "system": system, "land_type": land_type, "year": year,
                     **summary, "reported": summary["n"] >= minimum})
    return rows


def _table_d1(rows: list[dict]) -> None:
    long = pd.DataFrame(rows)
    long.loc[~long["reported"], "median"] = np.nan
    write_csv(long, OUT / "Table_D1_long.csv")
    long["table_value"] = long.apply(
        lambda row: f"{int(row['n']):,} ({row['median']:.3f})" if row["reported"] else "—", axis=1
    )
    wide = (long.pivot(index=["system", "land_type", "year"], columns="lulc_rule", values="table_value")
            .reset_index().rename_axis(columns=None))
    write_csv(wide[["system", "land_type", "year", "dominant", *ALTERNATIVE_RULES]], OUT / "Table_D1.csv")


def _regional(scenario: dict[str, pd.DataFrame], rule: str) -> list[dict]:
    rows = []
    for name in ("water", "food"):
        cropland = scenario[name]
        cropland = cropland[cropland["year"].eq(2020) & cropland["lulc_code"].eq(1)].copy()
        cropland["region_code"] = region_code(cropland)
        for label, region_name, *_ in REGIONS:
            rows.append({"lulc_rule": rule, "system": name, "year": 2020, "region_code": label,
                         "region_name": region_name,
                         **describe(cropland.loc[cropland["region_code"].eq(label), "group_score"])})
    return rows


def _sde_map(axis, parameters: pd.DataFrame, title: str) -> None:
    """One map panel: ellipses of 2010, 2015 and 2020 and the centre city."""
    draw_base(axis, reaches=True, cities=True, river=True)
    parameters = parameters.sort_values("year")
    for row in parameters.itertuples(index=False):
        colour = YEAR_COLORS[row.year]
        shape = dict(xy=(row.center_x_m, row.center_y_m), width=2 * row.major_semi_axis_m,
                     height=2 * row.minor_semi_axis_m, angle=row.orientation_deg_ccw_from_east)
        axis.add_patch(Ellipse(**shape, facecolor=colour, edgecolor="none", alpha=0.16, zorder=7))
        axis.add_patch(Ellipse(**shape, fill=False, edgecolor=colour, linewidth=1.7, zorder=8))
    axis.plot(parameters["center_x_m"], parameters["center_y_m"], "x-", color="#111111",
              linewidth=0.8, markersize=5, markeredgewidth=1.2, zorder=9)
    end = parameters.iloc[-1]  # label the 2020 centre; both centre cities are listed in panel (g)
    label = axis.annotate(end["center_city_name"], (end["center_x_m"], end["center_y_m"]), xytext=(4, 0),
                          textcoords="offset points", ha="left", va="center", fontsize=12, zorder=10)
    label.set_path_effects([patheffects.withStroke(linewidth=2.5, foreground="white")])
    label.set_in_layout(False)
    axis.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=2)


def _direction_label(move) -> str:
    word = DIRECTION_WORDS[move.movement_direction_8point]
    return f"{word} ({move.azimuth_deg_clockwise_from_north:.1f}°)" if SHOW_AZIMUTH else word


def _displacement_table(axis, grid_moves: pd.DataFrame, city_moves: pd.DataFrame) -> None:
    """Panel (g): three-line table of the 2010-2020 centre displacement."""
    axis.axis("off")
    columns = (("Scale", 0.10), ("System", 0.255), ("2010 Centre", 0.40), ("2020 Centre", 0.55),
               ("Movement direction", 0.735), ("Distance, km", 0.925))
    rows = []
    for scale, moves in (("Land-use-based", grid_moves), ("City-scale", city_moves)):
        moves = moves.set_index("system")
        for name in SYSTEMS:
            move = moves.loc[name]
            rows.append((scale, name.capitalize(), move.from_center_city, move.to_center_city,
                         _direction_label(move), f"{move.geodesic_displacement_km:.3f}"))
    top, header, first, step = 0.90, 0.80, 0.66, 0.115
    axis.text(0, 1.04, "(g) SDE centre displacement from 2010 to 2020", transform=axis.transAxes,
              fontsize=14, fontweight="bold", va="bottom")
    for y, width in ((top, 1.6), (header - 0.07, 0.8), (first - step * (len(rows) - 0.5), 1.6)):
        axis.plot([0, 1], [y, y], color="black", linewidth=width, transform=axis.transAxes, clip_on=False)
    for label, x in columns:
        axis.text(x, header, label, transform=axis.transAxes, ha="center", va="center",
                  fontsize=11.5, fontweight="bold")
    for index, row in enumerate(rows):
        for (_, x), value in zip(columns, row):
            axis.text(x, first - index * step, value, transform=axis.transAxes, ha="center",
                      va="center", fontsize=11.5)


def plot_sde_figure(grid_parameters: pd.DataFrame, grid_moves: pd.DataFrame,
                    city_parameters: pd.DataFrame, city_moves: pd.DataFrame, name: str) -> None:
    """Figs. D5-D6 (same layout as Fig. 7): grid and city SDEs, legend and panel (g)."""
    with plt.rc_context({"font.family": "serif", "font.serif": SERIF, "font.size": 12}):
        fig = plt.figure(figsize=(8.0, 10.5), layout="constrained")
        layout = fig.add_gridspec(5, 2, height_ratios=(1, 1, 1, 0.34, 0.70))
        panel = 0
        for row, system in enumerate(SYSTEMS):
            for column, (scale, parameters) in enumerate((("Grid", grid_parameters), ("City", city_parameters))):
                _sde_map(fig.add_subplot(layout[row, column]), parameters[parameters["system"].eq(system)],
                         f"({chr(97 + panel)}) {system.capitalize()} - {scale}")
                panel += 1
        legend = fig.add_subplot(layout[3, :])
        legend.axis("off")
        handles = [Patch(facecolor=matplotlib.colors.to_rgba(colour, 0.16), edgecolor=colour,
                         linewidth=1.7, label=str(year)) for year, colour in YEAR_COLORS.items()]
        handles += [Line2D([0], [0], marker="x", color="#111111", linestyle="none", markersize=8,
                           markeredgewidth=1.6, label="Ellipse centre"),
                    Line2D([0], [0], color="#00BFFF", linewidth=1.6, label="The Yellow River")]
        reaches = [Patch(facecolor=face, edgecolor=edge, label=label) for label, face, edge in REACH_LEGEND]
        first = legend.legend(handles=handles, loc="upper center", ncol=5, frameon=False, fontsize=12,
                              handlelength=1.8, columnspacing=1.6, bbox_to_anchor=(0.5, 1.0))
        legend.add_artist(first)
        second = legend.legend(handles=reaches, loc="lower center", ncol=3, frameon=False, fontsize=12,
                               handlelength=1.8, columnspacing=3.0, bbox_to_anchor=(0.5, 0.0))
        for item in (first, second):
            item.set_in_layout(False)
        _displacement_table(fig.add_subplot(layout[4, :]), grid_moves, city_moves)
        _save(fig, name)


def _figures(rule: str, scenario: dict[str, pd.DataFrame], city, parameters: pd.DataFrame,
             moves: pd.DataFrame, city_sde: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    scale, regions, ellipses = FIGURE_NUMBER[rule]
    plot_grid_city_figure(scenario, city, f"Fig_{scale}",
                          f"Resource performance under the dominant land-use share {RULE_LABEL[rule]} rule")
    systems = {}
    for name in ("water", "food"):
        frame = scenario[name]
        frame = frame[frame["year"].eq(2020) & frame["lulc_code"].eq(1)].copy()
        frame["region_code"] = region_code(frame)
        systems[name] = frame
    plot_section_4_2(systems, f"Fig_{regions}")
    plot_sde_figure(parameters, moves, *city_sde, ellipses)


def _table_d2(correlations: pd.DataFrame, regional: pd.DataFrame, moves: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rule in ("dominant", *ALTERNATIVE_RULES):
        row = {"lulc_rule": RULE_LABEL[rule]}
        if rule == "dominant":
            row["section_4_1_grid_rank_correlation"] = "Reference"
        else:
            rho = correlations.loc[correlations["lulc_rule"].eq(rule), "spearman_rho"]
            row["section_4_1_grid_rank_correlation"] = f"{rho.min():.3f}–{rho.max():.3f}"
        for name in ("water", "food"):
            subset = regional[regional["lulc_rule"].eq(rule) & regional["system"].eq(name)]
            row[f"section_4_2_{name}_cropland_order"] = _order(subset)
        for name in SYSTEMS:
            move = moves[moves["scenario"].eq(rule) & moves["system"].eq(name)].iloc[0]
            row[f"section_4_3_{name}_shift"] = (
                f"{move.geodesic_displacement_km:.3f} km {move.movement_direction_8point} "
                f"({move.azimuth_deg_clockwise_from_north:.1f}°)"
            )
        rows.append(row)
    return pd.DataFrame(rows)


def run_thresholds(minimum: int) -> None:
    grid = load_grid()
    baseline = _baseline_scores()
    city = {name: read_csv(ROOT / f"results/performance/city_{name}.csv", dtype={"city_id": str})
            for name in SYSTEMS}
    d1 = [row for name, frame in baseline.items() for row in _d1_rows(frame, "dominant", name, minimum)]
    regional = _regional(baseline, "dominant")
    _, dominant_moves = sde_tables(baseline, PRIMARY_SCORE, "dominant")
    # City totals do not depend on the land-use rule, so the city panels of Figs. D5-D6
    # are those of Fig. 7, based on the same 76 cities in every year (Appendix A.4).
    city_common = {name: common_city_score_sample(frame, "city_score") for name, frame in city.items()}
    city_sde = sde_tables(city_common, {name: "city_score" for name in SYSTEMS}, "city")
    write_csv(city_sde[0], OUT / "city_sde_parameters.csv")
    write_csv(city_sde[1], OUT / "city_sde_displacements.csv")
    correlations, all_parameters, all_moves, references = [], [], [dominant_moves], []

    for rule in ALTERNATIVE_RULES:
        scenario = {}
        for name, system in SYSTEMS.items():
            scores, reference = score_grid_frontiers(grid_candidates(grid, rule, name), system)
            scenario[name] = scores
            reference.insert(0, "lulc_rule", rule)
            references.append(reference)
            write_csv(scores, OUT / f"scores/{rule}_{name}.csv")
            d1.extend(_d1_rows(scores, rule, name, minimum))

            score = PRIMARY_SCORE[name]
            primary = baseline[name][["dmu_id", "year", "land_type", score]].rename(columns={score: "dominant_score"})
            alternative = scores[["dmu_id", score]].rename(columns={score: "scenario_score"})
            common = primary.merge(alternative, on="dmu_id", how="inner", validate="one_to_one")
            for (year, land_type), group in common.groupby(["year", "land_type"], observed=True):
                correlations.append({
                    "lulc_rule": rule, "system": name, "year": year, "land_type": land_type,
                    "n_common_dmus": len(group),
                    "spearman_rho": _spearman(group["dominant_score"], group["scenario_score"]),
                    "median_absolute_score_difference":
                        (group["scenario_score"] - group["dominant_score"]).abs().median(),
                })
        regional.extend(_regional(scenario, rule))
        parameters, moves = sde_tables(scenario, PRIMARY_SCORE, rule)
        all_parameters.append(parameters)
        all_moves.append(moves)
        _figures(rule, scenario, city, parameters, moves, city_sde)

    _table_d1(d1)
    correlations, regional = pd.DataFrame(correlations), pd.DataFrame(regional)
    write_csv(pd.concat(references, ignore_index=True), OUT / "frontier_references.csv")
    write_csv(correlations, OUT / "grid_rank_correlations.csv")
    write_csv(regional, OUT / "regional_medians_2020.csv")

    orders = []
    for name in ("water", "food"):
        reference = _order(regional[regional["lulc_rule"].eq("dominant") & regional["system"].eq(name)])
        for rule in ALTERNATIVE_RULES:
            order = _order(regional[regional["lulc_rule"].eq(rule) & regional["system"].eq(name)])
            orders.append({"system": name, "lulc_rule": rule, "dominant_order": reference,
                           "scenario_order": order, "ordering_retained": order == reference})
    write_csv(pd.DataFrame(orders), OUT / "regional_order_retention.csv")

    moves = pd.concat(all_moves, ignore_index=True)
    moves["movement_direction_16point"] = moves["azimuth_deg_clockwise_from_north"].map(
        lambda value: _direction(value, POINTS16))
    write_csv(pd.concat(all_parameters, ignore_index=True), OUT / "sde_parameters.csv")
    write_csv(moves, OUT / "sde_displacements.csv")
    base = moves[moves["scenario"].eq("dominant")].set_index("system")
    comparison = []
    for row in moves[moves["scenario"].isin(ALTERNATIVE_RULES)].itertuples(index=False):
        reference = base.loc[row.system]
        comparison.append({
            "system": row.system, "lulc_rule": row.scenario,
            "dominant_km": reference["geodesic_displacement_km"], "scenario_km": row.geodesic_displacement_km,
            "dominant_azimuth_deg": reference["azimuth_deg_clockwise_from_north"],
            "scenario_azimuth_deg": row.azimuth_deg_clockwise_from_north,
            "azimuth_difference_deg": _angle_difference(row.azimuth_deg_clockwise_from_north,
                                                        reference["azimuth_deg_clockwise_from_north"]),
            "dominant_direction_8point": reference["movement_direction_8point"],
            "scenario_direction_8point": row.movement_direction_8point,
            "direction_retained_8point": row.movement_direction_8point == reference["movement_direction_8point"],
        })
    write_csv(pd.DataFrame(comparison), OUT / "sde_vs_dominant.csv")
    write_csv(_table_d2(correlations, regional, moves), OUT / "Table_D2.csv")


# ----------------------------------------------------------------------------- D.4
_CACHE: dict = {}


def _load_subsample_inputs() -> None:
    grid = load_grid()
    _CACHE["candidates"] = {name: grid_candidates(grid, "dominant", name).query("analysis_eligible").copy()
                            for name in SYSTEMS}
    _CACHE["full"] = _baseline_scores()


def _sample(frame: pd.DataFrame, rng: np.random.Generator, fraction: float) -> pd.DataFrame:
    indices = []
    for _, group in frame.groupby(["year", "lulc_code"], observed=True, sort=True):
        n = math.floor(len(group) * fraction)
        indices.extend(rng.choice(group.index.to_numpy(), n, replace=False))
    return frame.loc[indices].copy()


def _replicate(task: tuple[str, int, np.random.SeedSequence, float]) -> list[dict]:
    name, replicate, seed, fraction = task
    if not _CACHE:
        _load_subsample_inputs()
    system = SYSTEMS[name]
    sampled = _sample(_CACHE["candidates"][name], np.random.default_rng(seed), fraction)
    scores, _ = score_grid_frontiers(sampled, system)
    full = _CACHE["full"][name]
    comparisons = [("group", code, LAND_TYPES[code], "group_score") for code in system.land_codes]
    if system.meta:
        comparisons.append(("meta", None, "Cropland,Forest,Grassland", "meta_score"))
    rows = []
    for frontier, code, label, column in comparisons:
        subsample = scores if code is None else scores[scores["lulc_code"].eq(code)]
        complete = full if code is None else full[full["lulc_code"].eq(code)]
        common = complete[["scenario_dmu_id", column]].merge(
            subsample[["scenario_dmu_id", column]], on="scenario_dmu_id", how="inner",
            suffixes=("_full", "_subsample"), validate="one_to_one").dropna()
        change = common[f"{column}_subsample"] - common[f"{column}_full"]
        rows.append({
            "replicate": replicate, "system": name, "frontier": frontier, "land_type_or_groups": label,
            "n_sampled_dmus": len(common),
            "spearman_rho": _spearman(common[f"{column}_full"], common[f"{column}_subsample"]),
            "full_median_on_sampled_dmus": common[f"{column}_full"].median(),
            "subsample_median": common[f"{column}_subsample"].median(),
            "median_score_change": common[f"{column}_subsample"].median() - common[f"{column}_full"].median(),
            "minimum_individual_score_change": change.min(),
        })
    return rows


def _figure_d7(metrics: pd.DataFrame) -> None:
    groups, labels = [], []
    for system, frontier, land_type, label in FRONTIERS:
        subset = metrics[metrics["system"].eq(system) & metrics["frontier"].eq(frontier)
                         & metrics["land_type_or_groups"].eq(land_type)]
        if len(subset):
            groups.append(subset)
            labels.append(label)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for axis, column, title in ((axes[0], "spearman_rho", "(a) Rank correlation"),
                                (axes[1], "median_score_change", "(b) Change in median TE")):
        boxes = axis.boxplot([group[column].dropna().to_numpy() for group in groups],
                             tick_labels=labels, showfliers=True, patch_artist=True)
        for box in boxes["boxes"]:
            box.set_facecolor("#78A6D1")
            box.set_alpha(0.75)
        axis.set_title(title)
        axis.tick_params(axis="x", labelrotation=45, labelsize=7)
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    axes[0].set_ylabel("Spearman rank correlation")
    axes[1].set_ylabel("Subsample minus full-sample median TE")
    axes[1].axhline(0, color="#555555", linewidth=0.8)
    _save(fig, "D7")


def run_subsample(replicates: int, fraction: float, seed: int, workers: int, systems: list[str]) -> None:
    # The task list and seeds are always built for all systems, so a run restricted
    # with --systems reproduces exactly the corresponding rows of the full run.
    tasks = [(name, replicate) for name in SYSTEMS for replicate in range(1, replicates + 1)]
    seeds = np.random.SeedSequence(seed).spawn(len(tasks))
    work = [(name, replicate, child, fraction) for (name, replicate), child in zip(tasks, seeds)
            if name in systems]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers, initializer=_load_subsample_inputs) as pool:
            results = list(pool.map(_replicate, work, chunksize=1))
    else:
        results = [_replicate(task) for task in work]
    metrics = pd.DataFrame([row for rows in results for row in rows])
    write_csv(metrics, OUT / "subsample_replicate_metrics.csv")
    summary = []
    keys = ["system", "frontier", "land_type_or_groups"]
    for key, group in metrics.groupby(keys, observed=True):
        row = dict(zip(keys, key))
        row["n_replicates"] = len(group)
        row["n_sampled_dmus"] = int(group["n_sampled_dmus"].median())
        for column in ("spearman_rho", "median_score_change"):
            values = group[column].dropna()
            row.update({f"{column}_minimum": values.min(), f"{column}_p05": values.quantile(0.05),
                        f"{column}_median": values.median(), f"{column}_p95": values.quantile(0.95),
                        f"{column}_maximum": values.max()})
        row["minimum_individual_score_change"] = group["minimum_individual_score_change"].min()
        summary.append(row)
    write_csv(pd.DataFrame(summary), OUT / "subsample_summary.csv")
    _figure_d7(metrics)


# ----------------------------------------------------------------------------- text values
def write_summary() -> None:
    lines = ["# Appendix D: values quoted in the text", ""]
    if (OUT / "Table_D1.csv").exists():
        d1 = read_csv(OUT / "Table_D1_long.csv")
        hidden = d1[~d1["reported"].astype(bool)]
        lines += ["## D.2 Number of DMUs and score levels", "",
                  "Not reported (fewer than nine DMUs in the year): "
                  + ("; ".join(f"{r.system} {r.land_type} {r.lulc_rule} {r.year} (n = {r.n})"
                               for r in hidden.itertuples()) or "none"), ""]
        rho = read_csv(OUT / "grid_rank_correlations.csv")
        lines += ["## D.3 Retention of the main spatial results", "", "### (1) Grid-scale rank correlations"]
        for rule, group in rho.groupby("lulc_rule"):
            low = group.loc[group["spearman_rho"].idxmin()]
            minima = {name: group.loc[group["system"].eq(name), "spearman_rho"].min() for name in SYSTEMS}
            lines.append(
                f"- {RULE_LABEL[rule]}: rho = {group['spearman_rho'].min():.3f}–{group['spearman_rho'].max():.3f}; "
                f"lowest {low.system} {low.land_type} {int(low.year)} (rho = {low.spearman_rho:.3f}, "
                f"n = {int(low.n_common_dmus)}); minimum by system: "
                + ", ".join(f"{name} {value:.3f}" for name, value in minima.items()))
        regional = read_csv(OUT / "regional_medians_2020.csv")
        lines += ["", "### (2) Regional ordering of 2020 cropland medians"]
        for (rule, name), group in regional.groupby(["lulc_rule", "system"]):
            ranked = group.sort_values("median", ascending=False)
            lines.append(f"- {RULE_LABEL[rule]}, {name}: "
                         + " > ".join(f"{r.region_code} ({r.median:.3f})" for r in ranked.itertuples()))
        retained = read_csv(OUT / "regional_order_retention.csv")
        lines.append(f"- Ordering retained in every case: {bool(retained['ordering_retained'].all())}")
        moves = read_csv(OUT / "sde_displacements.csv")
        lines += ["", "### (3) Grid-scale SDE centre shifts, 2010–2020"]
        for row in moves.itertuples():
            lines.append(f"- {RULE_LABEL[row.scenario]}, {row.system}: {row.geodesic_displacement_km:.3f} km, "
                         f"azimuth {row.azimuth_deg_clockwise_from_north:.1f}° "
                         f"({row.movement_direction_8point}; 16-point {row.movement_direction_16point})")
        counts = d1[d1["system"].eq("energy") & d1["lulc_rule"].eq("threshold_70")].sort_values("year")
        lines.append("- Energy DMUs under the ≥70% rule by year: "
                     + ", ".join(f"{int(r.year)}: {int(r.n)}" for r in counts.itertuples()))
        lines.append("")
    if (OUT / "subsample_summary.csv").exists():
        summary = read_csv(OUT / "subsample_summary.csv")
        metrics = read_csv(OUT / "subsample_replicate_metrics.csv")
        label = {(s, f, t): text for s, f, t, text in FRONTIERS}
        summary["label"] = [label[(r.system, r.frontier, r.land_type_or_groups)] for r in summary.itertuples()]
        low, high = summary.loc[summary["spearman_rho_median"].idxmin()], summary.loc[summary["spearman_rho_median"].idxmax()]
        lines += ["## D.4 Robustness to subsampled frontiers", "",
                  f"- Replicates per frontier: {int(summary['n_replicates'].min())}",
                  f"- Median rho: {low.spearman_rho_median:.4f} ({low.label}) to "
                  f"{high.spearman_rho_median:.4f} ({high.label})"]
        for row in summary.itertuples():
            lines.append(f"- {row.label}: n drawn = {row.n_sampled_dmus}; rho p05 = {row.spearman_rho_p05:.4f}; "
                         f"median change in median TE = {row.median_score_change_median:+.4f}; "
                         f"maximum = {row.median_score_change_maximum:+.4f}; "
                         f"smallest individual change = {row.minimum_individual_score_change:+.2e}")
        others = metrics[~metrics["system"].eq("energy")]
        if len(others):
            lines.append(f"- Change in median TE below +0.01 in every replicate for the water and food "
                         f"frontiers: {bool((others['median_score_change'] < 0.01).all())}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "appendix_d_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ----------------------------------------------------------------------------- CLI
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--stages", nargs="+", choices=("all", "thresholds", "subsample", "summary"),
                        default=["all"])
    parser.add_argument("--minimum-dmus", type=int, default=MINIMUM_DMUS)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    parser.add_argument("--sample-fraction", type=float, default=SAMPLE_FRACTION)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--systems", nargs="+", choices=tuple(SYSTEMS), default=list(SYSTEMS),
                        help="Restrict the subsampling stage, e.g. for a quick check.")
    args = parser.parse_args()
    stages = ("thresholds", "subsample", "summary") if "all" in args.stages else args.stages
    OUT.mkdir(parents=True, exist_ok=True)
    write_design()
    write_run_config(args)
    if "thresholds" in stages:
        run_thresholds(args.minimum_dmus)
    if "subsample" in stages:
        run_subsample(args.replicates, args.sample_fraction, args.seed, args.workers, args.systems)
    if "summary" in stages or "all" in args.stages:
        write_summary()


if __name__ == "__main__":
    main()
