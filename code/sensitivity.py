"""Appendix D sensitivity: dominant, 50%, and 70% land-share rules."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analysis_4_1 import plot_grid_city_figure
from analysis_4_2 import plot_section_4_2, region_code
from analysis_4_3 import _ellipse_panel, sde_legend_handles, sde_tables
from common import (
    PRIMARY_SCORE,
    REGIONS,
    ROOT,
    SYSTEMS,
    describe,
    grid_candidates,
    load_grid,
    read_csv,
    write_csv,
)
from dea_performance import score_grid_frontiers


FIGURE_NUMBER = {"threshold_50": ("D1", "D3", "D5"), "threshold_70": ("D2", "D4", "D6")}


def _correlation(first: pd.Series, second: pd.Series) -> float:
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


def _figure_scale(rule: str, scenario: dict[str, pd.DataFrame], city) -> None:
    """Figs. D1-D2: grid and unchanged city scores for all systems and years."""
    threshold = {"threshold_50": "≥50%", "threshold_70": "≥70%"}[rule]
    plot_grid_city_figure(
        scenario,
        city,
        f"Fig_{FIGURE_NUMBER[rule][0]}",
        f"Resource performance under the dominant land-use share {threshold} rule",
    )


def _figure_regions(rule: str, scenario: dict[str, pd.DataFrame]) -> None:
    """Figs. D3-D4: repeat Section 4.2 under one threshold rule."""
    systems = {}
    for name in ("water", "food"):
        frame = scenario[name]
        frame = frame[frame["year"].eq(2020) & frame["lulc_code"].eq(1)].copy()
        frame["region_code"] = region_code(frame)
        systems[name] = frame
    plot_section_4_2(systems, f"Fig_{FIGURE_NUMBER[rule][1]}")


def _figure_sde(rule: str, parameters: pd.DataFrame) -> None:
    """Figs. D5-D6: repeat the grid-scale SDE analysis under one rule."""
    fig, axes = plt.subplots(1, 3, figsize=(11, 4), layout="constrained")
    for index, name in enumerate(SYSTEMS):
        _ellipse_panel(
            axes[index],
            parameters[parameters["system"].eq(name)],
            f"({chr(97 + index)}) {name.capitalize()}",
        )
    fig.legend(
        handles=sde_legend_handles(),
        loc="outside upper center",
        ncol=5,
        frameon=False,
    )
    _save(fig, FIGURE_NUMBER[rule][2])


def _d1_rows(frame: pd.DataFrame, rule: str, system: str) -> list[dict]:
    rows = []
    score = PRIMARY_SCORE[system]
    for (year, land_type), group in frame.groupby(["year", "land_type"], observed=True):
        summary = describe(group[score])
        rows.append(
            {
                "lulc_rule": rule,
                "system": system,
                "land_type": land_type,
                "year": year,
                **summary,
                "reported": summary["n"] >= 9,
            }
        )
    return rows


def _table_d2(correlations, regional, moves) -> pd.DataFrame:
    rows = []
    for row in correlations:
        rows.append(
            {
                "lulc_rule": row["lulc_rule"],
                "result": "grid_rank_correlation",
                "system": row["system"],
                "year": row["year"],
                "group_or_region": row["land_type"],
                "value": row["spearman_rho"],
                "direction": pd.NA,
            }
        )
    for row in regional:
        rows.append(
            {
                "lulc_rule": row["lulc_rule"],
                "result": "regional_median",
                "system": row["system"],
                "year": row["year"],
                "group_or_region": row["region_code"],
                "value": row["median"],
                "direction": pd.NA,
            }
        )
    for frame in moves:
        for row in frame.itertuples(index=False):
            rows.append(
                {
                    "lulc_rule": row.scenario,
                    "result": "sde_displacement_2010_2020",
                    "system": row.system,
                    "year": 2020,
                    "group_or_region": "grid_scale",
                    "value": row.geodesic_displacement_km,
                    "direction": row.movement_direction_8point,
                }
            )
    return pd.DataFrame(rows)


def run_sensitivity() -> None:
    grid = load_grid()
    baseline = {
        name: read_csv(ROOT / f"results/performance/grid_{name}.csv") for name in SYSTEMS
    }
    city = {
        name: read_csv(ROOT / f"results/performance/city_{name}.csv") for name in SYSTEMS
    }
    d1 = [row for name, frame in baseline.items() for row in _d1_rows(frame, "dominant", name)]
    correlations, regional, all_parameters, all_moves = [], [], [], []

    for rule in ("threshold_50", "threshold_70"):
        scenario, references = {}, []
        for name, system in SYSTEMS.items():
            scores, reference = score_grid_frontiers(grid_candidates(grid, rule, name), system)
            scenario[name] = scores
            reference.insert(0, "lulc_rule", rule)
            references.append(reference)
            write_csv(scores, ROOT / f"results/sensitivity/{rule}_{name}.csv")
            d1.extend(_d1_rows(scores, rule, name))

            score = PRIMARY_SCORE[name]
            primary = baseline[name][["dmu_id", "year", "land_type", score]].rename(
                columns={score: "dominant_score"}
            )
            alternative = scores[["dmu_id", score]].rename(columns={score: "scenario_score"})
            common = primary.merge(alternative, on="dmu_id", how="inner", validate="one_to_one")
            for (year, land_type), group in common.groupby(
                ["year", "land_type"], observed=True
            ):
                correlations.append(
                    {
                        "lulc_rule": rule,
                        "system": name,
                        "year": year,
                        "land_type": land_type,
                        "n_common_dmus": len(group),
                        "spearman_rho": _correlation(
                            group["dominant_score"], group["scenario_score"]
                        ),
                    }
                )
        write_csv(
            pd.concat(references, ignore_index=True),
            ROOT / f"results/sensitivity/{rule}_references.csv",
        )

        for name in ("water", "food"):
            cropland = scenario[name]
            cropland = cropland[
                cropland["year"].eq(2020) & cropland["lulc_code"].eq(1)
            ].copy()
            cropland["region_code"] = region_code(cropland)
            for label, region_name, *_ in REGIONS:
                regional.append(
                    {
                        "lulc_rule": rule,
                        "system": name,
                        "year": 2020,
                        "region_code": label,
                        "region_name": region_name,
                        **describe(
                            cropland.loc[cropland["region_code"].eq(label), "group_score"]
                        ),
                    }
                )

        parameters, moves = sde_tables(scenario, PRIMARY_SCORE, rule)
        all_parameters.append(parameters)
        all_moves.append(moves)
        _figure_scale(rule, scenario, city)
        _figure_regions(rule, scenario)
        _figure_sde(rule, parameters)

    d1 = pd.DataFrame(d1)
    d1.loc[~d1["reported"], ["median"]] = np.nan
    write_csv(d1, ROOT / "results/sensitivity/Table_D1_long.csv")
    d1["table_value"] = d1.apply(
        lambda row: (
            f"{int(row['n']):,} ({row['median']:.3f})" if row["reported"] else "—"
        ),
        axis=1,
    )
    table_d1 = (
        d1.pivot(
            index=["system", "land_type", "year"],
            columns="lulc_rule",
            values="table_value",
        )
        .reset_index()
        .rename_axis(columns=None)
    )
    table_d1 = table_d1[
        ["system", "land_type", "year", "dominant", "threshold_50", "threshold_70"]
    ]
    write_csv(table_d1, ROOT / "results/sensitivity/Table_D1.csv")
    write_csv(pd.DataFrame(correlations), ROOT / "results/sensitivity/grid_rank_correlations.csv")
    write_csv(pd.DataFrame(regional), ROOT / "results/sensitivity/regional_scores.csv")
    write_csv(pd.concat(all_parameters, ignore_index=True), ROOT / "results/sensitivity/sde_parameters.csv")
    write_csv(pd.concat(all_moves, ignore_index=True), ROOT / "results/sensitivity/sde_displacements.csv")
    write_csv(_table_d2(correlations, regional, all_moves), ROOT / "results/sensitivity/Table_D2.csv")


if __name__ == "__main__":
    run_sensitivity()
