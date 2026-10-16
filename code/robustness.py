"""Appendix D robustness: 100 stratified 80% frontier subsamples."""

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import LAND_TYPES, ROOT, SYSTEMS, grid_candidates, load_grid, read_csv, write_csv
from dea_performance import score_grid_frontiers


SEED = 2602569
REPLICATES = 100
SAMPLE_FRACTION = 0.80
FRONTIERS = (
    ("water", "group", "Cropland", "Water–Cropland"),
    ("water", "group", "Forest", "Water–Forest"),
    ("water", "group", "Grassland", "Water–Grassland"),
    ("water", "meta", "Cropland,Forest,Grassland", "Water meta-frontier"),
    ("energy", "group", "Impervious", "Energy–Impervious"),
    ("food", "group", "Cropland", "Food–Cropland"),
)


def _sample(frame: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    indices = []
    for _, group in frame.groupby(["year", "lulc_code"], observed=True, sort=True):
        n = math.floor(len(group) * SAMPLE_FRACTION)
        indices.extend(rng.choice(group.index.to_numpy(), n, replace=False))
    return frame.loc[indices].copy()


def _rho(first: pd.Series, second: pd.Series) -> float:
    if len(first) < 2 or first.nunique() < 2 or second.nunique() < 2:
        return np.nan
    return float(spearmanr(first, second).statistic)


def _figure_d7(metrics: pd.DataFrame) -> None:
    groups = []
    labels = []
    for system, frontier, land_type, label in FRONTIERS:
        groups.append(
            metrics[
                metrics["system"].eq(system)
                & metrics["frontier"].eq(frontier)
                & metrics["land_type_or_groups"].eq(land_type)
            ]
        )
        labels.append(label)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for axis, column, title in (
        (axes[0], "spearman_rho", "(a) Rank correlation"),
        (axes[1], "median_score_change", "(b) Change in median TE"),
    ):
        boxes = axis.boxplot(
            [group[column].dropna().to_numpy() for group in groups],
            tick_labels=labels,
            showfliers=True,
            patch_artist=True,
        )
        for box in boxes["boxes"]:
            box.set_facecolor("#78A6D1")
            box.set_alpha(0.75)
        axis.set_title(title)
        axis.tick_params(axis="x", labelrotation=45, labelsize=7)
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    axes[0].set_ylabel("Spearman rank correlation")
    axes[1].set_ylabel("Subsample minus full-sample median TE")
    axes[1].axhline(0, color="#555555", linewidth=0.8)
    output = ROOT / "figures"
    output.mkdir(exist_ok=True)
    fig.savefig(output / "Fig_D7.png", dpi=300, bbox_inches="tight")
    fig.savefig(output / "Fig_D7.pdf", bbox_inches="tight")
    plt.close(fig)


def run_robustness() -> None:
    grid = load_grid()
    candidates = {
        name: grid_candidates(grid, "dominant", name).query("analysis_eligible").copy()
        for name in SYSTEMS
    }
    full_scores = {
        name: read_csv(ROOT / f"results/performance/grid_{name}.csv")
        for name in SYSTEMS
    }
    tasks = [(name, replicate) for name in SYSTEMS for replicate in range(1, REPLICATES + 1)]
    seeds = np.random.SeedSequence(SEED).spawn(len(tasks))
    metrics = []

    for (name, replicate), child_seed in zip(tasks, seeds):
        system = SYSTEMS[name]
        sampled = _sample(candidates[name], np.random.default_rng(child_seed))
        scores, _ = score_grid_frontiers(sampled, system)
        comparisons = [
            ("group", code, LAND_TYPES[code], "group_score") for code in system.land_codes
        ]
        if system.meta:
            comparisons.append(("meta", None, "Cropland,Forest,Grassland", "meta_score"))

        for frontier, code, label, score_column in comparisons:
            subsample = scores if code is None else scores[scores["lulc_code"].eq(code)]
            complete = (
                full_scores[name]
                if code is None
                else full_scores[name][full_scores[name]["lulc_code"].eq(code)]
            )
            common = complete[["scenario_dmu_id", score_column]].merge(
                subsample[["scenario_dmu_id", score_column]],
                on="scenario_dmu_id",
                how="inner",
                suffixes=("_full", "_subsample"),
                validate="one_to_one",
            ).dropna()
            metrics.append(
                {
                    "replicate": replicate,
                    "system": name,
                    "frontier": frontier,
                    "land_type_or_groups": label,
                    "n_sampled_dmus": len(common),
                    "spearman_rho": _rho(
                        common[f"{score_column}_full"], common[f"{score_column}_subsample"]
                    ),
                    "full_median_on_sampled_dmus": common[f"{score_column}_full"].median(),
                    "subsample_median": common[f"{score_column}_subsample"].median(),
                    "median_score_change": common[f"{score_column}_subsample"].median()
                    - common[f"{score_column}_full"].median(),
                }
            )

    metrics = pd.DataFrame(metrics)
    write_csv(metrics, ROOT / "results/robustness/replicate_metrics.csv")
    summary = []
    keys = ["system", "frontier", "land_type_or_groups"]
    for key, group in metrics.groupby(keys, observed=True):
        row = dict(zip(keys, key))
        row["n_replicates"] = len(group)
        for column in ("spearman_rho", "median_score_change"):
            values = group[column].dropna()
            row[f"{column}_minimum"] = values.min()
            row[f"{column}_p05"] = values.quantile(0.05)
            row[f"{column}_median"] = values.median()
            row[f"{column}_p95"] = values.quantile(0.95)
            row[f"{column}_maximum"] = values.max()
        summary.append(row)
    write_csv(pd.DataFrame(summary), ROOT / "results/robustness/summary_metrics.csv")
    _figure_d7(metrics)


if __name__ == "__main__":
    run_robustness()
