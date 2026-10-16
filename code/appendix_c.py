"""Appendix C: Table C1 and publication-style Figs. C1-C3."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd

from common import LAND_TYPES, ROOT, YEARS, read_csv, write_csv
from map_utils import add_grid_scores, draw_base


LAND_COLORS = {1: "#E6AB02", 2: "#1B7837", 4: "#7FBF7B", 8: "#D73027"}
YEAR_STYLES = {2010: ":", 2015: "--", 2020: "-"}


def _performance() -> dict[str, pd.DataFrame]:
    return {
        system: read_csv(ROOT / f"results/performance/grid_{system}.csv")
        for system in ("water", "energy", "food")
    }


def table_c1(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    specifications = (
        ("Water", frames["water"], "group_score", (1, 2, 4)),
        ("Energy", frames["energy"], "group_score", (8,)),
        ("Food", frames["food"], "group_score", (1,)),
    )
    for system, frame, score, codes in specifications:
        for code in codes:
            for year in YEARS:
                n = frame.loc[
                    frame["year"].eq(year) & frame["lulc_code"].eq(code), score
                ].notna().sum()
                rows.append(
                    {"system": system, "land_type": LAND_TYPES[code], "year": year, "n": n}
                )
    for year in YEARS:
        rows.extend(
            (
                {"system": "Water", "land_type": "Shrub", "year": year, "n": 0},
                {
                    "system": "Water",
                    "land_type": "All (meta-frontier)",
                    "year": year,
                    "n": frames["water"].loc[
                        frames["water"]["year"].eq(year), "meta_score"
                    ].notna().sum(),
                },
                {"system": "City scale", "land_type": "Each system", "year": year, "n": 76},
            )
        )
    result = pd.DataFrame(rows)
    result["pooled_frontier_n"] = result.groupby(
        ["system", "land_type"], observed=True
    )["n"].transform("sum")
    return result


def evaluated_cells(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    columns = [
        "dmu_id",
        "grid_id",
        "year",
        "x_center_m",
        "y_center_m",
        "lulc_code",
        "land_type",
        "lulc_dominant_share",
    ]
    water = frames["water"].loc[frames["water"]["meta_score"].notna(), columns]
    energy = frames["energy"].loc[frames["energy"]["group_score"].notna(), columns]
    return pd.concat([water, energy], ignore_index=True).drop_duplicates("dmu_id")


def figure_c1(frames: dict[str, pd.DataFrame]) -> None:
    evaluated = evaluated_cells(frames)
    rows = []
    for (year, code, land_type), group in evaluated.groupby(
        ["year", "lulc_code", "land_type"], observed=True
    ):
        share = group["lulc_dominant_share"].dropna()
        rows.append(
            {
                "year": year,
                "lulc_code": code,
                "land_type": land_type,
                "n": len(share),
                "median_share": share.median(),
                "n_share_ge_50": share.ge(0.5).sum(),
                "proportion_share_ge_50": share.ge(0.5).mean(),
                "n_share_ge_70": share.ge(0.7).sum(),
                "proportion_share_ge_70": share.ge(0.7).mean(),
            }
        )
    year_2020 = evaluated[evaluated["year"].eq(2020)]["lulc_dominant_share"].dropna()
    rows.append(
        {
            "year": 2020,
            "lulc_code": "All",
            "land_type": "All evaluated cells",
            "n": len(year_2020),
            "median_share": year_2020.median(),
            "n_share_ge_50": year_2020.ge(0.5).sum(),
            "proportion_share_ge_50": year_2020.ge(0.5).mean(),
            "n_share_ge_70": year_2020.ge(0.7).sum(),
            "proportion_share_ge_70": year_2020.ge(0.7).mean(),
        }
    )
    write_csv(pd.DataFrame(rows), ROOT / "results/appendix_c/Fig_C1_statistics.csv")

    fig, axes = plt.subplots(
        1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": (1.15, 1)}, layout="constrained"
    )
    map_data = evaluated[evaluated["year"].eq(2020)].copy()
    map_data["share_class"] = np.select(
        (
            map_data["lulc_dominant_share"].lt(0.5),
            map_data["lulc_dominant_share"].lt(0.7),
            map_data["lulc_dominant_share"].lt(0.9),
        ),
        (1, 2, 3),
        default=4,
    )
    class_colors = ("#DEEBF7", "#9ECAE1", "#4292C6", "#08519C")
    draw_base(axes[0])
    add_grid_scores(
        axes[0],
        map_data,
        "share_class",
        ListedColormap(class_colors),
        BoundaryNorm((0.5, 1.5, 2.5, 3.5, 4.5), 4),
    )
    axes[0].set_title("(a) Spatial distribution in 2020", loc="left", fontweight="bold")
    axes[0].text(
        0.01,
        0.03,
        f"n = {len(map_data):,} evaluated cells",
        transform=axes[0].transAxes,
        fontsize=9,
    )
    axes[0].legend(
        handles=[
            Patch(facecolor=color, label=label)
            for color, label in zip(class_colors, ("<0.5", "0.5–<0.7", "0.7–<0.9", "≥0.9"))
        ]
        + [Patch(facecolor="#ECECE9", edgecolor="#B0B0B0", label="Not evaluated")],
        title="Dominant-class share",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=5,
        frameon=False,
        fontsize=8,
    )

    for (code, year), group in evaluated.groupby(["lulc_code", "year"], observed=True):
        values = np.sort(group["lulc_dominant_share"].dropna().to_numpy())
        axes[1].step(
            values,
            np.arange(1, len(values) + 1) / len(values),
            where="post",
            color=LAND_COLORS[int(code)],
            linestyle=YEAR_STYLES[year],
            linewidth=1.4,
        )
    for threshold in (0.5, 0.7):
        axes[1].axvline(threshold, color="#BBBBBB", linestyle="--", linewidth=1)
        axes[1].text(threshold + 0.005, 0.97, f"{threshold:.0%}", color="#777777", va="top")
    axes[1].set(
        title="(b) Cumulative distribution by land-use type and year",
        xlabel="Share occupied by dominant land-use class",
        ylabel="Cumulative proportion of evaluated cells",
        xlim=(0.25, 1),
        ylim=(0, 1),
    )
    axes[1].grid(color="#E2E2E2", linewidth=0.5)
    land_handles = [
        Line2D([0], [0], color=LAND_COLORS[code], label=LAND_TYPES[code])
        for code in (1, 2, 4, 8)
    ]
    year_handles = [
        Line2D([0], [0], color="#333333", linestyle=YEAR_STYLES[year], label=str(year))
        for year in YEARS
    ]
    first = axes[1].legend(handles=year_handles, title="Year", loc="upper left", frameon=False)
    axes[1].add_artist(first)
    axes[1].legend(handles=land_handles, title="Land-use type", loc="lower right", frameon=False)
    _save(fig, "Fig_C1")


def _distribution_panel(axis, frame, score, codes, title, rng) -> None:
    values, labels, colors, positions = [], [], [], []
    position = 1
    for year in YEARS:
        for code in codes:
            group = frame.loc[
                frame["year"].eq(year) & frame["lulc_code"].eq(code), score
            ].dropna()
            values.append(group.to_numpy())
            labels.append(f"{LAND_TYPES[code]}\nn={len(group):,}")
            colors.append(LAND_COLORS[code])
            positions.append(position)
            axis.scatter(
                rng.normal(position, 0.07, len(group)),
                group,
                s=2,
                alpha=0.18,
                color=LAND_COLORS[code],
                linewidths=0,
            )
            position += 1
    boxes = axis.boxplot(
        values,
        positions=positions,
        tick_labels=labels,
        showfliers=False,
        patch_artist=True,
        widths=0.55,
    )
    for box, color in zip(boxes["boxes"], colors):
        box.set_facecolor(color)
        box.set_alpha(0.35)
        box.set_edgecolor(color)
    for separator in range(len(codes), len(values), len(codes)):
        axis.axvline(separator + 0.5, color="#CCCCCC", linestyle=":", linewidth=0.8)
    for year_index, year in enumerate(YEARS):
        centre = year_index * len(codes) + (len(codes) + 1) / 2
        axis.text(centre, 0.985, str(year), ha="center", va="top", fontweight="bold")
    axis.set(title=title, ylim=(0.485, 1.012))
    axis.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    axis.tick_params(axis="x", labelrotation=47, labelsize=7)
    axis.grid(axis="y", color="#DDDDDD", linewidth=0.5)


def figure_c2(frames: dict[str, pd.DataFrame]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharey=True, layout="constrained")
    rng = np.random.default_rng(2602569)
    _distribution_panel(
        axes[0, 0], frames["water"], "group_score", (1, 2, 4),
        r"(a) Water: group-frontier score ($TE^{G}$)", rng,
    )
    _distribution_panel(
        axes[0, 1], frames["water"], "meta_score", (1, 2, 4),
        r"(b) Water: meta-frontier score ($TE^{M}$)", rng,
    )
    _distribution_panel(
        axes[1, 0], frames["energy"], "group_score", (8,),
        r"(c) Energy: group-frontier score ($TE^{G}$), Impervious", rng,
    )
    _distribution_panel(
        axes[1, 1], frames["food"], "group_score", (1,),
        r"(d) Food: group-frontier score ($TE^{G}$), Cropland", rng,
    )
    axes[0, 0].set_ylabel("Resource-performance score")
    axes[1, 0].set_ylabel("Resource-performance score")
    fig.legend(
        handles=[
            Patch(facecolor=LAND_COLORS[code], alpha=0.35, label=LAND_TYPES[code])
            for code in (1, 2, 4, 8)
        ],
        title="Dominant land-use type",
        loc="outside lower center",
        ncol=4,
        frameon=False,
    )
    _save(fig, "Fig_C2")


def figure_c3(frames: dict[str, pd.DataFrame]) -> None:
    water = frames["water"]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7), layout="constrained")
    image = None
    panel = 0
    for row, code in enumerate((2, 4)):
        for column, year in enumerate((2010, 2020)):
            group = water[
                water["lulc_code"].eq(code) & water["year"].eq(year)
            ].dropna(subset=["group_score"])
            draw_base(axes[row, column])
            image = add_grid_scores(
                axes[row, column], group, "group_score", plt.get_cmap("Blues"), Normalize(0.5, 1)
            )
            axes[row, column].set_title(str(year), fontweight="bold")
            axes[row, column].text(
                0,
                1.02,
                f"({chr(97 + panel)})",
                transform=axes[row, column].transAxes,
                fontweight="bold",
                fontsize=12,
            )
            axes[row, column].text(
                0.02,
                0.90,
                f"n = {len(group):,}",
                transform=axes[row, column].transAxes,
                bbox={"facecolor": "white", "edgecolor": "#BBBBBB", "boxstyle": "round"},
                fontsize=8,
            )
            panel += 1
        axes[row, 0].text(
            -0.08,
            0.5,
            LAND_TYPES[code],
            transform=axes[row, 0].transAxes,
            rotation=90,
            va="center",
            fontweight="bold",
            fontsize=13,
        )
    fig.colorbar(
        image,
        ax=axes,
        orientation="horizontal",
        label=r"Water group-frontier score ($TE^{G}$)",
        shrink=0.55,
        pad=0.02,
    )
    _save(fig, "Fig_C3")


def _save(fig, name: str) -> None:
    output = ROOT / "figures"
    output.mkdir(exist_ok=True)
    fig.savefig(output / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def run_appendix_c() -> None:
    frames = _performance()
    write_csv(table_c1(frames), ROOT / "results/appendix_c/Table_C1.csv")
    figure_c1(frames)
    figure_c2(frames)
    figure_c3(frames)


if __name__ == "__main__":
    run_appendix_c()
