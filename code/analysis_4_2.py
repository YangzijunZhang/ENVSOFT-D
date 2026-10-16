"""Section 4.2: 2020 Cropland Water and Food performance in Regions A-D."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter
import pandas as pd

from common import REGIONS, ROOT, describe, read_csv, write_csv
from map_utils import add_grid_scores, draw_base, performance_cmap


REGION_TICKS = {
    "A": "A\nHetao Irrigation\nDistrict",
    "B": "B\nWestern Loess\nPlateau",
    "C": "C\nGuanzhong\nPlain",
    "D": "D\nLower Yellow\nRiver Plain",
}


def region_code(frame: pd.DataFrame) -> pd.Series:
    code = pd.Series(pd.NA, index=frame.index, dtype="string")
    for label, _, xmin, xmax, ymin, ymax in REGIONS:
        inside = frame["x_center_m"].between(xmin, xmax) & frame["y_center_m"].between(
            ymin, ymax
        )
        code.loc[inside] = label
    return code


def _draw_regions(axis) -> None:
    for label, _, xmin, xmax, ymin, ymax in REGIONS:
        axis.add_patch(
            Rectangle(
                (xmin, ymin),
                xmax - xmin,
                ymax - ymin,
                fill=False,
                edgecolor="#111111",
                linewidth=1.0,
                zorder=8,
            )
        )
        axis.text(
            xmin + 0.02 * (xmax - xmin),
            ymax - 0.04 * (ymax - ymin),
            label,
            va="top",
            fontsize=9,
            zorder=9,
        )


def plot_section_4_2(systems: dict[str, pd.DataFrame], output_name: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), layout="constrained")
    titles = {"water": "Water", "food": "Food"}
    box_colors = {"water": "#73A2CA", "food": "#82B997"}
    norm = Normalize(0.5, 1)
    for column, name in enumerate(("water", "food")):
        frame = systems[name]
        cmap = performance_cmap(name)
        draw_base(axes[0, column])
        image = add_grid_scores(axes[0, column], frame, "group_score", cmap, norm)
        _draw_regions(axes[0, column])
        axes[0, column].set_title(
            f"({chr(97 + column)}) {titles[name]} performance (2020)",
            loc="left",
            fontweight="bold",
        )
        fig.colorbar(
            image,
            ax=axes[0, column],
            orientation="horizontal",
            label="Group-frontier TE",
            fraction=0.045,
            pad=0.02,
        )

        values, labels = [], []
        for label, region_name, *_ in REGIONS:
            group = frame.loc[frame["region_code"].eq(label), "group_score"].dropna()
            values.append(group.to_numpy())
            labels.append(f"{REGION_TICKS[label]}\n(n = {len(group):,})")
        boxes = axes[1, column].boxplot(
            values, tick_labels=labels, patch_artist=True, showfliers=True
        )
        for box in boxes["boxes"]:
            box.set_facecolor(box_colors[name])
            box.set_alpha(0.75)
        basin_median = frame["group_score"].median()
        axes[1, column].axhline(
            basin_median, color="#777777", linestyle="--", linewidth=0.8
        )
        axes[1, column].text(
            4.32,
            basin_median,
            f"Basin median\n{basin_median:.3f}",
            va="center",
            fontsize=8,
            color="#666666",
        )
        axes[1, column].set(
            title=f"({chr(99 + column)}) {titles[name]} performance by agricultural region",
            ylabel="Group-frontier TE",
            ylim=(0.48, 1.02),
            xlim=(0.5, 4.85),
        )
        axes[1, column].yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        axes[1, column].tick_params(axis="x", labelsize=7)
        axes[1, column].grid(axis="y", color="#DDDDDD", linewidth=0.6)

    output = ROOT / "figures"
    output.mkdir(exist_ok=True)
    fig.savefig(output / f"{output_name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output / f"{output_name}.pdf", bbox_inches="tight")
    plt.close(fig)


def run_section_4_2() -> None:
    systems, summaries = {}, []
    for system in ("water", "food"):
        frame = read_csv(ROOT / f"results/performance/grid_{system}.csv")
        frame = frame[frame["year"].eq(2020) & frame["lulc_code"].eq(1)].copy()
        frame["region_code"] = region_code(frame)
        systems[system] = frame
        for label, name, *_ in REGIONS:
            summaries.append(
                {
                    "system": system,
                    "year": 2020,
                    "region_code": label,
                    "region_name": name,
                    **describe(frame.loc[frame["region_code"].eq(label), "group_score"]),
                }
            )
    write_csv(pd.DataFrame(summaries), ROOT / "results/analysis_4_2/regional_scores.csv")
    plot_section_4_2(systems, "Fig_6")


if __name__ == "__main__":
    run_section_4_2()
