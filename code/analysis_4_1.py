"""Section 4.1: grid- and city-scale Water, Energy and Food maps."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from common import ROOT, SYSTEMS, YEARS, read_csv
from map_utils import add_city_scores, add_grid_scores, draw_base, performance_cmap


def plot_grid_city_figure(grid, city, output_name: str, title: str | None = None) -> None:
    fig, axes = plt.subplots(6, 3, figsize=(12, 15), layout="constrained")
    norm = Normalize(0.5, 1)
    for system_index, name in enumerate(SYSTEMS):
        cmap = performance_cmap(name)
        grid_row, city_row = 2 * system_index, 2 * system_index + 1
        image = None
        for column, year in enumerate(YEARS):
            grid_year = grid[name][grid[name]["year"].eq(year)]
            draw_base(axes[grid_row, column])
            image = add_grid_scores(
                axes[grid_row, column], grid_year, "group_score", cmap, norm
            )
            city_year = city[name][city[name]["year"].eq(year)]
            add_city_scores(axes[city_row, column], city_year, "city_score", cmap, norm)
            if system_index == 0:
                axes[grid_row, column].set_title(str(year), fontsize=13)
        axes[grid_row, 0].text(
            -0.02,
            1.03,
            f"({chr(97 + system_index)}) {name.capitalize()}",
            transform=axes[grid_row, 0].transAxes,
            fontsize=13,
            fontweight="bold",
        )
        axes[grid_row, 0].text(
            0,
            0.97,
            "Grid-scale group frontier",
            transform=axes[grid_row, 0].transAxes,
            va="top",
            fontsize=9,
        )
        axes[city_row, 0].text(
            0,
            0.97,
            "City-scale pooled frontier",
            transform=axes[city_row, 0].transAxes,
            va="top",
            fontsize=9,
        )
        fig.colorbar(
            image,
            ax=axes[grid_row : city_row + 1, :],
            label=f"{name.capitalize()} TE",
            shrink=0.68,
            pad=0.01,
        )
    fig.suptitle(
        title
        or "Water, Energy and Food performance at grid and city scales "
        "(2010, 2015 and 2020)",
        fontsize=16,
    )
    output = ROOT / "figures"
    output.mkdir(exist_ok=True)
    fig.savefig(output / f"{output_name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output / f"{output_name}.pdf", bbox_inches="tight")
    plt.close(fig)


def run_section_4_1() -> None:
    grid = {
        name: read_csv(ROOT / f"results/performance/grid_{name}.csv") for name in SYSTEMS
    }
    city = {
        name: read_csv(
            ROOT / f"results/performance/city_{name}.csv", dtype={"city_id": str}
        )
        for name in SYSTEMS
    }
    plot_grid_city_figure(grid, city, "wef_grid_city_performance_3years_combined")


if __name__ == "__main__":
    run_section_4_1()
