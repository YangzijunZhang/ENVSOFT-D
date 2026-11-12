"""Appendix E: specification diagnostics (Tables E1-E6).

Reads data/ and the primary scores in results/performance/
(run code/dea_performance.py --scale both first). About 8 minutes.

Outputs (results/appendix_e/)
    Table_E1.csv      E.1  share of basin totals held by the evaluated energy cells
    Table_E2.csv      E.2  scores versus simple output-input ratios
    Table_E3.csv      E.2  soil erosion beyond NPP/ET (water meta-frontier)
    Table_E4.csv      E.3  constant returns to scale and strong disposability
    E3_zhengzhou.csv  E.3  Zhengzhou city scores under VRS and CRS
    Table_E5.csv      E.4  excluding grid cells < 10 km2 and cities < 100 km2 inside the basin
    Table_E6.csv      E.5  excluding water DMUs with > 20% non-vegetated land
"""

import numpy as np
import pandas as pd
from pyproj import Geod, Transformer
from scipy.optimize import linprog
from scipy.stats import pearsonr, spearmanr

from analysis_4_3 import add_center_context
from common import ALBERS_CRS, PRIMARY_SCORE, ROOT, SYSTEMS, YEARS, grid_candidates, load_grid, read_csv, write_csv
from dea_performance import build_reference, score_city_frontier, score_grid_frontiers

OUT = ROOT / "results" / "appendix_e"
FRONTIER = 1 - 1e-6  # score treated as on the frontier
NONVEGETATED = [f"lulc_fraction_{c}" for c in ("5_water", "6_snow_ice", "7_barren", "8_impervious", "9_wetland")]
GEOD = Geod(ellps="WGS84")
TO_WGS84 = Transformer.from_crs(ALBERS_CRS, "EPSG:4326", always_xy=True)


# ---------------------------------------------------------------- helpers
def eligible(frame):
    return frame[frame["analysis_eligible"].fillna(False).astype(bool)].copy()


def spearman(a, b):
    return float(spearmanr(a, b, nan_policy="omit")[0])


def log_pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = (x > 0) & (y > 0)
    return float(pearsonr(np.log(x[keep]), np.log(y[keep]))[0])


def frontiers(frame, system):
    """(label, score column, DMUs) for each grid frontier of a system."""
    for land_type, group in frame.groupby("land_type", observed=True):
        yield f"group: {land_type}", "group_score", group
    if system.meta:
        yield "meta", "meta_score", frame


def centre(frame, score):
    frame = frame[frame[score] > 0]
    return np.average(frame[["x_center_m", "y_center_m"]], axis=0, weights=frame[score])


def shift(frame, score):
    """Geodesic distance (km) and azimuth (degrees from north) of the centre, 2010-2020."""
    (x0, y0), (x1, y1) = centre(frame[frame["year"].eq(2010)], score), centre(frame[frame["year"].eq(2020)], score)
    (lon0, lon1), (lat0, lat1) = TO_WGS84.transform([x0, x1], [y0, y1])
    azimuth, _, distance = GEOD.inv(lon0, lat0, lon1, lat1)
    return distance / 1000, azimuth % 360


def variant_scores(frame, system, crs=False, strong=False, no_bad=False):
    """Observed-direction DDF scores under an alternative specification.

    crs: no convexity constraint; strong: inequality on the undesirable output;
    no_bad: undesirable output omitted. Failed solves are retried as in dea_performance.py.
    """
    reduced = build_reference(frame, system)
    attempts = [(reduced, "highs")] + [(reduced.full(), m) for m in ("highs", "highs-ds", "highs-ipm")]
    scores = []
    for x, y, b in frame[list(system.variables)].to_numpy(float):
        score = np.nan
        for ref, method in attempts:
            n = len(ref.ids)
            a_ub, b_ub = [np.r_[ref.x[:, 0], x], np.r_[-ref.y[:, 0], y]], [x, -y]
            a_eq, b_eq = [], []
            if not no_bad:
                (a_ub if strong else a_eq).append(np.r_[ref.b[:, 0], b])
                (b_ub if strong else b_eq).append(b)
            if not crs:
                a_eq.append(np.r_[np.ones(n), 0.0])
                b_eq.append(1.0)
            lp = linprog(np.r_[np.zeros(n), -1.0], A_ub=a_ub, b_ub=b_ub, A_eq=a_eq or None, b_eq=b_eq or None,
                         bounds=[(0, None)] * n + [(0, 1)], method=method)
            if lp.success:
                score = 1 / (1 + lp.x[-1])
                break
        scores.append(score)
    return np.array(scores)


def rescore_grid(panel, name, keep):
    """Re-estimate the grid frontiers of one system with only the DMUs where keep() is True."""
    candidates = grid_candidates(panel, "dominant", name)
    candidates["analysis_eligible"] = candidates["analysis_eligible"].fillna(False).astype(bool) & keep(candidates)
    return eligible(score_grid_frontiers(candidates, SYSTEMS[name])[0])


def compare(base, new, system, score_of, rule):
    """Rows of Tables E5/E6: excluded DMUs, rank correlation and SDE shift before and after."""
    main_km, main_az = shift(base, score_of)
    new_km, new_az = shift(new, score_of)
    new = new.set_index("dmu_id")
    rows = []
    for label, score, subset in frontiers(base, system) if score_of != "city_score" else [("city", "city_score", base)]:
        kept = pd.Index(subset["dmu_id"]).intersection(new.index)
        on_frontier = subset[subset[score] >= FRONTIER]["dmu_id"]
        rows.append({"rule": rule, "frontier": label, "n_main": len(subset), "n_excluded": len(subset) - len(kept),
                     "n_frontier_main": len(on_frontier), "frontier_dmus_excluded": int((~on_frontier.isin(new.index)).sum()),
                     "spearman_retained": spearman(subset.set_index("dmu_id").loc[kept, score], new.loc[kept, score]),
                     "sde_main_km": main_km, "sde_main_az": main_az, "sde_new_km": new_km, "sde_new_az": new_az})
    return rows


# ---------------------------------------------------------------- main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    panel = load_grid()
    grid = {name: eligible(read_csv(ROOT / f"results/performance/grid_{name}.csv", low_memory=False)) for name in SYSTEMS}
    city_panel = read_csv(ROOT / "data/city_panel_10km.csv", dtype={"city_id": str, "province_id": str})
    city = {}
    for name, system in SYSTEMS.items():
        result = score_city_frontier(city_panel, system)[0]
        city[name] = result[result["city_common_score_sample"]].assign(dmu_id=lambda f: f["scenario_dmu_id"])

    # E.1 Coverage of the evaluated energy cells
    energy = SYSTEMS["energy"]
    rows = []
    for year in YEARS:
        basin, cells = panel[panel["year"].eq(year)], grid["energy"][grid["energy"]["year"].eq(year)]
        rows.append({"year": year, "n_cells": len(cells),
                     **{f"{v}_share_pct": 100 * cells[c].sum() / basin[c].sum()
                        for v, c in (("electricity", energy.input), ("gdp", energy.good), ("ghg", energy.bad))}})
    write_csv(pd.DataFrame(rows), OUT / "Table_E1.csv")
    print("Table E1 done", flush=True)

    # E.2 Scores versus output-input ratios
    rows = []
    for name, system in SYSTEMS.items():
        units = [(label, score, subset) for label, score, subset in frontiers(grid[name], system)]
        units.append(("city", "city_score", city[name]))
        for label, score, frame in units:
            rows.append({"scale": "city" if label == "city" else "grid", "system": name, "frontier": label, "n": len(frame),
                         "spearman_score_vs_ratio": spearman(frame[score], frame[system.good] / frame[system.input]),
                         "pearson_log_input_vs_log_good": log_pearson(frame[system.input], frame[system.good])})
    write_csv(pd.DataFrame(rows), OUT / "Table_E2.csv")

    water, w = SYSTEMS["water"], grid["water"].copy()
    ratio, erosion = w[water.good] / w[water.input], w[water.bad] / w[water.input]
    no_erosion = pd.Series(variant_scores(w, water, no_bad=True), index=w.index)
    strong = pd.Series(variant_scores(w, water, strong=True), index=w.index)
    pct = lambda s: s.rank(pct=True)
    top = erosion > erosion.quantile(0.9)

    def coefficient(scores):  # rank regression: score rank ~ NPP/ET rank + erosion/ET rank
        design = np.c_[np.ones(len(w)), pct(ratio), pct(erosion)]
        return float(np.linalg.lstsq(design, pct(scores), rcond=None)[0][2])

    stats = {
        "Spearman meta score vs NPP/ET": spearman(w["meta_score"], ratio),
        "Spearman meta score without erosion vs NPP/ET": spearman(no_erosion, ratio),
        "Share of DMUs with percentile-rank gap > 10 points (%)": 100 * float(((pct(w["meta_score"]) - pct(ratio)).abs() > 0.1).mean()),
        "Rank-regression coefficient on erosion/ET, weak disposability": coefficient(w["meta_score"]),
        "Rank-regression coefficient on erosion/ET, strong disposability": coefficient(strong),
        "Mean percentile-rank gap of top-10% erosion/ET DMUs, weak (points)": 100 * float((pct(w["meta_score"]) - pct(ratio))[top].mean()),
        "Mean percentile-rank gap of top-10% erosion/ET DMUs, strong (points)": 100 * float((pct(strong) - pct(ratio))[top].mean()),
    }
    write_csv(pd.DataFrame(stats.items(), columns=["statistic", "value"]), OUT / "Table_E3.csv")
    print("Tables E2-E3 done", flush=True)

    # E.3 Constant returns to scale and strong disposability
    rows = []
    for name, system in SYSTEMS.items():
        units = list(frontiers(grid[name], system)) + [("city", "city_score", city[name])]
        for label, score, frame in units:
            crs, strong = variant_scores(frame, system, crs=True), variant_scores(frame, system, strong=True)
            rows.append({"scale": "city" if label == "city" else "grid", "system": name, "frontier": label, "n": len(frame),
                         "spearman_crs": spearman(frame[score], crs), "spearman_strong": spearman(frame[score], strong),
                         "frontier_pct_main": 100 * float((frame[score] >= FRONTIER).mean()),
                         "frontier_pct_crs": 100 * float((crs >= FRONTIER).mean()),
                         "frontier_pct_strong": 100 * float((strong >= FRONTIER).mean())})
            if name == "energy" and label == "city":
                zz = frame.assign(city_score_crs=crs)
                write_csv(zz[zz["city_name"].eq("Zhengzhou")][["city_name", "year", "city_score", "city_score_crs"]],
                          OUT / "E3_zhengzhou.csv")
    write_csv(pd.DataFrame(rows), OUT / "Table_E4.csv")
    print("Table E4 done", flush=True)

    # E.4 Small grid cells and cities
    rows = []
    for name, system in SYSTEMS.items():
        new = rescore_grid(panel, name, lambda f: f["area_inside_yrb_km2"].ge(10))
        rows += [{"scale": "grid", "system": name, **r} for r in compare(grid[name], new, system, PRIMARY_SCORE[name], "grid cells < 10 km2")]
    small = city_panel["city_area_inside_yrb_km2"] < 100
    for name, system in SYSTEMS.items():
        new = score_city_frontier(city_panel[~small], system)[0]
        new = new[new["city_common_score_sample"]].assign(dmu_id=lambda f: f["scenario_dmu_id"])
        rows += [{"scale": "city", "system": name, **r} for r in compare(city[name], new, system, "city_score", "cities < 100 km2")]
    table = pd.DataFrame(rows)
    table["excluded_cities"] = "; ".join(sorted(city_panel.loc[small, "city_name"].unique()))
    write_csv(table, OUT / "Table_E5.csv")
    print("Table E5 done", flush=True)

    # E.5 Water DMUs with more than 20% non-vegetated land
    new = rescore_grid(panel, "water", lambda f: f[NONVEGETATED].sum(axis=1).le(0.2))
    table = pd.DataFrame(compare(grid["water"], new, water, "meta_score", "non-vegetated share > 20%"))
    excluded = grid["water"]["dmu_id"].isin(new["dmu_id"]).eq(False)
    in_all_years = grid["water"].groupby("grid_id")["year"].transform("size").eq(3)
    partly = grid["water"][in_all_years].assign(x=excluded).groupby("grid_id")["x"].agg(lambda s: 0 < s.sum() < 3)
    table["cells_excluded_in_some_years_only"] = int(partly.sum())
    centres = pd.DataFrame([dict(zip(("center_x_m", "center_y_m"), centre(new[new["year"].eq(y)], "meta_score")))
                            for y in (2010, 2020)])
    cities = add_center_context(centres)["center_city_name"]
    table["sde_centre_city_2010"], table["sde_centre_city_2020"] = cities[0], cities[1]
    write_csv(table, OUT / "Table_E6.csv")
    print(f"Table E6 done. Outputs in {OUT}")


if __name__ == "__main__":
    main()
