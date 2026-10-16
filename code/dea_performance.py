"""Estimate grid- and city-scale performance with the manuscript VRS-DDF."""

import argparse
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.spatial import ConvexHull, QhullError

from common import LAND_TYPES, ROOT, SYSTEMS, YEARS, System, grid_candidates, load_grid, read_csv, write_csv


@dataclass(frozen=True)
class Reference:
    x: np.ndarray
    y: np.ndarray
    b: np.ndarray
    ids: np.ndarray
    full_x: np.ndarray
    full_y: np.ndarray
    full_b: np.ndarray
    full_ids: np.ndarray
    n_original: int
    n_unique: int

    def full(self):
        return replace(
            self,
            x=self.full_x,
            y=self.full_y,
            b=self.full_b,
            ids=self.full_ids,
        )


def _hull_vertices(values: np.ndarray) -> np.ndarray:
    """Remove only points that cannot change the exact convex technology."""
    if len(values) <= values.shape[1] + 1:
        return np.arange(len(values))
    centred = values - values[0]
    rank = int(np.linalg.matrix_rank(centred, tol=1e-11))
    if rank == 0:
        return np.array([0])
    if rank == 1:
        dimension = int(np.argmax(np.ptp(values, axis=0)))
        return np.unique([np.argmin(values[:, dimension]), np.argmax(values[:, dimension])])
    try:
        return np.sort(ConvexHull(values, qhull_options="QJ").vertices)
    except QhullError:
        return np.arange(len(values))


def build_reference(frame: pd.DataFrame, system: System) -> Reference:
    x = frame[[system.input]].to_numpy(float)
    y = frame[[system.good]].to_numpy(float)
    b = frame[[system.bad]].to_numpy(float)
    ids = frame["scenario_dmu_id"].astype(str).to_numpy(object)
    physical_values = np.hstack((x, y, b))
    _, unique_index = np.unique(physical_values, axis=0, return_index=True)
    unique_index = np.sort(unique_index)
    ux, uy, ub, uids = x[unique_index], y[unique_index], b[unique_index], ids[unique_index]
    vertices = _hull_vertices(physical_values[unique_index])
    return Reference(
        ux[vertices],
        uy[vertices],
        ub[vertices],
        uids[vertices],
        ux,
        uy,
        ub,
        uids,
        len(frame),
        len(unique_index),
    )


def _solve_once(reference: Reference, x: float, y: float, b: float, method: str):
    """Observed direction g=(x_o,y_o,b_o), VRS, weak bad-output disposability."""
    n = len(reference.ids)
    objective = np.r_[np.zeros(n), -1.0]
    inequalities = np.vstack(
        (
            np.r_[reference.x[:, 0], x],
            np.r_[-reference.y[:, 0], y],
        )
    )
    equalities = np.vstack(
        (
            np.r_[reference.b[:, 0], b],
            np.r_[np.ones(n), 0.0],
        )
    )
    solved = linprog(
        objective,
        A_ub=inequalities,
        b_ub=np.array([x, -y]),
        A_eq=equalities,
        b_eq=np.array([b, 1.0]),
        bounds=[(0, None)] * n + [(0, 1)],
        method=method,
        options={"presolve": True},
    )
    if not solved.success or solved.x is None:
        return np.nan, False, int(solved.status), np.nan
    beta = float(np.clip(solved.x[-1], 0, 1))
    lambdas = solved.x[:-1]
    relative_residual = max(
        abs(lambdas.sum() - 1),
        max(float(lambdas @ reference.x[:, 0] - (1 - beta) * x), 0)
        / max(abs((1 - beta) * x), 1),
        max(float((1 + beta) * y - lambdas @ reference.y[:, 0]), 0)
        / max(abs((1 + beta) * y), 1),
        abs(float(lambdas @ reference.b[:, 0] - (1 - beta) * b))
        / max(abs(b), 1),
    )
    return beta, True, int(solved.status), relative_residual


def solve(reference: Reference, x: float, y: float, b: float) -> dict[str, object]:
    beta, success, status, residual = _solve_once(reference, x, y, b, "highs")
    if not success:
        beta, success, status, residual = _solve_once(reference.full(), x, y, b, "highs-ds")
    return {
        "beta": beta,
        "score": 1 / (1 + beta) if success else np.nan,
        "solver_success": success,
        "solver_status": status,
        "maximum_relative_residual": residual,
    }


def _score_rows(
    result: pd.DataFrame,
    indices: pd.Index,
    reference: Reference,
    system: System,
    prefix: str,
) -> None:
    for index in indices:
        row = result.loc[index]
        solved = solve(
            reference,
            float(row[system.input]),
            float(row[system.good]),
            float(row[system.bad]),
        )
        for name, value in solved.items():
            result.at[index, f"{prefix}_{name}"] = value


def score_grid_frontiers(
    candidates: pd.DataFrame, system: System, minimum_dmus: int = 9
) -> tuple[pd.DataFrame, pd.DataFrame]:
    result = candidates.copy()
    for prefix in ("group", "meta"):
        for column in ("beta", "score", "solver_status", "maximum_relative_residual"):
            result[f"{prefix}_{column}"] = np.nan
        result[f"{prefix}_solver_success"] = False

    summaries, estimable = [], []
    eligible = result["analysis_eligible"].fillna(False)
    for land_code in system.land_codes:
        indices = result.index[eligible & result["lulc_code"].eq(land_code)]
        if len(indices) < minimum_dmus:
            summaries.append(
                {
                    "system": system.name,
                    "frontier": "group",
                    "land_type": LAND_TYPES[land_code],
                    "n_reference": len(indices),
                    "status": "not_estimated",
                }
            )
            continue
        reference = build_reference(result.loc[indices], system)
        _score_rows(result, indices, reference, system, "group")
        estimable.append(land_code)
        summaries.append(
            {
                "system": system.name,
                "frontier": "group",
                "land_type": LAND_TYPES[land_code],
                "n_reference": reference.n_original,
                "n_unique": reference.n_unique,
                "n_hull": len(reference.ids),
                "status": "estimated",
            }
        )

    if system.meta and len(estimable) >= 2:
        indices = result.index[eligible & result["lulc_code"].isin(estimable)]
        reference = build_reference(result.loc[indices], system)
        _score_rows(result, indices, reference, system, "meta")
        summaries.append(
            {
                "system": system.name,
                "frontier": "meta",
                "land_type": ",".join(LAND_TYPES[code] for code in estimable),
                "n_reference": reference.n_original,
                "n_unique": reference.n_unique,
                "n_hull": len(reference.ids),
                "status": "estimated",
            }
        )
    return result, pd.DataFrame(summaries)


def score_city_frontier(frame: pd.DataFrame, system: System) -> tuple[pd.DataFrame, dict]:
    values = frame[list(system.variables)].apply(pd.to_numeric, errors="coerce")
    eligible = (
        np.isfinite(values.to_numpy(float)).all(axis=1)
        & values[system.input].gt(0)
        & values[system.good].ge(0)
        & values[system.bad].ge(0)
    )
    result = frame.loc[eligible].copy()
    if len(result) != 76 * len(YEARS):
        raise ValueError(f"{system.name}: expected 228 eligible city-years")
    reference = build_reference(result, system)
    for column in ("city_beta", "city_score", "city_solver_status", "city_maximum_relative_residual"):
        result[column] = np.nan
    result["city_solver_success"] = False
    _score_rows(result, result.index, reference, system, "city")
    return result, {
        "system": system.name,
        "frontier": "city",
        "n_reference": reference.n_original,
        "n_unique": reference.n_unique,
        "n_hull": len(reference.ids),
    }


def run_grid() -> None:
    grid = load_grid()
    summaries = []
    for name, system in SYSTEMS.items():
        scores, summary = score_grid_frontiers(grid_candidates(grid, "dominant", name), system)
        write_csv(scores, ROOT / f"results/performance/grid_{name}.csv")
        summaries.append(summary)
    write_csv(pd.concat(summaries, ignore_index=True), ROOT / "results/performance/grid_references.csv")


def run_city() -> None:
    city = read_csv(
        ROOT / "data/city_panel_10km.csv",
        dtype={"city_id": str, "province_id": str},
    )
    summaries = []
    for name, system in SYSTEMS.items():
        scores, summary = score_city_frontier(city, system)
        scores.insert(0, "system", name)
        write_csv(scores, ROOT / f"results/performance/city_{name}.csv")
        summaries.append(summary)
    write_csv(pd.DataFrame(summaries), ROOT / "results/performance/city_references.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--scale", choices=("grid", "city", "both"), default="both")
    scale = parser.parse_args().scale
    if scale in ("grid", "both"):
        run_grid()
    if scale in ("city", "both"):
        run_city()
