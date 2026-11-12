"""Shared definitions used by every analysis script."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
YEARS = (2010, 2015, 2020)
RULES = ("dominant", "threshold_50", "threshold_70")
LAND_TYPES = {1: "Cropland", 2: "Forest", 4: "Grassland", 8: "Impervious"}
ALBERS_CRS = (
    "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 "
    "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)


@dataclass(frozen=True)
class System:
    name: str
    input: str
    good: str
    bad: str
    land_codes: tuple[int, ...]
    meta: bool = False

    @property
    def variables(self) -> tuple[str, str, str]:
        return self.input, self.good, self.bad


SYSTEMS = {
    "water": System(
        "water",
        "water_input_et_volume_m3_yr",
        "water_desirable_npp_t_c_yr",
        "water_undesirable_soil_erosion_t_yr",
        (1, 2, 4),
        True,
    ),
    "energy": System(
        "energy",
        "energy_input_electricity_kwh_yr",
        "energy_desirable_gdp_2017_int_usd_yr",
        "energy_undesirable_ghg_t_co2e_yr",
        (8,),
    ),
    "food": System(
        "food",
        "food_input_nitrogen_kg_n_yr",
        "food_desirable_crop_production_kt_yr",
        "food_undesirable_n2o_t_n2o_yr",
        (1,),
    ),
}
INDICATORS = tuple(column for system in SYSTEMS.values() for column in system.variables)
PRIMARY_SCORE = {"water": "meta_score", "energy": "group_score", "food": "group_score"}
RULE_COLUMNS = {
    "dominant": ("lulc_dominant_assigned", "lulc_dominant_code", "lulc_dominant_land_type"),
    "threshold_50": (
        "lulc_threshold_50_assigned",
        "lulc_threshold_50_code",
        "lulc_threshold_50_land_type",
    ),
    "threshold_70": (
        "lulc_threshold_70_assigned",
        "lulc_threshold_70_code",
        "lulc_threshold_70_land_type",
    ),
}
REGIONS = (
    ("A", "Hetao Irrigation District", 120_000, 380_000, 4_200_000, 4_450_000),
    ("B", "Western Loess Plateau", -110_000, 190_000, 3_680_000, 3_880_000),
    ("C", "Guanzhong Plain", 190_000, 470_000, 3_600_000, 3_720_000),
    ("D", "Lower Yellow River Plain", 670_000, 1_080_000, 3_680_000, 4_020_000),
)


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, na_values=["NA", ""], keep_default_na=True, **kwargs)


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        index=False,
        encoding="utf-8",
        na_rep="NA",
        float_format="%.12g",
        lineterminator="\n",
    )


def common_city_score_sample(
    frame: pd.DataFrame,
    score_column: str = "city_score",
    years: tuple[int, ...] = YEARS,
) -> pd.DataFrame:
    """Return the same successfully scored cities in every benchmark year.

    City-scale comparisons must not allow a city to enter one benchmark year
    but not another.  This helper retains the intersection of city IDs with a
    finite, positive score in every requested year.  When solver-success
    information is present, unsuccessful rows are excluded before the
    intersection is formed.

    The returned frame contains exactly one row per retained city and year.
    The original frame is not modified.
    """
    required = {"city_id", "year", score_column}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"City score table is missing: {missing}")

    sample = frame.copy()
    sample["city_id"] = sample["city_id"].astype(str)
    sample["year"] = pd.to_numeric(sample["year"], errors="coerce")
    if sample.duplicated(["city_id", "year"]).any():
        raise ValueError("City score table must contain one row per city and year")

    scores = pd.to_numeric(sample[score_column], errors="coerce")
    valid = (
        sample["year"].isin(years)
        & np.isfinite(scores.to_numpy(float))
        & scores.gt(0)
    )
    if "city_solver_success" in sample.columns:
        success = sample["city_solver_success"]
        if success.dtype != bool:
            success = success.astype(str).str.strip().str.lower().map(
                {"true": True, "false": False}
            )
        valid &= success.fillna(False).astype(bool)

    available = sample.loc[valid, ["city_id", "year"]]
    required_years = set(years)
    common_ids = [
        city_id
        for city_id, group in available.groupby("city_id", observed=True)
        if set(group["year"].astype(int)) == required_years
    ]
    if not common_ids:
        raise ValueError("No city has a valid score in every benchmark year")

    result = sample.loc[valid & sample["city_id"].isin(common_ids)].copy()
    counts = result.groupby("year", observed=True)["city_id"].nunique()
    expected = len(common_ids)
    if set(counts.index.astype(int)) != required_years or not counts.eq(expected).all():
        raise ValueError("The city comparison sample is not balanced across years")

    sort_columns = [column for column in ("year", "city_index", "city_id") if column in result]
    return result.sort_values(sort_columns, kind="stable").reset_index(drop=True)


def load_grid() -> pd.DataFrame:
    frame = read_csv(ROOT / "data/grid_panel_10km.csv", low_memory=False)
    required = {
        "dmu_id",
        "grid_id",
        "year",
        "grid_row",
        "grid_col",
        "x_center_m",
        "y_center_m",
        *INDICATORS,
        *(column for columns in RULE_COLUMNS.values() for column in columns),
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"grid_panel_10km.csv is missing: {missing}")
    if frame.duplicated(["grid_id", "year"]).any() or set(frame["year"]) != set(YEARS):
        raise ValueError("The grid input must be a unique 2010/2015/2020 panel")
    return frame


def grid_candidates(grid: pd.DataFrame, rule: str, system_name: str) -> pd.DataFrame:
    system = SYSTEMS[system_name]
    assigned_column, code_column, type_column = RULE_COLUMNS[rule]
    assigned = grid[assigned_column].fillna(False).astype(bool)
    frame = grid[assigned & grid[code_column].isin(system.land_codes)].copy()
    frame["lulc_rule"] = rule
    frame["lulc_code"] = pd.to_numeric(frame[code_column], errors="coerce").astype("Int64")
    frame["land_type"] = frame[type_column].astype("string")
    frame["scenario_dmu_id"] = frame["dmu_id"].astype(str) + "__" + rule
    values = frame[list(system.variables)].apply(pd.to_numeric, errors="coerce")
    frame["analysis_eligible"] = (
        np.isfinite(values.to_numpy(float)).all(axis=1)
        & values[system.input].gt(0)
        & values[system.good].gt(0)
        & values[system.bad].ge(0)
    )
    return frame


def describe(values: pd.Series) -> dict[str, float | int]:
    values = pd.to_numeric(values, errors="coerce").dropna()
    return {
        "n": len(values),
        "mean": values.mean(),
        "q25": values.quantile(0.25),
        "median": values.median(),
        "q75": values.quantile(0.75),
        "minimum": values.min(),
        "maximum": values.max(),
    }
