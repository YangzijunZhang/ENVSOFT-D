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
