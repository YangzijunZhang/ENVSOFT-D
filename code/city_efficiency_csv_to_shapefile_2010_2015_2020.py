#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert the ACTUAL city_water.csv, city_energy.csv and city_food.csv produced
by dea_performance.py to YRB city polygon Shapefiles (2010, 2015, 2020).

SOURCE OF SCORES (NO DEA RECALCULATION)
---------------------------------------
    results/performance/city_water.csv
    results/performance/city_energy.csv
    results/performance/city_food.csv

The city DEA results contain 1 row per city/year/system and the fields:
city_id, year, city_score, city_beta, city_solver_success,
city_common_input_sample and city_common_score_sample. City polygons are
joined by city_id = ct_adcode, NEVER by city names or record order.

Outputs (by default: results/performance/city_shapefiles):
  water/water_city_score_2010.shp      (and 2015, 2020)
  energy/energy_city_score_2010.shp    (and 2015, 2020)
  food/food_city_score_2010.shp        (and 2015, 2020)
  combined/city_wef_scores_2010_2015_2020.shp
  export_quality_control.csv, unmatched_city_ids.csv,
  field_dictionary.csv, run_metadata.json

The 9 individual layers preserve all columns of the original CSVs using
explicit ESRI Shapefile-safe field aliases. Each contains EffScore (the
original city_score, unchanged) and W2010/E2015/F2020-style aliases.
The combined layer has 9 score fields, as well as 9 beta fields.

Run from the project root:
    python code/city_efficiency_csv_to_shapefile_2010_2015_2020.py

Existing output shapefiles are overwritten by default (they are generated
products). Pass --no-overwrite to block overwriting.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

# May help when using QGIS's embedded Python on macOS.
if "QGIS.app" in sys.executable:
    os.environ.setdefault("PROJ_LIB", "/Applications/QGIS.app/Contents/Resources/proj")
    os.environ.setdefault("GDAL_DATA", "/Applications/QGIS.app/Contents/Resources/gdal")

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.ops import unary_union

try:
    from common import ROOT as PROJECT_ROOT
except ImportError:
    PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROJECT_ROOT = Path(PROJECT_ROOT)

OLD_GIS_ROOT = Path("/Volumes/ZX6 1TB/paper1")
YEARS = (2010, 2015, 2020)
SYSTEMS = ("water", "energy", "food")
PREFIX = {"water": "W", "energy": "E", "food": "F"}
SHAPEFILE_EXTS = (".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix", ".fix", ".sbn", ".sbx")

# These 36 columns match the three actual uploaded CSVs, not the old
# city_system_scores_all.csv. Every alias is unique and at most 10 characters.
# All scientific values are carried through from the CSV without modification.
CSV_TO_SHP = {
    "system": "System",
    "city_id": "CityID",
    "water_input_et_volume_m3_yr": "W_ET",
    "water_input_et_volume_m3_yr__coverage_fraction": "W_ETCov",
    "water_desirable_npp_t_c_yr": "W_NPP",
    "water_desirable_npp_t_c_yr__coverage_fraction": "W_NPPCov",
    "water_undesirable_soil_erosion_t_yr": "W_Erosion",
    "water_undesirable_soil_erosion_t_yr__coverage_fraction": "W_EroCov",
    "energy_input_electricity_kwh_yr": "E_Elec",
    "energy_input_electricity_kwh_yr__coverage_fraction": "E_EleCov",
    "energy_desirable_gdp_2017_int_usd_yr": "E_GDP",
    "energy_desirable_gdp_2017_int_usd_yr__coverage_fraction": "E_GDPCov",
    "energy_undesirable_ghg_t_co2e_yr": "E_GHG",
    "energy_undesirable_ghg_t_co2e_yr__coverage_fraction": "E_GHGCov",
    "food_input_nitrogen_kg_n_yr": "F_Nitrogen",
    "food_input_nitrogen_kg_n_yr__coverage_fraction": "F_NCov",
    "food_desirable_crop_production_kt_yr": "F_Crop",
    "food_desirable_crop_production_kt_yr__coverage_fraction": "F_CropCov",
    "food_undesirable_n2o_t_n2o_yr": "F_N2O",
    "food_undesirable_n2o_t_n2o_yr__coverage_fraction": "F_N2OCov",
    "year": "Year",
    "city_index": "SrcIndex",
    "city_name": "SrcName",
    "province_id": "SrcProvID",
    "province_name": "SrcProvNm",
    "city_area_inside_yrb_km2": "SrcArea",
    "x_center_m": "XCenter",
    "y_center_m": "YCenter",
    "scenario_dmu_id": "DMUID",
    "city_common_input_sample": "CommonInp",
    "city_beta": "Beta",
    "city_score": "EffScore",
    "city_solver_status": "SolvStat",
    "city_maximum_relative_residual": "MaxResid",
    "city_solver_success": "Solved",
    "city_common_score_sample": "CommonScr",
}

CORE_REQUIRED = {
    "system", "city_id", "year", "city_score", "city_beta",
    "city_common_input_sample", "city_common_score_sample",
    "city_solver_success", "scenario_dmu_id",
}
BOOL_COLS = ("city_common_input_sample", "city_common_score_sample", "city_solver_success")


def default_file(*relative_paths: str) -> Path | None:
    candidates = [root / rel for rel in relative_paths
                  for root in (PROJECT_ROOT, OLD_GIS_ROOT)]
    return next((path for path in candidates if path.is_file()), candidates[0] if candidates else None)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Convert actual DEA city CSVs to YRB polygon Shapefiles")
    for system in SYSTEMS:
        p.add_argument(f"--{system}-csv", type=Path,
                       default=PROJECT_ROOT / "results/performance" / f"city_{system}.csv",
                       help=f"Direct path to city_{system}.csv")
    p.add_argument("--city-shp", type=Path,
                   default=default_file("rawdata/adBoundary/YRB/YRB_City.shp"))
    p.add_argument("--basin-shp", type=Path,
                   default=default_file("rawdata/adBoundary/YRB/YRB.shp"))
    p.add_argument("--reference-raster", type=Path,
                   default=(next((root / "inputdata/grids/10km/YRB_boundary_fraction_10km.tif"
                                  for root in (PROJECT_ROOT, OLD_GIS_ROOT)
                                  if (root / "inputdata/grids/10km/YRB_boundary_fraction_10km.tif").is_file()), None)),
                   help="Optional; use the raster CRS as Shapefile output CRS")
    p.add_argument("--output-root", type=Path,
                   default=PROJECT_ROOT / "results/performance/city_shapefiles")
    p.add_argument("--years", type=int, nargs="+", choices=YEARS, default=list(YEARS))
    p.add_argument("--systems", nargs="+", choices=SYSTEMS, default=list(SYSTEMS))
    p.add_argument("--keep-unmatched-boundary-cities", action="store_true",
                   help="Also export boundary polygons without source CSV results (null scores)")
    p.add_argument("--include-partial-scores", action="store_true",
                   help="Allow successfully solved rows not in the all-years common-score sample")
    p.add_argument("--no-clip", action="store_true",
                   help="Do not intersect city polygons with the YRB polygon")
    p.add_argument("--no-combined", action="store_true")
    p.add_argument("--no-overwrite", action="store_true")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()
    if len(args.years) != len(set(args.years)) or len(args.systems) != len(set(args.systems)):
        p.error("Duplicate --years or --systems")
    return args


def require_files(paths) -> None:
    missing = [str(path) for path in paths if path is None or not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing input file(s):\n  " + "\n  ".join(missing)
                                + "\nSpecify a path using --water-csv / --city-shp / --basin-shp, etc.")


def normalize_id(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.replace(r"\.0$", "", regex=True)


def parse_bool(series: pd.Series, name: str, filename: Path) -> pd.Series:
    mapped = series.astype("string").str.strip().str.lower().map({
        "true": True, "false": False, "1": True, "0": False,
        "yes": True, "no": False,
    })
    if mapped.isna().any():
        raise ValueError(f"{filename.name}: invalid or missing values in {name}")
    return mapped.astype(bool)


def load_city_csv(path: Path, system: str, years: tuple[int, ...]) -> pd.DataFrame:
    # Read the exact CSV directly. No dependency on old output files or on
    # recomputing the DEA model.
    df = pd.read_csv(path, dtype={"city_id": "string", "province_id": "string"},
                     low_memory=False)
    df.columns = [str(col).strip() for col in df.columns]
    missing = CORE_REQUIRED - set(df.columns)
    if missing:
        raise KeyError(f"{path.name} is missing required DEA columns: {sorted(missing)}")
    unknown = sorted(set(df.columns) - set(CSV_TO_SHP))
    if unknown:
        raise KeyError(f"{path.name} contains unrecognised columns. Add explicit short aliases "
                       f"to CSV_TO_SHP rather than silently dropping them: {unknown}")
    df["city_id"] = normalize_id(df["city_id"])
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    if df["city_id"].isna().any() or df["city_id"].eq("").any() or df["year"].isna().any():
        raise ValueError(f"{path.name}: city_id / year is blank or invalid")
    if not df["system"].astype("string").str.strip().str.lower().eq(system).all():
        raise ValueError(f"{path.name}: system column does not match {system}")
    df["system"] = system
    if df.duplicated(["city_id", "year"]).any():
        raise ValueError(f"{path.name}: duplicate city_id/year rows; no averaging is permitted")
    if df.duplicated(["scenario_dmu_id"]).any():
        raise ValueError(f"{path.name}: duplicate scenario_dmu_id rows")

    for field in BOOL_COLS:
        df[field] = parse_bool(df[field], field, path)
    for field in ("city_score", "city_beta"):
        df[field] = pd.to_numeric(df[field], errors="coerce")

    if (df["city_common_score_sample"] & ~df["city_common_input_sample"]).any():
        raise ValueError(f"{path.name}: common score sample not in common input sample")
    if df.groupby("city_id")["city_common_score_sample"].nunique().gt(1).any():
        raise ValueError(f"{path.name}: inconsistent all-years sample flag across years")

    solved = df["city_solver_success"]
    invalid = solved & (
        ~np.isfinite(df["city_score"]) | ~np.isfinite(df["city_beta"])
        | ~df["city_score"].between(0.5 - 1e-7, 1.0 + 1e-7)
        | ~df["city_beta"].between(-1e-7, 1.0 + 1e-7)
        | ~pd.Series(np.isclose(df["city_score"], 1 / (1 + df["city_beta"]),
                                rtol=1e-7, atol=1e-7), index=df.index)
    )
    if invalid.any():
        examples = df.loc[invalid, ["city_id", "year", "city_score", "city_beta"]].head(5)
        raise ValueError(f"{path.name}: invalid solved DEA values:\n{examples.to_string(index=False)}")

    # Filter only requested years; values come straight from the corresponding
    # source rows, not from unrelated years or another system.
    df = df.loc[df["year"].isin(years)].copy()
    for year in years:
        if not df["year"].eq(year).any():
            raise ValueError(f"{path.name}: no source rows for {year}")
    logging.info("Read %s directly: %d rows; %s", path.name, len(df),
                 ", ".join(f"{year}: {int(df.year.eq(year).sum())} cities" for year in years))
    return df


def make_valid_polygon(geom):
    if geom is None or geom.is_empty:
        return None
    if not geom.is_valid:
        try:
            from shapely import make_valid
        except ImportError:
            geom = geom.buffer(0)
        else:
            geom = make_valid(geom)
    if geom is None or geom.is_empty:
        return None
    if geom.geom_type in ("Polygon", "MultiPolygon"):
        return geom
    fragments = []
    for item in getattr(geom, "geoms", []):
        if item.geom_type in ("Polygon", "MultiPolygon") and not item.is_empty:
            fragments.append(item)
    return unary_union(fragments) if fragments else None


def load_city_boundaries(city_shp: Path, basin_shp: Path, reference_raster: Path | None,
                         no_clip: bool) -> gpd.GeoDataFrame:
    cities = gpd.read_file(city_shp)
    if cities.crs is None or cities.empty:
        raise ValueError(f"City boundary is empty or has no CRS: {city_shp}")
    required = {"ct_adcode", "ct_name", "pr_adcode", "pr_name", "geometry"}
    missing = required - set(cities.columns)
    if missing:
        raise KeyError(f"Missing required city boundary columns: {sorted(missing)}")
    source_crs = cities.crs
    cities = cities[["ct_adcode", "ct_name", "pr_adcode", "pr_name", "geometry"]].copy()
    cities["CityID"] = normalize_id(cities["ct_adcode"])
    cities["CityName"] = cities["ct_name"].astype("string").fillna("")
    cities["ProvID"] = normalize_id(cities["pr_adcode"])
    cities["ProvName"] = cities["pr_name"].astype("string").fillna("")
    if cities["CityID"].isna().any() or cities["CityID"].eq("").any():
        raise ValueError("ct_adcode contains blank IDs")

    cities = cities[["CityID", "CityName", "ProvID", "ProvName", "geometry"]]
    cities["geometry"] = cities.geometry.map(make_valid_polygon)
    cities = gpd.GeoDataFrame(cities, geometry="geometry", crs=source_crs)
    cities = cities.loc[cities.geometry.notna() & ~cities.geometry.is_empty].copy()

    if reference_raster is not None:
        import rasterio
        with rasterio.open(reference_raster) as raster:
            if raster.crs is None:
                raise ValueError(f"Reference raster has no CRS: {reference_raster}")
            target_crs = raster.crs
    else:
        target_crs = source_crs
    cities = cities.to_crs(target_crs)

    if not no_clip:
        basin = gpd.read_file(basin_shp)
        if basin.crs is None or basin.empty:
            raise ValueError(f"Basin boundary is empty or has no CRS: {basin_shp}")
        basin = basin.to_crs(target_crs)
        basin_geom = make_valid_polygon(unary_union(basin.geometry.dropna()))
        if basin_geom is None:
            raise ValueError("Basin contains no valid polygon geometry")
        cities["geometry"] = cities.geometry.map(
            lambda geom: make_valid_polygon(geom.intersection(basin_geom)))
        cities = gpd.GeoDataFrame(cities, geometry="geometry", crs=target_crs)
        cities = cities.loc[cities.geometry.notna() & ~cities.geometry.is_empty].copy()

    if cities["CityID"].duplicated().any():
        cities = cities.dissolve(by="CityID", as_index=False,
                                aggfunc={"CityName": "first", "ProvID": "first", "ProvName": "first"})
    cities = cities[["CityID", "CityName", "ProvID", "ProvName", "geometry"]]
    cities = cities.sort_values("CityID").reset_index(drop=True)
    if cities["CityID"].duplicated().any():
        raise ValueError("CityID is not unique in prepared boundary")
    logging.info("Prepared %d city polygons (clip=%s, CRS=%s)",
                 len(cities), not no_clip, cities.crs)
    return cities


def shapefile_family(path: Path) -> list[Path]:
    return [path.with_suffix(ext) for ext in SHAPEFILE_EXTS]


def write_shape(gdf: gpd.GeoDataFrame, path: Path, overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    existing = [item for item in shapefile_family(path) if item.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"Output exists: {path} (omit --no-overwrite to replace it)")
    for item in existing:
        item.unlink()
    invalid_fields = [field for field in gdf.columns if field != "geometry" and len(field) > 10]
    if invalid_fields:
        raise ValueError(f"Shapefile fields exceed 10 characters: {invalid_fields}")
    gdf.to_file(path, driver="ESRI Shapefile", encoding="UTF-8", index=False)


def mapped_rows(df: pd.DataFrame, include_partial: bool) -> pd.DataFrame:
    df = df.copy()
    good = df["city_common_input_sample"] & df["city_solver_success"]
    if not include_partial:
        good &= df["city_common_score_sample"]
    # A failed/non-common city remains a polygon but has no mapped score.
    # The original CSV is NOT modified or re-estimated.
    df.loc[~good, ["city_score", "city_beta"]] = np.nan
    return df


def source_as_shapefile_fields(df: pd.DataFrame) -> pd.DataFrame:
    # Preserve all original columns, but create intentional <=10-character names.
    result = df.rename(columns=CSV_TO_SHP).copy()
    for col in ("CommonInp", "CommonScr", "Solved"):
        result[col] = result[col].astype("int16")
    result["Year"] = result["Year"].astype("int32")
    result["CityID"] = normalize_id(result["CityID"])
    result["SrcProvID"] = normalize_id(result["SrcProvID"])
    return result


def alias(system: str, year: int) -> str:
    return f"{PREFIX[system]}{year}"


def export_all(boundaries: gpd.GeoDataFrame, sources: dict[str, pd.DataFrame],
               out_root: Path, years: tuple[int, ...], systems: tuple[str, ...],
               keep_extra: bool, include_partial: bool, overwrite: bool,
               no_combined: bool):
    qc = []
    unmatched = []
    score_tables = {}
    shape_ids = set(boundaries["CityID"].astype(str))
    for system in systems:
        for year in years:
            csv_rows = sources[system].loc[sources[system]["year"].eq(year)].copy()
            csv_ids = set(csv_rows["city_id"].astype(str))
            missing_shape = sorted(csv_ids - shape_ids)
            if missing_shape:
                unmatched.extend({"system": system, "year": year, "city_id": city_id,
                                  "reason": "CSV_city_ID_missing_from_clipped_boundary"}
                                 for city_id in missing_shape)
                raise ValueError(f"{system} {year}: {len(missing_shape)} source CSV city IDs "
                                 f"have no YRB polygon: {missing_shape[:20]}. "
                                 "Check ct_adcode and boundary clipping before exporting.")
            extras = sorted(shape_ids - csv_ids)
            unmatched.extend({"system": system, "year": year, "city_id": city_id,
                              "reason": "Boundary_city_ID_has_no_CSV_row"}
                             for city_id in extras)

            cleaned = mapped_rows(csv_rows, include_partial)
            attrs = source_as_shapefile_fields(cleaned)
            field = alias(system, year)
            attrs[field] = attrs["EffScore"]
            attrs["Matched"] = np.int16(1)
            # By default export exactly the CSV city sample, not unrelated cities.
            boundary_set = boundaries if keep_extra else boundaries.loc[
                boundaries["CityID"].isin(csv_ids)].copy()
            merged = boundary_set.merge(attrs, on="CityID", how="left", validate="one_to_one")
            merged = gpd.GeoDataFrame(merged, geometry="geometry", crs=boundaries.crs)
            if keep_extra:
                merged["Matched"] = merged["Matched"].fillna(0).astype("int16")
                merged["Year"] = merged["Year"].fillna(year).astype("int32")
                merged["System"] = merged["System"].fillna(system)
            if int(merged["Matched"].sum()) != len(csv_rows):
                raise ValueError(f"{system} {year}: join has lost or duplicated source records")
            path = out_root / system / f"{system}_city_score_{year}.shp"
            write_shape(merged, path, overwrite)

            # Read back DBF attributes and verify every written score against
            # the corresponding CSV row. Compare by CityID, never by row order.
            reloaded = gpd.read_file(path)[["CityID", "EffScore", field]]
            original = cleaned[["city_id", "city_score"]].rename(
                columns={"city_id": "CityID", "city_score": "CSVScore"})
            check = original.merge(reloaded, on="CityID", how="left", validate="one_to_one")
            for out_col in ("EffScore", field):
                expected = pd.to_numeric(check["CSVScore"], errors="coerce").to_numpy(float)
                actual = pd.to_numeric(check[out_col], errors="coerce").to_numpy(float)
                if not np.allclose(expected, actual, atol=1e-9, rtol=1e-9, equal_nan=True):
                    raise ValueError(f"{system} {year}: output {out_col} differs from source CSV!")

            values = pd.to_numeric(merged["EffScore"], errors="coerce").dropna()
            qc.append({
                "system": system, "year": year, "csv_file": str(sources[system].attrs["source_path"]),
                "shape_file": str(path), "n_source_rows": len(csv_rows),
                "n_exported_polygons": len(merged), "n_missing_polygon_ids": len(missing_shape),
                "n_boundary_only_ids": len(extras), "n_mapped_scores": len(values),
                "n_null_scores": int(merged["EffScore"].isna().sum()),
                "n_common_score": int(cleaned["city_common_score_sample"].sum()),
                "n_solver_success": int(cleaned["city_solver_success"].sum()),
                "score_min": values.min() if len(values) else np.nan,
                "score_mean": values.mean() if len(values) else np.nan,
                "score_median": values.median() if len(values) else np.nan,
                "score_max": values.max() if len(values) else np.nan,
                "verified_against_csv": True,
            })
            score_tables[(system, year)] = cleaned[["city_id", "city_score", "city_beta"]].copy()
            logging.info("Wrote %s | CSV rows=%d | score=%d | null=%d | CSV match=PASS",
                         path, len(csv_rows), len(values), int(merged["EffScore"].isna().sum()))

    combined_path = None
    if not no_combined:
        combined = boundaries.copy()
        if not keep_extra:
            sample_ids = set().union(*(set(df["city_id"].astype(str)) for df in sources.values()))
            combined = combined.loc[combined["CityID"].isin(sample_ids)].copy()
        for system in systems:
            for year in years:
                score_field = alias(system, year)
                beta_field = "B" + score_field
                source = score_tables[(system, year)].rename(columns={
                    "city_id": "CityID", "city_score": score_field, "city_beta": beta_field})
                combined = combined.merge(source, on="CityID", how="left", validate="one_to_one")
        combined = gpd.GeoDataFrame(combined, geometry="geometry", crs=boundaries.crs)
        suffix = "_".join(map(str, years))
        combined_path = out_root / "combined" / f"city_wef_scores_{suffix}.shp"
        write_shape(combined, combined_path, overwrite)
        logging.info("Wrote combined layer %s", combined_path)

    return pd.DataFrame(qc), pd.DataFrame(
        unmatched, columns=["system", "year", "city_id", "reason"]), combined_path


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s | %(levelname)s | %(message)s", datefmt="%H:%M:%S")
    assert len(CSV_TO_SHP) == len(set(CSV_TO_SHP.values()))
    assert all(len(name) <= 10 for name in CSV_TO_SHP.values())
    years = tuple(args.years)
    systems = tuple(args.systems)
    csv_paths = {system: getattr(args, f"{system}_csv") for system in systems}
    require_files([*csv_paths.values(), args.city_shp]
                  + ([] if args.no_clip else [args.basin_shp])
                  + ([args.reference_raster] if args.reference_raster else []))

    sources = {system: load_city_csv(csv_paths[system], system, years) for system in systems}
    # The dataframe .attrs is a stable way to carry the file origin to QC.
    for system, frame in sources.items():
        frame.attrs["source_path"] = csv_paths[system]
    boundaries = load_city_boundaries(args.city_shp, args.basin_shp,
                                     args.reference_raster, args.no_clip)
    qc, unmatched, combined = export_all(
        boundaries, sources, args.output_root, years, systems,
        args.keep_unmatched_boundary_cities, args.include_partial_scores,
        not args.no_overwrite, args.no_combined,
    )

    args.output_root.mkdir(parents=True, exist_ok=True)
    qc.to_csv(args.output_root / "export_quality_control.csv", index=False, encoding="utf-8")
    unmatched.to_csv(args.output_root / "unmatched_city_ids.csv", index=False, encoding="utf-8")
    field_rows = [{"source_csv_field": key, "shapefile_field": value}
                  for key, value in CSV_TO_SHP.items()]
    field_rows += [
        {"source_csv_field": "boundary:ct_name", "shapefile_field": "CityName"},
        {"source_csv_field": "boundary:pr_adcode", "shapefile_field": "ProvID"},
        {"source_csv_field": "boundary:pr_name", "shapefile_field": "ProvName"},
        {"source_csv_field": "join_success_flag", "shapefile_field": "Matched"},
    ]
    for system in systems:
        for year in years:
            field_rows.append({"source_csv_field": f"{system} / {year}: city_score",
                               "shapefile_field": alias(system, year)})
            field_rows.append({"source_csv_field": f"{system} / {year}: city_beta (combined)",
                               "shapefile_field": "B" + alias(system, year)})
    pd.DataFrame(field_rows).to_csv(args.output_root / "field_dictionary.csv", index=False)
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_csvs": {k: str(v.resolve()) for k, v in csv_paths.items()},
        "city_boundary": str(args.city_shp.resolve()),
        "basin_boundary": None if args.no_clip else str(args.basin_shp.resolve()),
        "reference_raster": str(args.reference_raster) if args.reference_raster else None,
        "output_root": str(args.output_root.resolve()),
        "years": list(years), "systems": list(systems),
        "source_of_efficiency": "original CSV city_score; no DEA recomputation",
        "method": "city pooled VRS-DDF (per system); weak disposability",
        "score_definition": "city_score = 1/(1+city_beta)",
        "join_key": "source city_id = boundary ct_adcode",
        "sample_selection": "all-years common-score sample" if not args.include_partial_scores
                            else "successful city-years in common-input sample",
        "missing_city_handling": "fail if source city lacks polygon; other boundaries optional",
        "combined_shapefile": str(combined) if combined else None,
    }
    (args.output_root / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    logging.info("DONE: %d individual Shapefiles%s. Outputs: %s",
                 len(years) * len(systems), " + combined layer" if combined else "", args.output_root)


if __name__ == "__main__":
    main()
