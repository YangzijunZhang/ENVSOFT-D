#!/usr/bin/env python3
"""Build the repository 10 km panel from aligned TIFFs and city-grid weights.

Place this file in ENVSOFT-D-main/code/ and run from any working directory:

    python code/raster_10km_to_csv.py

Requirements: numpy, pandas, rasterio, pyproj.

Output (relative to the repository root):
    data/grid_panel_10km.csv

The 15-column city-grid weights are a *spatial input*, calculated from exact
city/basin/grid polygon intersections; they cannot be reconstructed from the
indicator TIFFs. By default this script reads and validates the existing
data/city_grid_weights_10km.csv as an input and does not modify that file.
Use --city-weights-source to supply the same crosswalk from another location.
Its rows provide the 8,622 in-basin cell locations and boundary fractions.
All 2010/2015/2020 LULC and nine indicator values come from TIFFs.

An existing grid_panel_10km.csv is compared against the TIFF-derived values.
It is never replaced if values differ unless --replace-existing is given; a
backup is created first. The script never rescales raster totals by boundary
fractions or allocation weights.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer


YEARS = (2010, 2015, 2020)
LULC = {
    1: "Cropland", 2: "Forest", 3: "Shrub", 4: "Grassland",
    5: "Water", 6: "Snow/Ice", 7: "Barren", 8: "Impervious", 9: "Wetland",
}
FRACTIONS = {
    1: "cropland", 2: "forest", 3: "shrub", 4: "grassland",
    5: "water", 6: "snow_ice", 7: "barren", 8: "impervious", 9: "wetland",
}
RULES = {
    "dominant": "LULC_dominant_{year}_10km.tif",
    "threshold_50": "LULC_dominant_threshold_50_{year}_10km.tif",
    "threshold_70": "LULC_dominant_threshold_70_{year}_10km.tif",
    "threshold_80": "LULC_dominant_threshold_80_{year}_10km.tif",
}
VARIABLES = {
    "water_input_et_volume_m3_yr": "water/ET_volume_{year}_10km.tif",
    "water_desirable_npp_t_c_yr": "water/NPP_total_{year}_10km.tif",
    "water_undesirable_soil_erosion_t_yr": "water/SoilErosion_total_{year}_10km.tif",
    "energy_input_electricity_kwh_yr": "energy/ElectricPowerConsumption_{year}_10km.tif",
    "energy_desirable_gdp_2017_int_usd_yr": "energy/GDP_total_{year}_10km.tif",
    "energy_undesirable_ghg_t_co2e_yr": "energy/GHG_total_{year}_10km.tif",
    "food_input_nitrogen_kg_n_yr": "food/NitrogenApplication_4crops_{year}_10km.tif",
    "food_desirable_crop_production_kt_yr": "food/CropProduction_4crops_{year}_10km.tif",
    "food_undesirable_n2o_t_n2o_yr": "food/N2O_3categories_{year}_10km.tif",
}
GRID_COLUMNS = (
    "dmu_id", "grid_id", "year", "resolution", "resolution_m", "grid_row",
    "grid_col", "x_center_m", "y_center_m", "longitude_wgs84",
    "latitude_wgs84", "boundary_fraction", "projected_cell_area_km2",
    "area_inside_yrb_km2",
)
WEIGHT_COLUMNS = (
    "city_index", "city_id", "city_name", "province_id", "province_name",
    "city_area_inside_yrb_m2", "city_area_inside_yrb_km2", "cell_id",
    "grid_row", "grid_col", "boundary_fraction", "full_cell_area_m2",
    "basin_cell_area_m2", "overlap_area_m2", "allocation_weight",
)
PANEL_COLUMNS = (*GRID_COLUMNS, "lulc_dominant_share", "lulc_valid_area_fraction", *(
    item for rule in RULES for item in (
        f"lulc_{rule}_code", f"lulc_{rule}_land_type", f"lulc_{rule}_assigned"
    )
), *(f"lulc_fraction_{band}_{name}" for band, name in FRACTIONS.items()),
    *VARIABLES)
NA_TEXT = "NA"


def arguments() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raster-root", type=Path, default=repo / "data" / "raster data")
    parser.add_argument("--city-weights-source", type=Path,
                        default=repo / "data" / "city_grid_weights_10km.csv",
                        help="Precomputed exact city/basin/grid polygon crosswalk (15 columns).")
    parser.add_argument("--output-dir", type=Path, default=repo / "data",
                        help="Directory for grid_panel_10km.csv (default: repository data/).")
    parser.add_argument("--replace-existing", action="store_true",
                        help="Back up and replace an existing CSV if its values differ.")
    return parser.parse_args()


def read_weights(path: Path) -> pd.DataFrame:
    """Read the exact spatial crosswalk; do not infer city geometry from TIFFs."""
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing city/grid intersection weights: {path}. "
            "The indicator TIFFs contain no city polygons or basin overlaps."
        )
    frame = pd.read_csv(path, dtype={"city_id": str, "province_id": str},
                        na_values=[NA_TEXT, ""], low_memory=False)
    if tuple(frame.columns) != WEIGHT_COLUMNS or frame.empty or frame.isna().any().any():
        raise ValueError("City weights must have exactly the 15 expected columns and no NA")
    if frame.duplicated(["city_id", "grid_row", "grid_col"]).any():
        raise ValueError("Repeated city/grid intersection in city weights")
    for col in ("grid_row", "grid_col", "cell_id", "city_index"):
        values = frame[col].to_numpy(float)
        if not (np.isfinite(values) & (values >= 0) & (values == np.rint(values))).all():
            raise ValueError(f"Invalid integer values in {col}")
    fraction = frame.boundary_fraction.to_numpy(float)
    full = frame.full_cell_area_m2.to_numpy(float)
    basin = frame.basin_cell_area_m2.to_numpy(float)
    overlap = frame.overlap_area_m2.to_numpy(float)
    alloc = frame.allocation_weight.to_numpy(float)
    if not (np.isfinite(fraction) & (fraction > 0) & (fraction <= 1.000001)).all():
        raise ValueError("boundary_fraction must be in (0, 1]")
    # The polygon intersection was calculated separately from the basin-cell
    # area. At some boundary cells its area exceeds the basin area by a few
    # square metres (~1e-8 relative), and the published weight is capped at 1.
    area_tolerance = np.maximum(1.0, basin * 1e-6)
    valid_overlap = (np.isfinite(overlap) & (overlap > 0)
                     & (overlap <= basin + area_tolerance))
    if not valid_overlap.all():
        i = int(np.flatnonzero(~valid_overlap)[0])
        raise ValueError(f"Invalid city/grid overlap at CSV row {i + 2}: "
                         f"overlap={overlap[i]}, basin area={basin[i]}")
    small_excess = (overlap > basin).sum()
    if small_excess:
        logging.info("Accepted %s sub-ppm polygon area differences at basin boundary cells",
                     small_excess)
    if not (np.isclose(basin, full * fraction, rtol=1e-8, atol=1)
            & np.isclose(alloc, np.minimum(1.0, overlap / basin),
                         rtol=1e-6, atol=1e-6)).all():
        raise ValueError("Boundary fractions, areas and allocation weights disagree")
    cells = frame.groupby(["grid_row", "grid_col"], sort=False)
    if (cells[["cell_id", "boundary_fraction", "full_cell_area_m2"]]
            .nunique().gt(1).any().any()):
        raise ValueError("A grid cell has conflicting spatial attributes")
    if cells.allocation_weight.sum().gt(1.000001).any():
        raise ValueError("City weights exceed the available area of a grid cell")
    city = frame.groupby("city_id", sort=False)
    if (city[["city_index", "city_name", "province_id", "province_name",
               "city_area_inside_yrb_m2", "city_area_inside_yrb_km2"]]
            .nunique().gt(1).any().any()):
        raise ValueError("City metadata vary among its intersected grid cells")
    city_area = city.city_area_inside_yrb_m2.first()
    if not np.allclose(city.overlap_area_m2.sum(), city_area, rtol=1e-8, atol=1):
        raise ValueError("City overlap areas do not sum to the reported city area")
    if not np.allclose(city.city_area_inside_yrb_km2.first() * 1e6,
                       city_area, rtol=1e-8, atol=1):
        raise ValueError("City km2 and m2 area columns disagree")
    return frame


class RasterReader:
    def __init__(self, root: Path):
        self.root = root
        first = root / "LULC" / "LULC_dominant_2010_10km.tif"
        if not first.is_file():
            raise FileNotFoundError(f"Missing reference raster: {first}")
        with rasterio.open(first) as src:
            self.shape, self.transform, self.crs = src.shape, src.transform, src.crs
        if self.crs is None or not (
            np.isclose(abs(self.transform.a), 10000, atol=0.01)
            and np.isclose(abs(self.transform.e), 10000, atol=0.01)
            and np.isclose(self.transform.b, 0)
            and np.isclose(self.transform.d, 0)
        ):
            raise ValueError(f"Reference TIFF is not an aligned 10 km projected grid: {first}")
        self.gdp_units: set[str] = set()

    def read(self, relative: str, bands: int = 1) -> np.ndarray:
        path = self.root / relative
        if not path.is_file():
            raise FileNotFoundError(f"Missing 10 km raster: {path}")
        with rasterio.open(path) as src:
            if src.count != bands or src.crs is None:
                raise ValueError(f"Expected {bands} band(s) and a valid CRS: {path}")
            if (src.shape != self.shape or src.crs != self.crs
                    or not np.allclose(tuple(src.transform), tuple(self.transform),
                                       rtol=0, atol=1e-8)):
                raise ValueError(f"Grid/CRS/transform mismatch; no automatic resampling: {path}")
            if bands == 9:
                for band in LULC:
                    desc = src.descriptions[band - 1] or ""
                    if f"class_{band}_" not in desc:
                        raise ValueError(f"LULC band {band} has unexpected description {desc!r}: {path}")
            values = src.read(masked=True).astype("float64").filled(np.nan)
            if relative.startswith("energy/GDP_total_"):
                self.gdp_units.add(src.tags().get("unit", ""))
            return values if bands > 1 else values[0]

    def validate_reference_geometry(self, frame: pd.DataFrame) -> None:
        row, col = frame.grid_row.to_numpy(int), frame.grid_col.to_numpy(int)
        if (row.min() < 0 or col.min() < 0 or row.max() >= self.shape[0]
                or col.max() >= self.shape[1]):
            raise ValueError("Reference panel has row/column indices outside the TIFF grid")
        x, y = rasterio.transform.xy(self.transform, row, col, offset="center")
        if not (np.allclose(x, frame.x_center_m, atol=1e-5, rtol=0)
                and np.allclose(y, frame.y_center_m, atol=1e-5, rtol=0)):
            raise ValueError("Reference panel grid centres do not match TIFF transform")
        lon, lat = Transformer.from_crs(self.crs, "EPSG:4326", always_xy=True).transform(x, y)
        if not (np.allclose(lon, frame.longitude_wgs84, atol=1e-6, rtol=0)
                and np.allclose(lat, frame.latitude_wgs84, atol=1e-6, rtol=0)):
            raise ValueError("Reference panel WGS84 coordinates do not match raster CRS")
        area = abs(self.transform.a * self.transform.e - self.transform.b * self.transform.d) / 1e6
        if not np.allclose(frame.projected_cell_area_km2, area, atol=1e-8, rtol=0):
            raise ValueError("Reference panel projected grid area does not match raster transform")
        expected_id = [f"10km_r{r:04d}_c{c:04d}" for r, c in zip(row, col)]
        if expected_id != frame.grid_id.to_list():
            raise ValueError("Reference grid IDs disagree with raster row/column indices")
        if [f"{g}_{y}" for g, y in zip(expected_id, frame.year)] != frame.dmu_id.to_list():
            raise ValueError("Reference DMU IDs disagree with grid IDs and years")


def make_grid_index(weights: pd.DataFrame, reader: RasterReader) -> pd.DataFrame:
    """Derive the stable raster index from the exact city-grid crosswalk."""
    cells = (weights[["grid_row", "grid_col", "cell_id", "boundary_fraction",
                      "full_cell_area_m2", "basin_cell_area_m2"]]
             .drop_duplicates(["grid_row", "grid_col"])
             .sort_values(["grid_row", "grid_col"], kind="stable")
             .reset_index(drop=True))
    row, col = cells.grid_row.to_numpy(int), cells.grid_col.to_numpy(int)
    height, width = reader.shape
    if ((row < 0) | (col < 0) | (row >= height) | (col >= width)).any():
        raise ValueError("City weight row/col indices exceed the raster extent")
    if not np.array_equal(cells.cell_id.to_numpy(int), row * width + col):
        raise ValueError("cell_id does not equal raster row * width + column")
    area_m2 = abs(reader.transform.a * reader.transform.e
                  - reader.transform.b * reader.transform.d)
    if not np.allclose(cells.full_cell_area_m2, area_m2, rtol=0, atol=1):
        raise ValueError("Weights and raster have different cell areas")
    x, y = rasterio.transform.xy(reader.transform, row, col, offset="center")
    lon, lat = Transformer.from_crs(reader.crs, "EPSG:4326", always_xy=True).transform(x, y)
    grid_id = [f"10km_r{r:04d}_c{c:04d}" for r, c in zip(row, col)]
    one_year = pd.DataFrame({
        "grid_id": grid_id, "resolution": "10km", "resolution_m": 10000,
        "grid_row": row, "grid_col": col, "x_center_m": x, "y_center_m": y,
        "longitude_wgs84": lon, "latitude_wgs84": lat,
        "boundary_fraction": cells.boundary_fraction.to_numpy(float),
        "projected_cell_area_km2": area_m2 / 1e6,
        "area_inside_yrb_km2": cells.basin_cell_area_m2.to_numpy(float) / 1e6,
    })
    frames = []
    for year in YEARS:
        frame = one_year.copy()
        frame.insert(1, "year", year)
        frame.insert(0, "dmu_id", frame.grid_id + f"_{year}")
        frames.append(frame[list(GRID_COLUMNS)])
    return pd.concat(frames, ignore_index=True)


def export_year(reader: RasterReader, index: pd.DataFrame, year: int) -> pd.DataFrame:
    out = index.loc[index.year.eq(year), list(GRID_COLUMNS)].copy()
    row, col = out.grid_row.to_numpy(int), out.grid_col.to_numpy(int)
    def selected(relative: str, bands: int = 1) -> np.ndarray:
        arr = reader.read(relative, bands)
        return arr[:, row, col] if bands > 1 else arr[row, col]

    out["lulc_dominant_share"] = selected(f"LULC/LULC_dominant_share_{year}_10km.tif")
    out["lulc_valid_area_fraction"] = selected(f"LULC/LULC_valid_area_fraction_{year}_10km.tif")
    fraction_stack = selected(f"LULC/LULC_class_fractions_{year}_10km.tif", bands=9)
    for band, land in FRACTIONS.items():
        out[f"lulc_fraction_{band}_{land}"] = fraction_stack[band - 1]
    for rule, pattern in RULES.items():
        raw = selected(f"LULC/{pattern.format(year=year)}")
        code = np.where(np.isfinite(raw), raw, 0)
        if (not np.isin(code, list(range(10))).all()
                or not np.allclose(code, np.rint(code), rtol=0, atol=1e-8)):
            raise ValueError(f"Invalid LULC code for {rule}, {year}")
        code = code.astype(np.int16)
        out[f"lulc_{rule}_code"] = code
        out[f"lulc_{rule}_land_type"] = [LULC.get(int(v), "Mixed_or_below_threshold") for v in code]
        out[f"lulc_{rule}_assigned"] = code > 0
    for field, pattern in VARIABLES.items():
        out[field] = selected(pattern.format(year=year))
    for field in ["lulc_dominant_share", "lulc_valid_area_fraction", *(
        f"lulc_fraction_{k}_{name}" for k, name in FRACTIONS.items()
    )]:
        values = out[field].to_numpy(float)
        if ((np.isfinite(values) & ((values < -1e-6) | (values > 1.000001))).any()):
            raise ValueError(f"Fraction outside [0, 1] in {field}, {year}")
    return out


def compare(current: pd.DataFrame, previous: pd.DataFrame) -> pd.DataFrame:
    """Field-by-field audit; NaN versus finite is always a mismatch."""
    if (tuple(current.columns) != tuple(previous.columns)
            or current.dmu_id.to_list() != previous.dmu_id.to_list()):
        raise ValueError("New export and existing panel have different columns or DMU rows/order")
    records = []
    for year in YEARS:
        a, b = current.loc[current.year.eq(year)], previous.loc[previous.year.eq(year)]
        for name in current.columns:
            if name not in previous:
                records.append(dict(year=year, field=name, changed_cells=len(a),
                                    new_missing=int(a[name].isna().sum()),
                                    existing_missing=len(a), maximum_absolute_difference=np.nan))
                continue
            if pd.api.types.is_numeric_dtype(a[name]) and not pd.api.types.is_bool_dtype(a[name]):
                x = pd.to_numeric(a[name], errors="coerce").to_numpy(float)
                y = pd.to_numeric(b[name], errors="coerce").to_numpy(float)
                same = np.isclose(x, y, rtol=1e-8, atol=1e-6, equal_nan=True)
                both = np.isfinite(x) & np.isfinite(y)
                maximum = float(np.max(np.abs(x[both] - y[both]))) if both.any() else np.nan
            else:
                x, y = a[name].astype("string"), b[name].astype("string")
                same = x.eq(y).fillna(False) | (x.isna() & y.isna())
                maximum = np.nan
            records.append(dict(year=year, field=name,
                                changed_cells=int(np.count_nonzero(~np.asarray(same))),
                                new_missing=int(a[name].isna().sum()),
                                existing_missing=int(b[name].isna().sum()),
                                maximum_absolute_difference=maximum))
    return pd.DataFrame.from_records(records)


def main() -> int:
    args = arguments()
    weights = read_weights(args.city_weights_source)
    reader = RasterReader(args.raster_root)
    index = make_grid_index(weights, reader)
    frames = [export_year(reader, index, year) for year in YEARS]
    panel = pd.concat(frames, ignore_index=True)
    reader.validate_reference_geometry(panel)
    if set(panel.columns) != set(PANEL_COLUMNS) or len(panel.columns) != len(PANEL_COLUMNS):
        raise ValueError("Internal panel schema does not match the expected 46-column layout")
    panel = panel.loc[:, list(PANEL_COLUMNS)]
    if panel.duplicated(["grid_id", "year"]).any():
        raise ValueError("Repeated grid/year combination in exported panel")
    output = args.output_dir.resolve()
    panel_target = output / "grid_panel_10km.csv"
    if panel_target.is_file():
        previous = pd.read_csv(panel_target, na_values=[NA_TEXT, ""], low_memory=False)
        try:
            audit = compare(panel, previous)
            differences = audit.loc[audit.changed_cells.gt(0)]
            for record in differences.itertuples(index=False):
                logging.warning("%s %s: %s TIFF-derived values differ from existing panel",
                                record.year, record.field, record.changed_cells)
            panel_differs = not differences.empty
        except ValueError as exc:
            logging.warning("Existing panel schema/index differs: %s", exc)
            panel_differs = True
    else:
        panel_differs = True
    if not args.replace_existing and panel_target.exists() and panel_differs:
        raise SystemExit(
            "Existing grid_panel_10km.csv differs from raster-derived output. "
            "Review the reported differences, then rerun with --replace-existing "
            "to back up and replace it. Nothing was overwritten."
        )
    output.mkdir(parents=True, exist_ok=True)
    if panel_differs:
        if panel_target.exists():
            backup(panel_target)
        with tempfile.NamedTemporaryFile(mode="w", dir=output, prefix=".panel_",
                                         suffix=".csv", encoding="utf-8",
                                         newline="", delete=False) as tmp:
            temp = Path(tmp.name)
        try:
            panel.to_csv(temp, index=False, na_rep=NA_TEXT,
                         float_format="%.12g", encoding="utf-8")
            temp.replace(panel_target)
        finally:
            temp.unlink(missing_ok=True)
    logging.info("Grid panel: %s (%s rows, %s columns)",
                 panel_target, len(panel), len(panel.columns))
    if any("2021" in unit for unit in reader.gdp_units):
        logging.warning(
            "GDP TIFF unit tag says 2021 international USD PPP, but the existing "
            "panel header says 2017. Retained the existing header for schema "
            "compatibility; verify source units and rename consistently before publication."
        )
    return 0


def backup(path: Path) -> None:
    destination = path.with_name(path.name + ".before_raster_export.bak")
    if destination.exists():
        raise FileExistsError(f"Backup already exists; move it before replacing {path}: {destination}")
    shutil.copy2(path, destination)
    logging.warning("Backed up existing file to %s", destination)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    raise SystemExit(main())
