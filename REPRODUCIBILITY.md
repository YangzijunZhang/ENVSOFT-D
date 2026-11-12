# Reproducibility guide

Run each command from the repository root in the order shown in [README.md](README.md). All paths are relative; no script reads files outside the repository.

## Workflow

1. `raster_to_panel.py` rebuilds `grid_panel_10km.csv`, `city_grid_weights_10km.csv` and `city_reference.csv` from the rasters and boundaries in `data/raster data/`. It compares them with the files in `data/` and writes the comparison to `results/raster_to_panel/`.
2. `city_aggregation.py` creates `data/city_panel_10km.csv` from the grid panel, the city-grid weights and the city reference table.
3. `dea_performance.py` estimates the grid and city scores and writes one CSV per resource system and scale to `results/performance/`.
4. The remaining scripts read these CSVs and have narrow responsibilities:
   - `appendix_c.py`: Table C1 and Figs. C1-C3.
   - `analysis_4_1.py`: grid and city maps for all systems and benchmark years (Figs. 3-5).
   - `analysis_4_2.py`: regional cropland summary and Fig. 6.
   - `analysis_4_3.py`: grid and city SDE parameters, displacements and Fig. 7.
   - `specification_checks.py`: Appendix B.3, the ET coverage statistics in Appendix C.1 and Appendix E.
   - `appendix_d.py`: re-estimation under the 50% and 70% rules (Tables D1-D2, Figs. D1-D6) and the subsampled frontiers (Fig. D7).

## Checks built into the scripts

- **Rebuilt inputs.** `results/raster_to_panel/comparison_with_reference.csv` compares every column of the rebuilt CSVs with the files in `data/`. `qc_summary.csv` reports grid and coverage checks, the consistency of the LULC assignment rules and the sum of the city weights in each cell. `input_manifest.csv` lists every input file with its SHA-256 checksum.
- **Solver.** Every score row records whether the linear programme was solved and its largest relative constraint residual.
- **Appendix D.** `results/appendix_d/run_config.json` records the random seed, the number of replicates, the sampling fraction, the Python and package versions and the SHA-256 checksums of the inputs. The subsampling uses a separate seed for each system and replicate, so the results do not depend on the number of workers.

## Expected numerical differences

The rebuilt CSVs match the files in `data/` to the precision stored, with two exceptions:

- `latitude_wgs84` in the grid panel can differ in the twelfth significant digit (about 10⁻¹⁰ degrees) for a few cells, because coordinate transformations round differently between PROJ builds.
- The polygon areas and allocation weights in the city tables can differ by up to 10⁻⁷ in relative terms, because exact intersections differ slightly between GEOS builds.

Both differences are far below the precision of the data and do not change any reported result. The comparison file labels such columns "equal within" the stated tolerance; any other difference is labelled `DIFFERENT`.

## Formats and environment

Generated CSVs encode missing values as `NA` and numbers with 12 significant digits. Package versions are pinned in `requirements.txt` and `environment.yml`. The results were produced with Python 3.14.3 on macOS (Apple silicon).

## Runtime

On one core, `raster_to_panel.py` takes about 10 seconds, `specification_checks.py` about 8 minutes and the threshold stage of `appendix_d.py` a few minutes. The subsampling stage of `appendix_d.py` solves the frontiers 100 times and takes about 2.5 hours; `--workers` runs the replicates in parallel.

## Map layers

The figures use the projected WKT geometries in `data/cartography_layers.csv`: the basin boundary, the upper, middle and lower reaches, the Yellow River and the 76 city boundaries clipped to the basin. These geometries are simplified and are used only for drawing; all spatial calculations use the shapefiles in `data/raster data/adBoundary/YRB/`.
