# GitHub and Zenodo release

Each GitHub release is archived by Zenodo, which assigns it a permanent DOI. A published DOI cannot be withdrawn, so complete steps 1-3 before creating the release.

## Before the release

1. Run the full workflow in [README.md](README.md) from a clean environment and check that:
   - `results/raster_to_panel/comparison_with_reference.csv` labels no column `DIFFERENT`;
   - `results/`, `figures/` and the values in the manuscript and SI agree, including `results/appendix_d/subsample_summary.csv` and `figures/Fig_D7.*`.
2. Check the repository contents:
   - it contains the files listed in [MANIFEST.md](MANIFEST.md) and the source rasters and boundaries in `data/raster data/`;
   - it does not contain `.venv/`, `__pycache__/`, `.DS_Store` or `*.aux.xml` files;
   - no single file exceeds 50 MB (GitHub warns above 50 MB and rejects files above 100 MB).
3. Check the metadata:
   - `.zenodo.json`: title, authors, licence and description. Zenodo uses `.zenodo.json` and ignores `CITATION.cff` when both are present.
   - `CITATION.cff`: version `1.0.0` and `date-released` set to the release date.

## Release

4. Make the GitHub repository `YangzijunZhang/ENVSOFT-D` public.
5. Sign in to <https://zenodo.org> with the GitHub account. Under the user menu, open **GitHub**, click **Sync now** and switch `ENVSOFT-D` on. Releases created before this step are not archived.
6. On GitHub, open **Releases > Draft a new release**. Create the tag `v1.0.0` on `main`, give the release a title (for example "v1.0.0 - revised manuscript ENVSOFT-D-26-02569") and publish it.
7. After a few minutes, open the new record on Zenodo and check the files and metadata. Zenodo assigns a version DOI (this release) and a concept DOI (all versions). Cite the version DOI in the manuscript.

## After the release

8. Add the version DOI to the manuscript's Software and data availability statement, to `README.md` (for example as a DOI badge), to the GitHub release notes and to the repository description.
9. These edits come after the release and are not part of the v1.0.0 archive. A new release (for example `v1.0.1`) is needed only if files in the archive itself must change; it receives its own version DOI.

Do not insert a provisional or invented DOI. Zenodo supplies the DOI when the release is archived.
