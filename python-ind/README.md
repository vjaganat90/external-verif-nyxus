# Python Independent Verification

This folder is intentionally separate from `../matlab` and `../matlab-ind`.

The scripts here use installable third-party Python packages as independent
oracles for Nyxus 2D feature tests. They do not reimplement Nyxus formulas.

The PyRadiomics verifier uses PyRadiomics as the oracle for comparable
first-order and radiomics texture features. PyRadiomics does not expose Nyxus
GLDZM or NGLDM feature classes in version 3.0.1, so those families remain
outside this specific verifier.

The scikit-image verifier is a single package-level script. It combines direct
`skimage.measure.regionprops` measurements for remaining 2D shape rows and
direct `regionprops` moment properties for the added 2D tests. Rows that would
require normalized raw-moment composition, Nyxus contour-distance weighting, or
morphology logic built from lower-level primitives are not mixed into the
comparison output; they are listed in `skimage_incomparable_no_oracle.csv`.

The MIRP verifier uses MIRP as a third-party oracle for GLDZM and NGLDM rows
that MIRP exposes as direct named features. MIRP 2.6.0 does not expose the
Nyxus mean-style rows `GLDZM_GLM`, `GLDZM_ZDM`, `NGLDM_GLM`, or `NGLDM_DCM` as
direct feature columns, so those rows are listed in
`mirp_incomparable_no_oracle.csv` rather than counted as oracle comparisons.

The Centrosome verifier uses CellProfiler's image-processing support library
as an independent oracle for the remaining shape and neighbor candidates where
Centrosome exposes direct measurements. `minimum_enclosing_circle` agrees with
Nyxus; `find_neighbors` is logged as a convention mismatch because Centrosome
uses 8-connected adjacency and therefore counts diagonal contacts in the
neighbor fixture.

The cp-measure verifier uses a package-level CellProfiler measurement port for
direct edge-intensity and neighbor measurements. It supersedes the earlier
Centrosome `NUM_NEIGHBORS` mismatch in the tracker because cp-measure's
`MeasureObjectNeighbors` adjacent mode matches the Nyxus neighbor-count
fixture. It also documents direct convention mismatches for mass displacement,
compactness, percent touching, and second closest-neighbor distance.

The DIPlib verifier uses DIPlib's `MeasurementTool` direct shape measurements.
It exposes named `Roundness` and `Circularity` measurements, but both are
convention mismatches for the current Nyxus shape fixture.

The imea verifier uses `imea.shape_measurements_2d` direct dataframe columns
for the remaining 2D shape long-tail rows where imea exposes feature-level
measurements. It does not use imea's optional raw chord/statistical-length
arrays, because reducing those arrays into Nyxus scalar values would be formula
recomputation. On the current fixture, 37 rows are comparable: `ALLCHORDS_MIN`
agrees within 5%, and the other 36 are convention mismatches.

The `remaining_2d_unique_incomparable_triage.csv` file is the strict audit for
features still marked `Not available`: it lists rows where the packages tried
so far expose only related primitives, a different convention, or no direct
built-in scalar at all.

## Environment

PyRadiomics was installed into the existing conda environment:

```bash
conda run -n nyxus_ci python -m pip install pyradiomics==3.0.1
```

scikit-image was installed into the same environment:

```bash
conda install -n nyxus_ci -y -c conda-forge scikit-image
```

MIRP was installed into the same environment:

```bash
conda run -n nyxus_ci python -m pip install mirp==2.6.0
```

Centrosome was installed into the same environment:

```bash
conda run -n nyxus_ci python -m pip install centrosome==1.3.3
```

cp-measure and DIPlib were installed into the same environment:

```bash
conda run -n nyxus_ci python -m pip install cp-measure==0.1.19 diplib==3.6.0
```

imea was installed into the same environment:

```bash
conda run -n nyxus_ci python -m pip install imea==0.3.5
```

DIPlib's pip wheel loads a separate OpenMP runtime in this macOS arm64 conda
environment. The DIPlib verifier sets `KMP_DUPLICATE_LIB_OK=TRUE` explicitly
and logs that caveat. Attempts to install `diplib`/`diplib-python` from
conda-forge in `nyxus_ci` did not find an osx-arm64 package.

On Python 3.12, the published PyRadiomics 3.0.1 source distribution needs an
install-time compatibility shim for legacy `configparser` APIs. That shim is
only needed while building the package; the verifier itself imports
PyRadiomics normally from `nyxus_ci`.

The mahotas verifier compares the Nyxus ZERNIKE2D regression vector against
`mahotas.features.zernike_moments`. Both the expected coefficients and the ROI
fixture are parsed from the Nyxus repo (no values are hard-coded in the script).
The Zernike normalization/order convention differs from Nyxus, so it is counted
as comparable built-in coverage rather than agreement.

## Nyxus repository source

None of these scripts hard-code a path to a local Nyxus checkout. They call the
shared helper `nyxus_repo.py`, which looks for an existing checkout in this
order:

- Path passed as the script's first argument.
- `NYXUS_REPO_DIR`.
- Common local locations: `../nyxus`, the current directory, and
  `../.cache/nyxus`.

The helper does not fetch, pull, checkout, or reset an existing checkout. For
fresh machines only, set `NYXUS_ALLOW_CLONE=1` to allow a one-time clone into
`../.cache/nyxus`. Optional clone knobs are `NYXUS_REPO_URL` and
`NYXUS_BRANCH`; if `NYXUS_BRANCH` is unset, git uses the remote default branch.

## Run

```bash
# Use an existing local checkout:
conda run -n nyxus_ci python python-ind/run_all_pyradiomics_independent_2d_checks.py
conda run -n nyxus_ci python python-ind/run_skimage_independent_2d_checks.py
conda run -n nyxus_ci python python-ind/run_mirp_independent_2d_texture_checks.py
conda run -n nyxus_ci python python-ind/run_centrosome_independent_2d_shape_neighbor_checks.py
conda run -n nyxus_ci python python-ind/run_cp_measure_independent_2d_shape_neighbor_checks.py
conda run -n nyxus_ci python python-ind/run_diplib_independent_2d_shape_checks.py
conda run -n nyxus_ci python python-ind/run_imea_independent_2d_shape_checks.py
conda run -n nyxus_ci python python-ind/run_mahotas_remaining2d_zernike_check.py

# Use an existing local checkout (no network):
conda run -n nyxus_ci python python-ind/run_skimage_independent_2d_checks.py /path/to/local/nyxus

# Allow a one-time clone into .cache/nyxus, optionally selecting a branch:
NYXUS_ALLOW_CLONE=1 NYXUS_BRANCH=main conda run -n nyxus_ci python python-ind/run_skimage_independent_2d_checks.py
```

## Outputs

The runner writes CSV reports next to these scripts:

- `pyradiomics_independent_2d_results.csv`
- `pyradiomics_independent_2d_summary.csv`
- `skimage_independent_2d_results.csv`
- `skimage_independent_2d_summary.csv`
- `skimage_incomparable_no_oracle.csv`
- `mirp_independent_2d_texture_results.csv`
- `mirp_independent_2d_texture_summary.csv`
- `mirp_incomparable_no_oracle.csv`
- `centrosome_independent_2d_shape_neighbor_results.csv`
- `centrosome_independent_2d_shape_neighbor_summary.csv`
- `cp_measure_independent_2d_shape_neighbor_results.csv`
- `cp_measure_independent_2d_shape_neighbor_summary.csv`
- `diplib_independent_2d_shape_results.csv`
- `diplib_independent_2d_shape_summary.csv`
- `imea_independent_2d_shape_results.csv`
- `imea_independent_2d_shape_summary.csv`
- `mahotas_remaining2d_zernike_results.csv`
- `mahotas_remaining2d_zernike_summary.csv`
- `remaining_2d_unique_incomparable_triage.csv`

## Interpretation

Rows marked `PASS` mean the current Nyxus test constant agrees within 5% with
the named third-party oracle under the settings reported in the CSV.

Rows marked `MISMATCH` mean the third-party package exposes a comparable
feature, but the current Nyxus asserted value differs by more than 5%. These
usually point to different discretization, ROI masking, direction aggregation,
or convention choices, not automatically to a Nyxus bug.

Every printed row and CSV row includes absolute difference, relative
difference, and tolerance fields so convention drift remains visible even when
the row passes.
