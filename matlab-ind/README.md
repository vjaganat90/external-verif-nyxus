# MATLAB Independent Verification

This folder is intentionally separate from `../matlab`.

The older `matlab` scripts are MATLAB reimplementations of Nyxus or IBSI-style
feature formulas. They are useful drift guards, but they are not always
third-party oracle checks.

The single verifier script in this folder uses MATLAB built-in functions as the
oracle wherever MATLAB exposes a comparable feature or core algorithm. The
comparison CSVs contain only actual oracle comparisons. Rows with no direct
MATLAB oracle are listed separately in `matlab_incomparable_no_oracle.csv`.

## Nyxus repository source

The verifier does not hard-code a path to a local Nyxus checkout. The entry
point calls `ensure_nyxus_repo`, which looks for an existing checkout in this
order:

- Path passed to `run_all_matlab_independent_2d_checks('/path/to/nyxus')`.
- `NYXUS_REPO_DIR`.
- Common local locations: `../nyxus`, the current directory, and
  `../.cache/nyxus`.

The helper does not fetch, pull, checkout, or reset an existing checkout. For
fresh machines only, set `NYXUS_ALLOW_CLONE=1` to allow a one-time clone into
`../.cache/nyxus`. Optional clone knobs are `NYXUS_REPO_URL` and
`NYXUS_BRANCH`; if `NYXUS_BRANCH` is unset, git uses the remote default branch.

## Run

From MATLAB (run from the external-verif repo root, or `cd` there first):

```matlab
addpath('matlab-ind')
run_all_matlab_independent_2d_checks            % find an existing local checkout
run_all_matlab_independent_2d_checks('/path/to/local/nyxus')  % use a local checkout
```

`run_all_matlab_independent_2d_checks.m` is a single-file program; its first
function is the entry point and a FUNCTION MAP comment at the top lists every
helper grouped by `%%` section.

Two first-order rows use genuine independent MATLAB built-ins rather than a
re-derivation of the Nyxus formula: `MEDIAN_ABSOLUTE_DEVIATION` is checked
against `mad(x,1)` (textbook median-of-deviations; Nyxus uses the mean of the
deviations, so the small residual difference is a real convention gap, not a
self-consistency artifact), and `QCOD` is composed transparently from MATLAB
`prctile` percentiles since MATLAB exposes no single QCOD built-in.

## Outputs

The runner writes CSV reports next to these scripts:

- `first_order_builtin_independent_results.csv`
- `shape2d_builtin_independent_results.csv`
- `glcm_builtin_independent_results.csv`
- `remaining2d_feret_builtin_independent_results.csv`
- `matlab_incomparable_no_oracle.csv`
- `matlab_independent_2d_summary.csv`

All MATLAB comparisons are launched from the one mega script
`run_all_matlab_independent_2d_checks.m`; the output CSVs are split only to
keep reviewer-facing feature families readable.

## Interpretation

Rows marked `PASS` mean the current Nyxus test constant agrees within 5% with
a MATLAB built-in or a direct composition of MATLAB built-ins under the reported
convention. The CSV outputs and MATLAB command-window log include `AbsDiff`,
`RelDiffPct`, and `TolerancePct` so convention drift remains visible even when
the row passes.

Rows marked `MISMATCH` are important. They usually mean MATLAB and Nyxus use
different documented conventions, or the Nyxus test constant is stale.

Rows with no direct MATLAB oracle are not mixed into the comparison outputs;
look at `matlab_incomparable_no_oracle.csv` for that list.
