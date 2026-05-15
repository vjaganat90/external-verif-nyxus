# MATLAB Independent Verification

This folder is intentionally separate from `../matlab`.

The older `matlab` scripts are MATLAB reimplementations of Nyxus or IBSI-style
feature formulas. They are useful drift guards, but they are not always
third-party oracle checks.

Scripts in this folder use MATLAB built-in functions as the oracle wherever
MATLAB exposes a comparable feature or core algorithm. When MATLAB has no
matching built-in, the scripts report that explicitly instead of pretending
that a formula reimplementation is third-party verification.

## Run

From MATLAB:

```matlab
addpath('/Users/jaganathv2/Code/Codex/external-verif-nyxus/matlab-ind')
run_all_matlab_independent_2d_checks('/Users/jaganathv2/Code/Codex/nyxus')
```

If `repoRoot` is omitted, the scripts use `/Users/jaganathv2/Code/Codex/nyxus`.

## Outputs

The runner writes CSV reports next to these scripts:

- `first_order_builtin_independent_results.csv`
- `shape2d_builtin_independent_results.csv`
- `glcm_builtin_independent_results.csv`
- `unsupported_builtin_families.csv`
- `matlab_independent_2d_summary.csv`

## Interpretation

Rows marked `PASS` mean the current Nyxus test constant agrees within 5% with
a MATLAB built-in or a direct composition of MATLAB built-ins under the reported
convention. The CSV outputs and MATLAB command-window log include `AbsDiff`,
`RelDiffPct`, and `TolerancePct` so convention drift remains visible even when
the row passes.

Rows marked `MISMATCH` are important. They usually mean MATLAB and Nyxus use
different documented conventions, or the Nyxus test constant is stale.

Rows marked `NO_BUILTIN` mean MATLAB does not provide a comparable high-level
feature oracle for that feature family.
