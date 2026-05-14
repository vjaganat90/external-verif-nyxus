# Python Independent Verification

This folder is intentionally separate from `../matlab` and `../matlab-ind`.

The scripts here use installable third-party Python packages as independent
oracles for Nyxus 2D feature tests. They do not reimplement Nyxus formulas.

The first verifier uses PyRadiomics as the oracle for comparable first-order
and radiomics texture features. PyRadiomics does not expose Nyxus GLDZM or
NGLDM feature classes in version 3.0.1, so those families remain outside this
specific verifier.

## Environment

PyRadiomics was installed into the existing conda environment:

```bash
conda run -n nyxus_ci python -m pip install pyradiomics==3.0.1
```

On Python 3.12, the published PyRadiomics 3.0.1 source distribution needs an
install-time compatibility shim for legacy `configparser` APIs. That shim is
only needed while building the package; the verifier itself imports
PyRadiomics normally from `nyxus_ci`.

## Run

```bash
conda run -n nyxus_ci python /Users/jaganathv2/Code/Codex/external-verif-nyxus/python-ind/run_all_pyradiomics_independent_2d_checks.py /Users/jaganathv2/Code/Codex/nyxus
```

If `repoRoot` is omitted, the script uses
`/Users/jaganathv2/Code/Codex/nyxus`.

## Outputs

The runner writes CSV reports next to these scripts:

- `pyradiomics_independent_2d_results.csv`
- `pyradiomics_independent_2d_summary.csv`

## Interpretation

Rows marked `PASS` mean the current Nyxus test constant agrees within 5% with
the PyRadiomics value under the settings reported in the CSV.

Rows marked `MISMATCH` mean PyRadiomics exposes a comparable feature, but the
current Nyxus asserted value differs by more than 5%. These usually point to
different discretization, ROI masking, direction aggregation, or convention
choices, not automatically to a Nyxus bug.

Every printed row and CSV row includes `AbsDiff`, `RelDiffPct`, and
`TolerancePct` so convention drift remains visible even when the row passes.
