#!/usr/bin/env python3
"""Run DIPlib independent checks for Nyxus 2D shape gaps."""

from __future__ import annotations

import csv
import importlib.metadata
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path

# The published macOS arm64 DIPlib wheel loads a separate OpenMP runtime in this
# conda environment. Conda-forge does not currently provide a compatible DIPlib
# package here, so keep the workaround explicit and logged in the output.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import diplib as dip  # noqa: E402
import diplib.PyDIP_bin.MeasurementTool as dip_measure  # noqa: E402
import numpy as np  # noqa: E402


from nyxus_repo import ensure_nyxus_repo  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent
TOLERANCE_PCT = 5.0


@dataclass(frozen=True)
class Row:
    family: str
    feature: str
    expected: float
    actual: float
    abs_diff: float
    rel_diff_pct: float
    tolerance_pct: float
    status: str
    oracle_feature: str
    image_source: str
    settings: str
    note: str


def relative_diff_pct(actual: float, expected: float) -> float:
    if math.isnan(actual) or math.isnan(expected):
        return math.inf
    if expected == 0:
        return 0.0 if abs(actual) <= 1e-12 else math.inf
    return 100.0 * abs(actual - expected) / abs(expected)


def result_row(feature: str, expected: float, actual: float, oracle_feature: str, note: str) -> Row:
    abs_diff = abs(actual - expected) if not math.isnan(actual) else math.nan
    rel_diff_pct = relative_diff_pct(actual, expected)
    status = "PASS" if rel_diff_pct <= TOLERANCE_PCT else "MISMATCH"
    return Row(
        family="Shape and morphology subset",
        feature=feature,
        expected=expected,
        actual=actual,
        abs_diff=abs_diff,
        rel_diff_pct=rel_diff_pct,
        tolerance_pct=TOLERANCE_PCT,
        status=status,
        oracle_feature=oracle_feature,
        image_source="tests/test_data.h shape2d_morphology_mask",
        settings="DIPlib MeasurementTool on 2D label image; KMP_DUPLICATE_LIB_OK=TRUE required for the pip wheel in nyxus_ci",
        note=note,
    )


def parse_pixel_array(test_data_text: str, array_name: str) -> list[tuple[int, int, float]]:
    pattern = (
        r"const\s+static\s+NyxusPixel\s+"
        + re.escape(array_name)
        + r"\[\]\s*=\s*\{([\s\S]*?)\};"
    )
    match = re.search(pattern, test_data_text)
    if not match:
        raise ValueError(f"Could not find {array_name}")
    triples = []
    for x, y, intensity in re.findall(
        r"\{\s*(?:/\*x\*/)?\s*([0-9.+\-eE]+)\s*,\s*(?:/\*y\*/)?\s*([0-9.+\-eE]+)\s*,\s*([0-9.+\-eE]+)\s*\}",
        match.group(1),
    ):
        triples.append((int(float(x)), int(float(y)), float(intensity)))
    return triples


def parse_shape_truth(header_text: str) -> dict[str, float]:
    match = re.search(
        r"static\s+std::unordered_map<std::string,\s*double>\s+shape2d_truth\s*\{([\s\S]*?)\};",
        header_text,
    )
    if not match:
        raise ValueError("Could not find shape2d_truth")
    return {
        feature: float(value)
        for feature, value in re.findall(
            r'\{\s*"([^"]+)"\s*,\s*([0-9.+\-eE]+)\s*\}', match.group(1)
        )
    }


def to_image(triples: list[tuple[int, int, float]], dtype: type[np.integer] | type[np.floating]) -> np.ndarray:
    width = max(x for x, _, _ in triples) + 1
    height = max(y for _, y, _ in triples) + 1
    image = np.zeros((height, width), dtype=dtype)
    for x, y, value in triples:
        image[y, x] = value
    return image


def diplib_version() -> str:
    return importlib.metadata.version("diplib")


def main(argv: list[str]) -> int:
    repo_root = ensure_nyxus_repo(argv[1] if len(argv) > 1 else None)
    test_data_text = (repo_root / "tests" / "test_data.h").read_text(encoding="utf-8")
    shape_truth = parse_shape_truth((repo_root / "tests" / "test_shape_morphology_2d.h").read_text(encoding="utf-8"))

    labels = (to_image(parse_pixel_array(test_data_text, "shape2d_morphology_mask"), np.uint32) > 0).astype(np.uint32)
    measurements = dip_measure.Measure(
        dip.Image(labels),
        features=["Roundness", "Circularity"],
        objectIDs=[1],
    )

    rows = [
        result_row(
            "ROUNDNESS",
            shape_truth["ROUNDNESS"],
            float(measurements[1]["Roundness"][0]),
            "diplib.MeasurementTool Roundness",
            "Direct DIPlib roundness measurement; convention differs from Nyxus' ellipse-fitting roundness.",
        ),
        result_row(
            "CIRCULARITY",
            shape_truth["CIRCULARITY"],
            float(measurements[1]["Circularity"][0]),
            "diplib.MeasurementTool Circularity",
            "Direct DIPlib circularity measurement; convention differs from Nyxus' sqrt(4*pi*area/perimeter^2) value using Nyxus contour perimeter.",
        ),
    ]

    result_path = OUT_DIR / "diplib_independent_2d_shape_results.csv"
    summary_path = OUT_DIR / "diplib_independent_2d_shape_summary.csv"

    with result_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(Row.__dataclass_fields__.keys())
        for row in rows:
            writer.writerow([getattr(row, field) for field in Row.__dataclass_fields__])

    pass_count = sum(row.status == "PASS" for row in rows)
    mismatch_count = len(rows) - pass_count
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["Library", "Version", "Rows", "PASS", "MISMATCH", "Note"])
        writer.writerow(
            [
                "DIPlib",
                diplib_version(),
                len(rows),
                pass_count,
                mismatch_count,
                "Conda-forge package unavailable for this osx-arm64 env; pip wheel requires KMP_DUPLICATE_LIB_OK=TRUE here.",
            ]
        )

    for row in rows:
        print(
            f"{row.family:30s} {row.feature:42s} "
            f"expected={row.expected:.15g} diplib={row.actual:.15g} "
            f"absdiff={row.abs_diff:.15g} reldiff={row.rel_diff_pct:.6g}% "
            f"tol={row.tolerance_pct:.2f}% status={row.status}"
        )
    print(
        f"DIPlib independent 2D shape rows: {pass_count} pass, {mismatch_count} mismatch. "
        f"Wrote {result_path} and {summary_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
