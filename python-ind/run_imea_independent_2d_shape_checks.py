#!/usr/bin/env python3
"""Run imea independent checks for Nyxus 2D shape long-tail rows.

This verifier intentionally uses only direct columns returned by
imea.shape_measurements_2d. Optional raw statistical-length/chord arrays are
not used, because reducing those arrays into Nyxus scalar features would be
formula recomputation rather than a third-party built-in oracle.
"""

from __future__ import annotations

import csv
import importlib.metadata
import math
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from imea import shape_measurements_2d


from nyxus_repo import ensure_nyxus_repo

OUT_DIR = Path(__file__).resolve().parent
RESULT_CSV = OUT_DIR / "imea_independent_2d_shape_results.csv"
SUMMARY_CSV = OUT_DIR / "imea_independent_2d_shape_summary.csv"
TOLERANCE_PCT = 5.0
VERIFIABLE_BUCKET = "VERIFIABLE_WITH_3P_BUILTIN_ORACLE"


@dataclass(frozen=True)
class CheckSpec:
    family: str
    feature: str
    truth_map: str
    imea_column: str
    cxx_test_name: str
    note: str


@dataclass(frozen=True)
class Row:
    family: str
    feature: str
    review_name: str
    verification_bucket: str
    cxx_test_name: str
    expected: float
    actual: float
    abs_diff: float
    rel_diff_pct: float
    tolerance_pct: float
    status: str
    oracle_tool: str
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


def result_row(spec: CheckSpec, expected: float, actual: float) -> Row:
    abs_diff = abs(actual - expected)
    rel_diff_pct = relative_diff_pct(actual, expected)
    status = "PASS" if rel_diff_pct <= TOLERANCE_PCT else "MISMATCH"
    return Row(
        family=spec.family,
        feature=spec.feature,
        review_name=f"{VERIFIABLE_BUCKET}__{spec.feature}",
        verification_bucket=VERIFIABLE_BUCKET,
        cxx_test_name=spec.cxx_test_name,
        expected=expected,
        actual=actual,
        abs_diff=abs_diff,
        rel_diff_pct=rel_diff_pct,
        tolerance_pct=TOLERANCE_PCT,
        status=status,
        oracle_tool="imea",
        oracle_feature=f"shape_measurements_2d.{spec.imea_column}",
        image_source="tests/test_data.h shape2d_morphology_mask",
        settings="imea.shape_measurements_2d(binary_mask, spatial_resolution_xy=1, dalpha=9)",
        note=spec.note,
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


def to_binary_image(triples: list[tuple[int, int, float]]) -> np.ndarray:
    width = max(x for x, _, _ in triples) + 1
    height = max(y for _, y, _ in triples) + 1
    image = np.zeros((height, width), dtype=bool)
    for x, y, value in triples:
        image[y, x] = value > 0
    return image


def parse_truth_map(header_text: str, map_name: str) -> dict[str, float]:
    match = re.search(
        r"static\s+std::unordered_map<std::string,\s*double>\s+"
        + re.escape(map_name)
        + r"\s*\{([\s\S]*?)\};",
        header_text,
    )
    if not match:
        raise ValueError(f"Could not find {map_name}")
    return {
        feature: float(value)
        for feature, value in re.findall(
            r'\{\s*"([^"]+)"\s*,\s*([0-9.+\-eE]+)\s*\}', match.group(1)
        )
    }


def imea_version() -> str:
    return importlib.metadata.version("imea")


def check_specs() -> list[CheckSpec]:
    shape_fractal_circle_test = "TEST_SHAPE2D_VERIFIABLE_WITH_3P_BUILTIN_ORACLE_FRACTAL_CIRCLE_FEATURES"
    shape_geodetic_test = "TEST_SHAPE2D_VERIFIABLE_WITH_3P_BUILTIN_ORACLE_GEODETIC_THICKNESS_EROSION"
    caliper_test = "TEST_REMAINING2D_VERIFIABLE_WITH_3P_BUILTIN_ORACLE_CALIPER_FEATURES"
    chord_test = "TEST_REMAINING2D_VERIFIABLE_WITH_3P_BUILTIN_ORACLE_CHORD_STAT_FEATURES"
    return [
        CheckSpec(
            "Shape and morphology subset",
            "DIAMETER_EQUAL_PERIMETER",
            "shape2d_truth",
            "diameter_equal_perimeter",
            "TEST_SHAPE2D_VERIFIABLE_WITH_3P_BUILTIN_ORACLE_CONTOUR_DIAMETER_EQUAL_PERIMETER",
            "Direct imea equal-perimeter diameter column.",
        ),
        CheckSpec(
            "Shape and morphology subset",
            "FRACT_DIM_BOXCOUNT",
            "shape2d_truth",
            "fractal_dimension_boxcounting_method",
            shape_fractal_circle_test,
            "Direct imea box-counting fractal-dimension column.",
        ),
        CheckSpec(
            "Shape and morphology subset",
            "FRACT_DIM_PERIMETER",
            "shape2d_truth",
            "fractal_dimension_perimeter_method",
            shape_fractal_circle_test,
            "Direct imea perimeter-method fractal-dimension column.",
        ),
        CheckSpec(
            "Shape and morphology subset",
            "DIAMETER_CIRCUMSCRIBING_CIRCLE",
            "shape2d_truth",
            "diameter_circumscribing_circle",
            shape_fractal_circle_test,
            "Direct imea circumscribing-circle diameter column.",
        ),
        CheckSpec(
            "Shape and morphology subset",
            "DIAMETER_INSCRIBING_CIRCLE",
            "shape2d_truth",
            "diameter_inscribing_circle",
            shape_fractal_circle_test,
            "Direct imea inscribing-circle diameter column.",
        ),
        CheckSpec(
            "Added 2D shape morphology extras",
            "GEODETIC_LENGTH",
            "shape2d_truth",
            "geodeticlength",
            shape_geodetic_test,
            "Direct imea geodetic-length column.",
        ),
        CheckSpec(
            "Added 2D shape morphology extras",
            "THICKNESS",
            "shape2d_truth",
            "thickness",
            shape_geodetic_test,
            "Direct imea thickness column.",
        ),
        CheckSpec(
            "Added 2D shape morphology extras",
            "EROSIONS_2_VANISH",
            "shape2d_truth",
            "n_erosions",
            shape_geodetic_test,
            "Direct imea erosion-count column.",
        ),
        CheckSpec(
            "Remaining 2D erosion complement",
            "EROSIONS_2_VANISH_COMPLEMENT",
            "remaining2d_truth",
            "n_erosions_complement",
            "TEST_REMAINING2D_VERIFIABLE_WITH_3P_BUILTIN_ORACLE_EROSION_COMPLEMENT",
            "Direct imea complement erosion-count column.",
        ),
        CheckSpec(
            "Remaining 2D Feret calipers",
            "STAT_FERET_DIAM_MEAN",
            "remaining2d_truth",
            "feret_mean",
            caliper_test,
            "Direct imea Feret mean column.",
        ),
        CheckSpec(
            "Remaining 2D Feret calipers",
            "STAT_FERET_DIAM_MEDIAN",
            "remaining2d_truth",
            "feret_median",
            caliper_test,
            "Direct imea Feret median column.",
        ),
        CheckSpec(
            "Remaining 2D Feret calipers",
            "STAT_FERET_DIAM_STDDEV",
            "remaining2d_truth",
            "feret_std",
            caliper_test,
            "Direct imea Feret standard-deviation column.",
        ),
        CheckSpec(
            "Remaining 2D Feret calipers",
            "STAT_FERET_DIAM_MODE",
            "remaining2d_truth",
            "feret_mode",
            caliper_test,
            "Direct imea Feret mode column.",
        ),
        CheckSpec("Remaining 2D Martin calipers", "STAT_MARTIN_DIAM_MIN", "remaining2d_truth", "martin_min", caliper_test, "Direct imea Martin minimum column."),
        CheckSpec("Remaining 2D Martin calipers", "STAT_MARTIN_DIAM_MAX", "remaining2d_truth", "martin_max", caliper_test, "Direct imea Martin maximum column."),
        CheckSpec("Remaining 2D Martin calipers", "STAT_MARTIN_DIAM_MEAN", "remaining2d_truth", "martin_mean", caliper_test, "Direct imea Martin mean column."),
        CheckSpec("Remaining 2D Martin calipers", "STAT_MARTIN_DIAM_MEDIAN", "remaining2d_truth", "martin_median", caliper_test, "Direct imea Martin median column."),
        CheckSpec("Remaining 2D Martin calipers", "STAT_MARTIN_DIAM_STDDEV", "remaining2d_truth", "martin_std", caliper_test, "Direct imea Martin standard-deviation column."),
        CheckSpec("Remaining 2D Martin calipers", "STAT_MARTIN_DIAM_MODE", "remaining2d_truth", "martin_mode", caliper_test, "Direct imea Martin mode column."),
        CheckSpec("Remaining 2D Nassenstein calipers", "STAT_NASSENSTEIN_DIAM_MIN", "remaining2d_truth", "nassenstein_min", caliper_test, "Direct imea Nassenstein minimum column."),
        CheckSpec("Remaining 2D Nassenstein calipers", "STAT_NASSENSTEIN_DIAM_MAX", "remaining2d_truth", "nassenstein_max", caliper_test, "Direct imea Nassenstein maximum column."),
        CheckSpec("Remaining 2D Nassenstein calipers", "STAT_NASSENSTEIN_DIAM_MEAN", "remaining2d_truth", "nassenstein_mean", caliper_test, "Direct imea Nassenstein mean column."),
        CheckSpec("Remaining 2D Nassenstein calipers", "STAT_NASSENSTEIN_DIAM_MEDIAN", "remaining2d_truth", "nassenstein_median", caliper_test, "Direct imea Nassenstein median column."),
        CheckSpec("Remaining 2D Nassenstein calipers", "STAT_NASSENSTEIN_DIAM_STDDEV", "remaining2d_truth", "nassenstein_std", caliper_test, "Direct imea Nassenstein standard-deviation column."),
        CheckSpec("Remaining 2D Nassenstein calipers", "STAT_NASSENSTEIN_DIAM_MODE", "remaining2d_truth", "nassenstein_mode", caliper_test, "Direct imea Nassenstein mode column."),
        CheckSpec("Remaining 2D chords", "MAXCHORDS_MAX", "remaining2d_truth", "maxchords_max", chord_test, "Direct imea maximum-chords maximum column."),
        CheckSpec("Remaining 2D chords", "MAXCHORDS_MIN", "remaining2d_truth", "maxchords_min", chord_test, "Direct imea maximum-chords minimum column."),
        CheckSpec("Remaining 2D chords", "MAXCHORDS_MEDIAN", "remaining2d_truth", "maxchords_median", chord_test, "Direct imea maximum-chords median column."),
        CheckSpec("Remaining 2D chords", "MAXCHORDS_MEAN", "remaining2d_truth", "maxchords_mean", chord_test, "Direct imea maximum-chords mean column."),
        CheckSpec("Remaining 2D chords", "MAXCHORDS_MODE", "remaining2d_truth", "maxchords_mode", chord_test, "Direct imea maximum-chords mode column."),
        CheckSpec("Remaining 2D chords", "MAXCHORDS_STDDEV", "remaining2d_truth", "maxchords_std", chord_test, "Direct imea maximum-chords standard-deviation column."),
        CheckSpec("Remaining 2D chords", "ALLCHORDS_MAX", "remaining2d_truth", "allchords_max", chord_test, "Direct imea all-chords maximum column."),
        CheckSpec("Remaining 2D chords", "ALLCHORDS_MIN", "remaining2d_truth", "allchords_min", chord_test, "Direct imea all-chords minimum column."),
        CheckSpec("Remaining 2D chords", "ALLCHORDS_MEDIAN", "remaining2d_truth", "allchords_median", chord_test, "Direct imea all-chords median column."),
        CheckSpec("Remaining 2D chords", "ALLCHORDS_MEAN", "remaining2d_truth", "allchords_mean", chord_test, "Direct imea all-chords mean column."),
        CheckSpec("Remaining 2D chords", "ALLCHORDS_MODE", "remaining2d_truth", "allchords_mode", chord_test, "Direct imea all-chords mode column."),
        CheckSpec("Remaining 2D chords", "ALLCHORDS_STDDEV", "remaining2d_truth", "allchords_std", chord_test, "Direct imea all-chords standard-deviation column."),
    ]


def write_results(rows: list[Row]) -> tuple[int, int]:
    with RESULT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(Row.__dataclass_fields__.keys())
        for row in rows:
            writer.writerow([getattr(row, field) for field in Row.__dataclass_fields__])

    pass_count = sum(row.status == "PASS" for row in rows)
    mismatch_count = sum(row.status == "MISMATCH" for row in rows)
    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["Library", "Version", "Rows", "PASS", "MISMATCH", "Note"])
        writer.writerow(
            [
                "imea",
                imea_version(),
                len(rows),
                pass_count,
                mismatch_count,
                "Only direct columns from shape_measurements_2d are counted; optional raw chord/statistical-length arrays are not used.",
            ]
        )
    return pass_count, mismatch_count


def main(argv: list[str]) -> int:
    repo_root = ensure_nyxus_repo(argv[1] if len(argv) > 1 else None)
    test_data_text = (repo_root / "tests" / "test_data.h").read_text(encoding="utf-8")
    shape_header_text = (repo_root / "tests" / "test_shape_morphology_2d.h").read_text(encoding="utf-8")
    remaining_header_text = (repo_root / "tests" / "test_2d_remaining_features.h").read_text(encoding="utf-8")

    truth_maps = {
        "shape2d_truth": parse_truth_map(shape_header_text, "shape2d_truth"),
        "remaining2d_truth": parse_truth_map(remaining_header_text, "remaining2d_truth"),
    }

    mask = to_binary_image(parse_pixel_array(test_data_text, "shape2d_morphology_mask"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        frame = shape_measurements_2d(mask, spatial_resolution_xy=1, dalpha=9)

    if len(frame) != 1:
        raise ValueError(f"Expected one imea measurement row, got {len(frame)}")

    rows: list[Row] = []
    for spec in check_specs():
        if spec.imea_column not in frame.columns:
            raise ValueError(f"imea did not return expected column {spec.imea_column}")
        expected = truth_maps[spec.truth_map][spec.feature]
        actual = float(frame.loc[frame.index[0], spec.imea_column])
        rows.append(result_row(spec, expected, actual))

    pass_count, mismatch_count = write_results(rows)

    for row in rows:
        print(
            f"{row.family:38s} {row.feature:36s} "
            f"expected={row.expected:.15g} imea={row.actual:.15g} "
            f"absdiff={row.abs_diff:.15g} reldiff={row.rel_diff_pct:.6g}% "
            f"tol={row.tolerance_pct:.2f}% status={row.status}"
        )
    print(
        f"imea independent 2D shape rows: {pass_count} pass, {mismatch_count} mismatch. "
        f"Wrote {RESULT_CSV} and {SUMMARY_CSV}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
