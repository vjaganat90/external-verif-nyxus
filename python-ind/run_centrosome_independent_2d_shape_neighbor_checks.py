#!/usr/bin/env python3
"""Run Centrosome independent checks for Nyxus 2D shape/neighbor gaps."""

from __future__ import annotations

import csv
import importlib.metadata
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import centrosome.cpmorphology as cpm
import numpy as np


from nyxus_repo import ensure_nyxus_repo

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


def result_row(
    family: str,
    feature: str,
    expected: float,
    actual: float,
    oracle_feature: str,
    image_source: str,
    settings: str,
    note: str,
) -> Row:
    abs_diff = abs(actual - expected) if not math.isnan(actual) else math.nan
    rel_diff_pct = relative_diff_pct(actual, expected)
    status = "PASS" if rel_diff_pct <= TOLERANCE_PCT else "MISMATCH"
    return Row(
        family=family,
        feature=feature,
        expected=expected,
        actual=actual,
        abs_diff=abs_diff,
        rel_diff_pct=rel_diff_pct,
        tolerance_pct=TOLERANCE_PCT,
        status=status,
        oracle_feature=oracle_feature,
        image_source=image_source,
        settings=settings,
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


def parse_neighbor_truth(header_text: str) -> dict[int, dict[str, float]]:
    truth: dict[int, dict[str, float]] = {}
    current_label: int | None = None
    for line in header_text.splitlines():
        label_match = re.match(r"\s*\{(\d+),\s*\{", line)
        if label_match:
            current_label = int(label_match.group(1))
            truth[current_label] = {}
            continue
        if current_label is None:
            continue
        value_match = re.search(r'\{\s*"([^"]+)"\s*,\s*([0-9.+\-eE]+)\s*\}', line)
        if value_match:
            truth[current_label][value_match.group(1)] = float(value_match.group(2))
        if "}}" in line:
            current_label = None
    return truth


def to_image(triples: list[tuple[int, int, float]], dtype: type[np.integer] | type[np.floating]) -> np.ndarray:
    width = max(x for x, _, _ in triples) + 1
    height = max(y for _, y, _ in triples) + 1
    image = np.zeros((height, width), dtype=dtype)
    for x, y, value in triples:
        image[y, x] = value
    return image


def centrosome_version() -> str:
    return importlib.metadata.version("centrosome")


def check_minimum_enclosing_circle(test_data_text: str, shape_truth: dict[str, float]) -> Row:
    mask = to_image(parse_pixel_array(test_data_text, "shape2d_morphology_mask"), np.int32)
    labels = (mask > 0).astype(np.int32)
    _centers, radii = cpm.minimum_enclosing_circle(labels, np.array([1], dtype=np.int32))
    actual = float(2.0 * radii[0])
    expected = shape_truth["DIAMETER_MIN_ENCLOSING_CIRCLE"]
    return result_row(
        family="Shape and morphology subset",
        feature="DIAMETER_MIN_ENCLOSING_CIRCLE",
        expected=expected,
        actual=actual,
        oracle_feature="centrosome.cpmorphology.minimum_enclosing_circle",
        image_source="tests/test_data.h shape2d_morphology_mask",
        settings="Centrosome 2D label matrix; diameter = 2 * returned radius",
        note="Direct Centrosome minimum-enclosing-circle oracle on the 2D ROI mask.",
    )


def check_num_neighbors(test_data_text: str, neighbor_truth: dict[int, dict[str, float]]) -> list[Row]:
    labels = to_image(parse_pixel_array(test_data_text, "neighborhood2d_scene_labels"), np.int32)
    counts, _indexes, _neighbors = cpm.find_neighbors(labels)
    rows: list[Row] = []
    for label in sorted(neighbor_truth):
        expected = neighbor_truth[label]["NUM_NEIGHBORS"]
        actual = float(counts[label - 1])
        rows.append(
            result_row(
                family="Neighbor / touching features",
                feature=f"NUM_NEIGHBORS[label={label}]",
                expected=expected,
                actual=actual,
                oracle_feature="centrosome.cpmorphology.find_neighbors",
                image_source="tests/test_data.h neighborhood2d_scene_labels",
                settings="Centrosome 8-connected touching adjacency on the 2D label matrix",
                note=(
                    "Centrosome uses 8-connected adjacency. This agrees for the central label but counts "
                    "diagonal neighbor-neighbor contacts that Nyxus does not count in this fixture."
                ),
            )
        )
    return rows


def write_results(rows: list[Row]) -> None:
    result_path = OUT_DIR / "centrosome_independent_2d_shape_neighbor_results.csv"
    summary_path = OUT_DIR / "centrosome_independent_2d_shape_neighbor_summary.csv"

    with result_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(Row.__dataclass_fields__.keys())
        for row in rows:
            writer.writerow([getattr(row, field) for field in Row.__dataclass_fields__])

    min_circle_pass = next(row.status == "PASS" for row in rows if row.feature == "DIAMETER_MIN_ENCLOSING_CIRCLE")
    neighbor_rows = [row for row in rows if row.feature.startswith("NUM_NEIGHBORS")]
    neighbor_pass = all(row.status == "PASS" for row in neighbor_rows)
    feature_pass = int(min_circle_pass) + int(neighbor_pass)
    feature_mismatch = 2 - feature_pass

    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["Library", "Version", "FeatureRows", "DetailRows", "PASS", "MISMATCH", "Note"])
        writer.writerow(
            [
                "Centrosome",
                centrosome_version(),
                2,
                len(rows),
                feature_pass,
                feature_mismatch,
                "Feature-level PASS requires all detail rows for that feature to pass.",
            ]
        )


def main(argv: list[str]) -> int:
    repo_root = ensure_nyxus_repo(argv[1] if len(argv) > 1 else None)
    test_data_text = (repo_root / "tests" / "test_data.h").read_text(encoding="utf-8")
    shape_header_text = (repo_root / "tests" / "test_shape_morphology_2d.h").read_text(encoding="utf-8")
    neighbor_header_text = (repo_root / "tests" / "test_neighbors_2d.h").read_text(encoding="utf-8")

    rows = [
        check_minimum_enclosing_circle(test_data_text, parse_shape_truth(shape_header_text)),
        *check_num_neighbors(test_data_text, parse_neighbor_truth(neighbor_header_text)),
    ]
    write_results(rows)

    for row in rows:
        print(
            f"{row.family:30s} {row.feature:42s} "
            f"expected={row.expected:.15g} centrosome={row.actual:.15g} "
            f"absdiff={row.abs_diff:.15g} reldiff={row.rel_diff_pct:.6g}% "
            f"tol={row.tolerance_pct:.2f}% status={row.status}"
        )

    feature_pass = 0
    feature_mismatch = 0
    for feature in ["DIAMETER_MIN_ENCLOSING_CIRCLE", "NUM_NEIGHBORS"]:
        feature_rows = [row for row in rows if row.feature == feature or row.feature.startswith(feature + "[")]
        if all(row.status == "PASS" for row in feature_rows):
            feature_pass += 1
        else:
            feature_mismatch += 1

    print(
        "Centrosome independent 2D feature rows: "
        f"{feature_pass} pass, {feature_mismatch} mismatch. "
        f"Wrote {OUT_DIR / 'centrosome_independent_2d_shape_neighbor_results.csv'} and "
        f"{OUT_DIR / 'centrosome_independent_2d_shape_neighbor_summary.csv'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
