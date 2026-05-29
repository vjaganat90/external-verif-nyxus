#!/usr/bin/env python3
"""Run cp-measure independent checks for Nyxus 2D shape, edge, and neighbor gaps."""

from __future__ import annotations

import csv
import importlib.metadata
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from cp_measure.featurizer import featurize, make_featurizer_config
import cp_measure.multimask.measureobjectneighbors as cp_neighbors


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


def cp_measure_version() -> str:
    return importlib.metadata.version("cp-measure")


def shape_feature_table(test_data_text: str) -> dict[str, float]:
    image2d = to_image(parse_pixel_array(test_data_text, "shape2d_morphology_intensity"), np.float32)
    mask2d = (to_image(parse_pixel_array(test_data_text, "shape2d_morphology_mask"), np.uint16) > 0).astype(np.uint16)
    config = make_featurizer_config(
        channels=["intensity"],
        objects=["roi"],
        intensity=True,
        texture=False,
        granularity=False,
        radial_distribution=False,
        radial_zernikes=False,
        sizeshape=True,
        zernike=False,
        feret=True,
        correlation_pearson=False,
        correlation_costes=False,
        correlation_manders_fold=False,
        correlation_rwc=False,
    )
    data, columns, rows = featurize(image2d[None, :, :], mask2d[None, :, :], config, image_id="shape2d")
    if data.shape[0] != 1 or len(rows) != 1:
        raise ValueError(f"Expected one cp-measure shape row, got data={data.shape}, rows={rows}")
    return {column: float(data[0, index]) for index, column in enumerate(columns)}


def add_shape_rows(rows: list[Row], feature_values: dict[str, float], truth: dict[str, float]) -> None:
    mapping = {
        "MASS_DISPLACEMENT": (
            "Intensity_MassDisplacement__intensity",
            "cp_measure.featurize Intensity_MassDisplacement",
            "Direct CellProfiler-style mass-displacement measurement; differs from Nyxus because Nyxus applies a one-based weighted-centroid compatibility offset.",
        ),
        "COMPACTNESS": (
            "Compactness",
            "cp_measure.featurize Compactness",
            "Direct CellProfiler-style compactness measurement; this is a different shape convention than Nyxus' std(distance-to-centroid)/area statistic.",
        ),
        "EDGE_INTEGRATED_INTENSITY": (
            "Intensity_IntegratedIntensityEdge__intensity",
            "cp_measure.featurize Intensity_IntegratedIntensityEdge",
            "Direct CellProfiler-style edge integrated intensity measurement.",
        ),
        "EDGE_MEAN_INTENSITY": (
            "Intensity_MeanIntensityEdge__intensity",
            "cp_measure.featurize Intensity_MeanIntensityEdge",
            "Direct CellProfiler-style edge mean intensity measurement.",
        ),
        "EDGE_STDDEV_INTENSITY": (
            "Intensity_StdIntensityEdge__intensity",
            "cp_measure.featurize Intensity_StdIntensityEdge",
            "Direct CellProfiler-style edge standard-deviation intensity measurement.",
        ),
        "EDGE_MAX_INTENSITY": (
            "Intensity_MaxIntensityEdge__intensity",
            "cp_measure.featurize Intensity_MaxIntensityEdge",
            "Direct CellProfiler-style edge maximum intensity measurement.",
        ),
        "EDGE_MIN_INTENSITY": (
            "Intensity_MinIntensityEdge__intensity",
            "cp_measure.featurize Intensity_MinIntensityEdge",
            "Direct CellProfiler-style edge minimum intensity measurement.",
        ),
        "CONVEX_HULL_AREA": (
            "ConvexArea",
            "cp_measure.featurize ConvexArea",
            "Direct CellProfiler-style convex-area measurement.",
        ),
        "SOLIDITY": (
            "Solidity",
            "cp_measure.featurize Solidity",
            "Direct CellProfiler-style solidity measurement.",
        ),
    }
    for feature, (column, oracle_feature, note) in mapping.items():
        rows.append(
            result_row(
                family="Shape and morphology subset",
                feature=feature,
                expected=truth[feature],
                actual=feature_values[column],
                oracle_feature=oracle_feature,
                image_source="tests/test_data.h shape2d_morphology_intensity + shape2d_morphology_mask",
                settings="cp-measure single-channel 2D image; intensity and sizeshape feature groups enabled",
                note=note,
            )
        )


def add_neighbor_rows(rows: list[Row], test_data_text: str, truth: dict[int, dict[str, float]]) -> None:
    labels = to_image(parse_pixel_array(test_data_text, "neighborhood2d_scene_labels"), np.uint16)
    measurements = cp_neighbors.measureobjectneighbors(
        labels,
        labels,
        distance_method=cp_neighbors.D_ADJACENT,
        distance=1,
    )
    mapping = {
        "NUM_NEIGHBORS": (
            "Neighbors_NumberOfNeighbors_Adjacent",
            "cp_measure.measureobjectneighbors NumberOfNeighbors_Adjacent",
            "Direct CellProfiler-style adjacent-neighbor count. This supersedes the earlier Centrosome mismatch for the tracker because it matches Nyxus' fixture convention.",
        ),
        "PERCENT_TOUCHING": (
            "Neighbors_PercentTouching_Adjacent",
            "cp_measure.measureobjectneighbors PercentTouching_Adjacent",
            "Direct CellProfiler-style adjacent percent-touching measurement; convention differs from Nyxus on this fixture.",
        ),
        "CLOSEST_NEIGHBOR1_DIST": (
            "Neighbors_FirstClosestDistance_Adjacent",
            "cp_measure.measureobjectneighbors FirstClosestDistance_Adjacent",
            "Direct CellProfiler-style first closest-neighbor centroid distance.",
        ),
        "CLOSEST_NEIGHBOR2_DIST": (
            "Neighbors_SecondClosestDistance_Adjacent",
            "cp_measure.measureobjectneighbors SecondClosestDistance_Adjacent",
            "Direct CellProfiler-style second closest-neighbor distance. cp-measure reports a second nearest object even for labels with one adjacent neighbor; Nyxus stores 0 for those labels.",
        ),
    }
    for feature, (column, oracle_feature, note) in mapping.items():
        values = measurements[column]
        for label in sorted(truth):
            rows.append(
                result_row(
                    family="Neighbor / touching features",
                    feature=f"{feature}[label={label}]",
                    expected=truth[label][feature],
                    actual=float(values[label - 1]),
                    oracle_feature=oracle_feature,
                    image_source="tests/test_data.h neighborhood2d_scene_labels",
                    settings="cp-measure MeasureObjectNeighbors equivalent; distance_method=Adjacent",
                    note=note,
                )
            )


def feature_key(feature: str) -> str:
    return feature.split("[", 1)[0]


def write_results(rows: list[Row]) -> tuple[int, int]:
    result_path = OUT_DIR / "cp_measure_independent_2d_shape_neighbor_results.csv"
    summary_path = OUT_DIR / "cp_measure_independent_2d_shape_neighbor_summary.csv"

    with result_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(Row.__dataclass_fields__.keys())
        for row in rows:
            writer.writerow([getattr(row, field) for field in Row.__dataclass_fields__])

    feature_statuses = {}
    for key in sorted({feature_key(row.feature) for row in rows}):
        grouped = [row for row in rows if feature_key(row.feature) == key]
        feature_statuses[key] = "PASS" if all(row.status == "PASS" for row in grouped) else "MISMATCH"
    feature_pass = sum(status == "PASS" for status in feature_statuses.values())
    feature_mismatch = sum(status == "MISMATCH" for status in feature_statuses.values())

    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["Library", "Version", "FeatureRows", "DetailRows", "PASS", "MISMATCH", "Note"])
        writer.writerow(
            [
                "cp-measure",
                cp_measure_version(),
                len(feature_statuses),
                len(rows),
                feature_pass,
                feature_mismatch,
                "Feature-level PASS requires all detail rows for that feature to pass.",
            ]
        )
    return feature_pass, feature_mismatch


def main(argv: list[str]) -> int:
    repo_root = ensure_nyxus_repo(argv[1] if len(argv) > 1 else None)
    test_data_text = (repo_root / "tests" / "test_data.h").read_text(encoding="utf-8")
    shape_truth = parse_shape_truth((repo_root / "tests" / "test_shape_morphology_2d.h").read_text(encoding="utf-8"))
    neighbor_truth = parse_neighbor_truth((repo_root / "tests" / "test_neighbors_2d.h").read_text(encoding="utf-8"))

    rows: list[Row] = []
    add_shape_rows(rows, shape_feature_table(test_data_text), shape_truth)
    add_neighbor_rows(rows, test_data_text, neighbor_truth)
    feature_pass, feature_mismatch = write_results(rows)

    for row in rows:
        print(
            f"{row.family:30s} {row.feature:42s} "
            f"expected={row.expected:.15g} cp_measure={row.actual:.15g} "
            f"absdiff={row.abs_diff:.15g} reldiff={row.rel_diff_pct:.6g}% "
            f"tol={row.tolerance_pct:.2f}% status={row.status}"
        )
    print(
        "cp-measure independent 2D feature rows: "
        f"{feature_pass} pass, {feature_mismatch} mismatch. "
        f"Wrote {OUT_DIR / 'cp_measure_independent_2d_shape_neighbor_results.csv'} and "
        f"{OUT_DIR / 'cp_measure_independent_2d_shape_neighbor_summary.csv'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
