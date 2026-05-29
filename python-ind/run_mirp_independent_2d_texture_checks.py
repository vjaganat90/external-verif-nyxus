#!/usr/bin/env python3
"""Run MIRP independent checks for Nyxus 2D GLDZM/NGLDM tests."""

from __future__ import annotations

import csv
import importlib.metadata
import logging
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

import mirp


from nyxus_repo import ensure_nyxus_repo

OUT_DIR = Path(__file__).resolve().parent
TOLERANCE_PCT = 5.0
INCOMPARABLE_CSV = OUT_DIR / "mirp_incomparable_no_oracle.csv"


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


GLDZM_MAP = {
    "GLDZM_SDE": "dzm_sde_2d",
    "GLDZM_LDE": "dzm_lde_2d",
    "GLDZM_LGLZE": "dzm_lgze_2d",
    "GLDZM_HGLZE": "dzm_hgze_2d",
    "GLDZM_SDLGLE": "dzm_sdlge_2d",
    "GLDZM_SDHGLE": "dzm_sdhge_2d",
    "GLDZM_LDLGLE": "dzm_ldlge_2d",
    "GLDZM_LDHGLE": "dzm_ldhge_2d",
    "GLDZM_GLNU": "dzm_glnu_2d",
    "GLDZM_GLNUN": "dzm_glnu_norm_2d",
    "GLDZM_ZDNU": "dzm_zdnu_2d",
    "GLDZM_ZDNUN": "dzm_zdnu_norm_2d",
    "GLDZM_ZP": "dzm_z_perc_2d",
    "GLDZM_GLV": "dzm_gl_var_2d",
    "GLDZM_ZDV": "dzm_zd_var_2d",
    "GLDZM_ZDE": "dzm_zd_entr_2d",
}

NGLDM_MAP = {
    "NGLDM_LDE": "ngl_lde_d1_a0.0_2d",
    "NGLDM_HDE": "ngl_hde_d1_a0.0_2d",
    "NGLDM_LGLCE": "ngl_lgce_d1_a0.0_2d",
    "NGLDM_HGLCE": "ngl_hgce_d1_a0.0_2d",
    "NGLDM_LDLGLE": "ngl_ldlge_d1_a0.0_2d",
    "NGLDM_LDHGLE": "ngl_ldhge_d1_a0.0_2d",
    "NGLDM_HDLGLE": "ngl_hdlge_d1_a0.0_2d",
    "NGLDM_HDHGLE": "ngl_hdhge_d1_a0.0_2d",
    "NGLDM_GLNU": "ngl_glnu_d1_a0.0_2d",
    "NGLDM_GLNUN": "ngl_glnu_norm_d1_a0.0_2d",
    "NGLDM_DCNU": "ngl_dcnu_d1_a0.0_2d",
    "NGLDM_DCNUN": "ngl_dcnu_norm_d1_a0.0_2d",
    "NGLDM_DCP": "ngl_dc_perc_d1_a0.0_2d",
    "NGLDM_GLV": "ngl_gl_var_d1_a0.0_2d",
    "NGLDM_DCV": "ngl_dc_var_d1_a0.0_2d",
    "NGLDM_DCENT": "ngl_dc_entr_d1_a0.0_2d",
    "NGLDM_DCENE": "ngl_dc_energy_d1_a0.0_2d",
}


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
    if math.isnan(actual):
        raise ValueError(f"{feature} has no MIRP oracle value; put it in {INCOMPARABLE_CSV.name}")
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


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


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


def parse_unordered_map(header_text: str, map_name: str) -> dict[str, float]:
    pattern = (
        r"static\s+std::unordered_map<std::string,\s*double>\s+"
        + re.escape(map_name)
        + r"\s*\{([\s\S]*?)\};"
    )
    match = re.search(pattern, header_text)
    if not match:
        raise ValueError(f"Could not find {map_name}")
    return {
        key: float(value)
        for key, value in re.findall(r'\{\s*"([^"]+)"\s*,\s*([^}/,\n]+)', match.group(1))
    }


def image_from_triples(
    intensity_triples: Iterable[tuple[int, int, float]],
    mask_triples: Iterable[tuple[int, int, float]],
) -> tuple[np.ndarray, np.ndarray]:
    intensities = {(x, y): value for x, y, value in intensity_triples}
    mask_lookup = {(x, y): value for x, y, value in mask_triples}

    xs = [x for x, _y in mask_lookup]
    ys = [y for _x, y in mask_lookup]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    image = np.zeros((max_y - min_y + 1, max_x - min_x + 1), dtype=float)
    mask = np.zeros_like(image, dtype=np.uint8)

    for (x, y), mask_value in mask_lookup.items():
        row = y - min_y
        col = x - min_x
        if mask_value != 0:
            image[row, col] = intensities[(x, y)]
            mask[row, col] = 1

    return image, mask


def load_ibsi_slices(test_data_text: str) -> list[tuple[np.ndarray, np.ndarray]]:
    slices = []
    for z in range(1, 5):
        slices.append(
            image_from_triples(
                parse_pixel_array(test_data_text, f"ibsi_phantom_z{z}_intensity"),
                parse_pixel_array(test_data_text, f"ibsi_phantom_z{z}_mask"),
            )
        )
    return slices


def mirp_version() -> str:
    return importlib.metadata.version("mirp")


def extract_mirp_feature_table(slices: list[tuple[np.ndarray, np.ndarray]]) -> dict[str, float]:
    values_by_column: dict[str, list[float]] = {}
    extract_kwargs = {
        "base_feature_families": ["gldzm", "ngldm"],
        "base_discretisation_method": "none",
        "by_slice": True,
        "gldzm_spatial_method": "2d",
        "ngldm_spatial_method": "2d",
        "ngldm_distance": 1.0,
        "ngldm_difference_level": 0.0,
        "ibsi_compliant": False,
        "export_features": True,
    }

    logging.disable(logging.CRITICAL)
    try:
        for index, (image, mask) in enumerate(slices, start=1):
            result = mirp.extract_features(
                image=image,
                mask=mask,
                sample_name=f"nyxus_ibsi_z{index}",
                **extract_kwargs,
            )
            if len(result) != 1:
                raise ValueError(f"Expected one MIRP feature table for slice {index}, got {len(result)}")
            table = result[0]
            if len(table.index) != 1:
                raise ValueError(f"Expected one MIRP feature row for slice {index}, got {len(table.index)}")
            row = table.iloc[0]
            for column, value in row.items():
                if isinstance(value, (int, float, np.integer, np.floating)) and not math.isnan(float(value)):
                    values_by_column.setdefault(column, []).append(float(value))
    finally:
        logging.disable(logging.NOTSET)

    return {column: float(np.mean(values)) for column, values in values_by_column.items()}


def add_rows(
    rows: list[Row],
    incomparable: list[dict[str, str]],
    *,
    family: str,
    truth: dict[str, float],
    feature_map: dict[str, str],
    feature_values: dict[str, float],
    image_source: str,
    settings: str,
    note: str,
) -> None:
    for feature, expected in truth.items():
        oracle_feature = feature_map.get(feature)
        if oracle_feature is None or oracle_feature not in feature_values:
            incomparable.append(
                {
                    "Oracle": "MIRP",
                    "Family": family,
                    "Feature": feature,
                    "Reason": f"MIRP {mirp_version()} does not expose this Nyxus row as a direct feature column.",
                }
            )
            continue
        rows.append(
            result_row(
                family,
                feature,
                expected,
                feature_values[oracle_feature],
                f"mirp.{oracle_feature}",
                image_source,
                settings,
                note,
            )
        )


def run(repo_root: Path) -> tuple[list[Row], list[dict[str, str]]]:
    tests_dir = repo_root / "tests"
    test_data_text = read_text(tests_dir / "test_data.h")
    slices = load_ibsi_slices(test_data_text)
    feature_values = extract_mirp_feature_table(slices)

    settings = (
        f"MIRP {mirp_version()}; base_feature_families=['gldzm','ngldm']; "
        "base_discretisation_method='none'; by_slice=True; "
        "gldzm_spatial_method='2d'; ngldm_spatial_method='2d'; "
        "ngldm_distance=1.0; ngldm_difference_level=0.0; ibsi_compliant=False"
    )
    image_source = "IBSI 2D phantom slices z1-z4 with corresponding Nyxus masks"

    rows: list[Row] = []
    incomparable: list[dict[str, str]] = []
    add_rows(
        rows,
        incomparable,
        family="GLDZM",
        truth=parse_unordered_map(read_text(tests_dir / "test_ibsi_gldzm.h"), "ibsi_gldzm_gtruth"),
        feature_map=GLDZM_MAP,
        feature_values=feature_values,
        image_source=image_source,
        settings=settings,
        note="MIRP GLDZM value averaged over the four Nyxus 2D IBSI phantom slices.",
    )
    add_rows(
        rows,
        incomparable,
        family="NGLDM",
        truth=parse_unordered_map(read_text(tests_dir / "test_ibsi_ngldm.h"), "ngldm_regression_gtruth"),
        feature_map=NGLDM_MAP,
        feature_values=feature_values,
        image_source=image_source,
        settings=settings,
        note=(
            "MIRP NGLDM value averaged over the four Nyxus 2D IBSI phantom slices; "
            "compared against the active Nyxus regression constants."
        ),
    )
    return rows, incomparable


def fmt(value: float) -> str:
    if math.isnan(value):
        return ""
    if math.isinf(value):
        return "Inf"
    return f"{value:.15g}"


def write_results(rows: list[Row], out_dir: Path) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / "mirp_independent_2d_texture_results.csv"
    summary_path = out_dir / "mirp_independent_2d_texture_summary.csv"

    with results_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "Family",
                "Feature",
                "Expected",
                "MIRPValue",
                "AbsDiff",
                "RelDiffPct",
                "TolerancePct",
                "Status",
                "OracleFeature",
                "ImageSource",
                "Settings",
                "Note",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.family,
                    row.feature,
                    fmt(row.expected),
                    fmt(row.actual),
                    fmt(row.abs_diff),
                    fmt(row.rel_diff_pct),
                    f"{row.tolerance_pct:.6g}",
                    row.status,
                    row.oracle_feature,
                    row.image_source,
                    row.settings,
                    row.note,
                ]
            )

    summary: dict[str, dict[str, int]] = {}
    for row in rows:
        summary.setdefault(row.family, {"Rows": 0, "PASS": 0, "MISMATCH": 0})
        summary[row.family]["Rows"] += 1
        summary[row.family][row.status] += 1

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Family", "ComparedRows", "PASS", "MISMATCH"])
        for family in sorted(summary):
            item = summary[family]
            writer.writerow([family, item["Rows"], item["PASS"], item["MISMATCH"]])
        writer.writerow(
            [
                "TOTAL",
                len(rows),
                sum(1 for row in rows if row.status == "PASS"),
                sum(1 for row in rows if row.status == "MISMATCH"),
            ]
        )

    return results_path, summary_path


def write_incomparable(rows: list[dict[str, str]]) -> Path:
    with INCOMPARABLE_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Oracle", "Family", "Feature", "Reason"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return INCOMPARABLE_CSV


def print_rows(rows: list[Row]) -> None:
    print(f"MIRP version: {mirp_version()}")
    for row in rows:
        actual = f"{row.actual:.15g}"
        abs_diff = f"{row.abs_diff:.15g}"
        rel = f"{row.rel_diff_pct:.6g}%" if math.isfinite(row.rel_diff_pct) else "Inf"
        print(
            f"{row.family:8s} {row.feature:24s} expected={row.expected:.15g} "
            f"mirp={actual} absdiff={abs_diff} reldiff={rel} "
            f"tol={row.tolerance_pct:.2f}% status={row.status}"
        )


def main(argv: list[str]) -> int:
    repo_root = ensure_nyxus_repo(argv[1] if len(argv) > 1 else None)
    rows, incomparable = run(repo_root)
    print_rows(rows)
    results_path, summary_path = write_results(rows, OUT_DIR)
    incomparable_path = write_incomparable(incomparable)
    pass_count = sum(1 for row in rows if row.status == "PASS")
    mismatch_count = sum(1 for row in rows if row.status == "MISMATCH")
    print(
        f"MIRP independent 2D texture comparisons: {pass_count} pass, "
        f"{mismatch_count} mismatch, {len(incomparable)} incomparable. "
        f"Wrote {results_path}, {summary_path}, and {incomparable_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
