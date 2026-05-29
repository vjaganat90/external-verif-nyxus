#!/usr/bin/env python3
"""Single scikit-image oracle verifier for Nyxus 2D feature tests.

The comparisons are grouped by feature family. This script uses direct
``skimage.measure.regionprops`` properties only. Rows that would require
recomputing Nyxus formulas are kept out of the comparison CSV and written to
``skimage_incomparable_no_oracle.csv``.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import skimage
from skimage.measure import label, regionprops


from nyxus_repo import ensure_nyxus_repo

OUT_DIR = Path(__file__).resolve().parent
RESULT_CSV = OUT_DIR / "skimage_independent_2d_results.csv"
SUMMARY_CSV = OUT_DIR / "skimage_independent_2d_summary.csv"
INCOMPARABLE_CSV = OUT_DIR / "skimage_incomparable_no_oracle.csv"
TOLERANCE_PCT = 5.0

FIELDNAMES = [
    "family",
    "feature",
    "expected",
    "actual",
    "abs_diff",
    "rel_diff_pct",
    "tolerance_pct",
    "status",
    "oracle_tool",
    "oracle_feature",
    "image_source",
    "settings",
    "note",
]

INCOMPARABLE_FIELDNAMES = ["Oracle", "Family", "Feature", "Reason"]

FAMILY_ORDER = {
    "Shape and morphology subset": 10,
    "Added 2D shape morphology extras": 20,
    "Added 2D shape geometric moments": 30,
    "Added 2D intensity geometric moments": 40,
}


def rel_diff_pct(expected: float, actual: float) -> float:
    if math.isnan(expected) or math.isnan(actual):
        return math.nan
    if expected == 0:
        return 0.0 if abs(actual) <= 1e-12 else math.inf
    return abs(expected - actual) / abs(expected) * 100.0


def status(expected: float, actual: float) -> str:
    return "PASS" if rel_diff_pct(expected, actual) <= TOLERANCE_PCT else "MISMATCH"


def make_result_row(
    *,
    family: str,
    feature: str,
    expected: float,
    actual: float | None,
    oracle_feature: str,
    image_source: str,
    settings: str,
    note: str,
) -> dict[str, Any]:
    if actual is None or not math.isfinite(actual):
        raise ValueError(f"{feature} has no scikit-image oracle value; put it in {INCOMPARABLE_CSV.name}")

    row_status = status(expected, actual)
    return {
        "family": family,
        "feature": feature,
        "expected": expected,
        "actual": actual,
        "abs_diff": abs(expected - actual),
        "rel_diff_pct": rel_diff_pct(expected, actual),
        "tolerance_pct": TOLERANCE_PCT,
        "status": row_status,
        "oracle_tool": "scikit-image",
        "oracle_feature": oracle_feature,
        "image_source": image_source,
        "settings": settings,
        "note": note,
    }


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


def parse_geomoment_expectations(nyxus_root: Path) -> list[tuple[str, float]]:
    text = (nyxus_root / "tests" / "test_2d_geometric_moments.h").read_text()
    rows = re.findall(
        r"\{Nyxus::Feature2D::([A-Z0-9_]+),\s*\"[A-Z0-9_]+\",\s*([^}\n]+)\}",
        text,
    )
    return [(feature, float(value.rstrip(","))) for feature, value in rows]


def parse_shape_extra_expectations(nyxus_root: Path) -> list[tuple[str, float]]:
    text = (nyxus_root / "tests" / "test_shape_morphology_2d.h").read_text()
    truth = {
        feature: float(value)
        for feature, value in re.findall(r'\{"([A-Z0-9_]+)",\s*([-+0-9.eE]+)\}', text)
    }
    return [
        (feature, truth[feature])
        for feature in ("GEODETIC_LENGTH", "THICKNESS", "EROSIONS_2_VANISH")
    ]


def to_image(triples: list[tuple[int, int, float]], dtype=float) -> np.ndarray:
    width = max(x for x, _, _ in triples) + 1
    height = max(y for _, y, _ in triples) + 1
    image = np.zeros((height, width), dtype=dtype)
    for x, y, value in triples:
        image[y, x] = value
    return image


def geomoment_fixture() -> tuple[np.ndarray, np.ndarray]:
    height, width = 40, 48
    mask = np.ones((height, width), dtype=np.uint8)
    intensity = np.fromfunction(
        lambda y, x: 10.0 + 3.0 * x + 5.0 * y + ((x * y) % 7),
        (height, width),
        dtype=float,
    )
    return mask, intensity


def pq_index(feature: str) -> tuple[int, int] | None:
    match = re.search(r"_(\d)(\d)$", feature)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def family_for(feature: str) -> str:
    if feature.startswith(("SPAT_", "CENTRAL_", "NORM_", "HU_", "WEIGHTED_", "WT_")):
        return "Added 2D shape geometric moments"
    if feature.startswith("IMOM_"):
        return "Added 2D intensity geometric moments"
    return "Added 2D shape morphology extras"


def na_note(feature: str) -> str:
    if feature.startswith("NORM_SPAT_MOMENT_"):
        return "scikit-image exposes raw spatial moments and normalized central moments, but no direct normalized raw spatial moment property."
    if feature.startswith(("WEIGHTED_", "WT_NORM_CTR_MOM_")):
        return "Nyxus weighted shape moments use contour-distance weighting; considered tools do not expose that as a direct built-in feature."
    if feature.startswith("IMOM_NRM_"):
        return "scikit-image exposes intensity-weighted raw moments, central moments, normalized central moments, and Hu moments, but no direct normalized raw moment property."
    if feature.startswith("IMOM_W"):
        return "Nyxus weighted intensity moments add contour-distance weighting; considered tools do not expose that as a direct built-in feature."
    return "No direct feature-level built-in oracle is exposed by scikit-image for this Nyxus row."


def oracle_for_moment(feature: str, props) -> tuple[float | None, str]:
    pq = pq_index(feature)
    if feature.startswith("SPAT_MOMENT_") and pq:
        p, q = pq
        return float(props.moments[q, p]), "regionprops.moments"
    if feature.startswith("CENTRAL_MOMENT_") and pq:
        p, q = pq
        return float(props.moments_central[q, p]), "regionprops.moments_central"
    if feature.startswith("NORM_CENTRAL_MOMENT_") and pq:
        p, q = pq
        return float(props.moments_normalized[q, p]), "regionprops.moments_normalized"
    if feature.startswith("HU_M"):
        idx = int(feature.split("HU_M", 1)[1]) - 1
        return float(props.moments_hu[idx]), "regionprops.moments_hu"
    if feature.startswith("IMOM_RM_") and pq:
        p, q = pq
        return float(props.moments_weighted[q, p]), "regionprops.moments_weighted"
    if feature.startswith("IMOM_CM_") and pq:
        p, q = pq
        return float(props.moments_weighted_central[q, p]), "regionprops.moments_weighted_central"
    if feature.startswith("IMOM_NCM_") and pq:
        p, q = pq
        return float(props.moments_weighted_normalized[q, p]), "regionprops.moments_weighted_normalized"
    if feature.startswith("IMOM_HU"):
        idx = int(feature.split("IMOM_HU", 1)[1]) - 1
        return float(props.moments_weighted_hu[idx]), "regionprops.moments_weighted_hu"
    return None, ""


def shape_subset_rows(repo_root: Path) -> list[dict[str, Any]]:
    test_data_text = (repo_root / "tests" / "test_data.h").read_text(encoding="utf-8")
    shape_header_text = (repo_root / "tests" / "test_shape_morphology_2d.h").read_text(encoding="utf-8")
    mask = to_image(parse_pixel_array(test_data_text, "shape2d_morphology_mask"), dtype=np.uint8)
    truth = parse_shape_truth(shape_header_text)

    labeled = label(mask > 0, connectivity=2)
    props = regionprops(labeled, spacing=(2.0, 2.0))[0]
    props_unit = regionprops(labeled)[0]
    source = "tests/test_data.h shape2d_morphology_mask"
    settings = "scikit-image regionprops; spacing=(2.0, 2.0) where supported"

    return [
        make_result_row(
            family="Shape and morphology subset",
            feature="DIAMETER_EQUAL_AREA",
            expected=truth["DIAMETER_EQUAL_AREA"],
            actual=float(props.equivalent_diameter_area),
            oracle_feature="skimage.measure.regionprops.equivalent_diameter_area",
            image_source=source,
            settings=settings,
            note="Direct scikit-image equivalent-area diameter measurement on the 2D ROI mask.",
        ),
        make_result_row(
            family="Shape and morphology subset",
            feature="PERIMETER",
            expected=truth["PERIMETER"],
            actual=float(props.perimeter),
            oracle_feature="skimage.measure.regionprops.perimeter",
            image_source=source,
            settings=settings,
            note="Direct scikit-image perimeter measurement on the 2D ROI mask.",
        ),
        make_result_row(
            family="Shape and morphology subset",
            feature="CONVEX_HULL_AREA",
            expected=truth["CONVEX_HULL_AREA"],
            actual=float(props_unit.area_convex),
            oracle_feature="skimage.measure.regionprops.area_convex",
            image_source=source,
            settings="scikit-image regionprops without spacing for pixel-count area_convex.",
            note="Direct scikit-image convex-area measurement on the 2D ROI mask.",
        ),
        make_result_row(
            family="Shape and morphology subset",
            feature="SOLIDITY",
            expected=truth["SOLIDITY"],
            actual=float(props_unit.solidity),
            oracle_feature="skimage.measure.regionprops.solidity",
            image_source=source,
            settings="scikit-image regionprops without spacing for dimensionless solidity.",
            note="Direct scikit-image solidity measurement on the 2D ROI mask.",
        ),
    ]


def added_test_rows(repo_root: Path) -> list[dict[str, Any]]:
    mask, intensity = geomoment_fixture()
    props = regionprops(label(mask), intensity_image=intensity)[0]

    rows: list[dict[str, Any]] = []
    for feature, expected in parse_geomoment_expectations(repo_root):
        actual, metric = oracle_for_moment(feature, props)
        if actual is None:
            continue
        rows.append(
            make_result_row(
                family=family_for(feature),
                feature=feature,
                expected=expected,
                actual=actual,
                oracle_feature=metric,
                image_source="tests/test_2d_geometric_moments.h",
                settings="scikit-image regionprops on the added 48x40 2D fixture.",
                note="Direct scikit-image regionprops property.",
            )
        )

    return rows


def incomparable_rows(repo_root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    mask, intensity = geomoment_fixture()
    props = regionprops(label(mask), intensity_image=intensity)[0]
    for feature, _expected in parse_geomoment_expectations(repo_root):
        actual, _metric = oracle_for_moment(feature, props)
        if actual is None:
            rows.append(
                {
                    "Oracle": "scikit-image",
                    "Family": family_for(feature),
                    "Feature": feature,
                    "Reason": na_note(feature),
                }
            )

    for feature, _expected in parse_shape_extra_expectations(repo_root):
        rows.append(
            {
                "Oracle": "scikit-image",
                "Family": family_for(feature),
                "Feature": feature,
                "Reason": na_note(feature),
            }
        )
    return rows


def family_sort_key(row: dict[str, Any]) -> tuple[int, str]:
    return (FAMILY_ORDER.get(str(row["family"]), 999), str(row["feature"]))


def write_results(rows: list[dict[str, Any]]) -> None:
    with RESULT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(rows: list[dict[str, Any]]) -> None:
    by_family: dict[str, dict[str, int]] = {}
    for row in rows:
        family = str(row["family"])
        by_family.setdefault(
            family,
            {"Rows": 0, "PASS": 0, "PASS_gt_1pct": 0, "MISMATCH": 0},
        )
        item = by_family[family]
        item["Rows"] += 1
        item[str(row["status"])] += 1
        rel = row["rel_diff_pct"]
        if row["status"] == "PASS" and isinstance(rel, float) and rel > 1.0:
            item["PASS_gt_1pct"] += 1

    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["Library", "Version", "Family", "ComparedRows", "PASS", "PASS_gt_1pct", "MISMATCH"])
        for family in sorted(by_family, key=lambda name: FAMILY_ORDER.get(name, 999)):
            item = by_family[family]
            writer.writerow(
                [
                    "scikit-image",
                    skimage.__version__,
                    family,
                    item["Rows"],
                    item["PASS"],
                    item["PASS_gt_1pct"],
                    item["MISMATCH"],
                ]
            )


def write_incomparable(rows: list[dict[str, str]]) -> None:
    with INCOMPARABLE_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=INCOMPARABLE_FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def print_rows(rows: list[dict[str, Any]]) -> None:
    current_family = None
    for row in rows:
        if row["family"] != current_family:
            current_family = row["family"]
            print(f"\n[{current_family}]")
        actual = row["actual"] if row["actual"] != "" else math.nan
        rel = row["rel_diff_pct"] if row["rel_diff_pct"] != "" else math.nan
        abs_diff = row["abs_diff"] if row["abs_diff"] != "" else math.nan
        print(
            f"{row['feature']:42s} expected={float(row['expected']):.15g} "
            f"skimage={float(actual):.15g} absdiff={float(abs_diff):.15g} "
            f"reldiff={float(rel):.6g}% tol={TOLERANCE_PCT:.2f}% status={row['status']}"
        )


def main(argv: list[str]) -> int:
    repo_root = ensure_nyxus_repo(argv[1] if len(argv) > 1 else None)
    rows = sorted(shape_subset_rows(repo_root) + added_test_rows(repo_root), key=family_sort_key)
    no_oracle_rows = sorted(incomparable_rows(repo_root), key=lambda row: (FAMILY_ORDER.get(row["Family"], 999), row["Feature"]))
    write_results(rows)
    write_summary(rows)
    write_incomparable(no_oracle_rows)
    print_rows(rows)

    counts = {name: sum(row["status"] == name for row in rows) for name in ("PASS", "MISMATCH")}
    print(
        f"\nscikit-image {skimage.__version__}: compared_rows={len(rows)} "
        f"PASS={counts['PASS']} MISMATCH={counts['MISMATCH']} "
        f"incomparable_rows={len(no_oracle_rows)}. Wrote {RESULT_CSV}, {SUMMARY_CSV}, and {INCOMPARABLE_CSV}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
