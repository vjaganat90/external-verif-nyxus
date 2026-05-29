#!/usr/bin/env python3
"""Compare Nyxus ZERNIKE2D regression values to mahotas' built-in Zernike moments.

Both the expected Nyxus ZERNIKE2D coefficients and the ROI fixture are parsed
directly from the Nyxus repository so this verifier cannot silently drift from
the C++ test it mirrors:

- expected values  -> tests/test_2d_remaining_features.h
                      (remaining2d_vector_truth["ZERNIKE2D"])
- ROI fixture      -> tests/test_data.h
                      (shape2d_morphology_intensity / shape2d_morphology_mask)

The oracle itself is mahotas.features.zernike_moments; no Nyxus formula is
reimplemented here.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import mahotas as mh
import numpy as np

from nyxus_repo import ensure_nyxus_repo


OUT_DIR = Path(__file__).resolve().parent


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


def parse_zernike_truth(header_text: str) -> np.ndarray:
    """Parse remaining2d_vector_truth["ZERNIKE2D"] from test_2d_remaining_features.h."""
    match = re.search(r'\{\s*"ZERNIKE2D"\s*,\s*\{([\s\S]*?)\}\s*\}', header_text)
    if not match:
        raise ValueError('Could not find remaining2d_vector_truth["ZERNIKE2D"]')
    values = [float(v) for v in re.findall(r"[-+0-9.eE]+", match.group(1))]
    if not values:
        raise ValueError("ZERNIKE2D truth block parsed to zero coefficients")
    return np.array(values, dtype=float)


def to_image(triples: list[tuple[int, int, float]]) -> np.ndarray:
    width = max(x for x, _, _ in triples) + 1
    height = max(y for _, y, _ in triples) + 1
    image = np.zeros((height, width), dtype=float)
    for x, y, value in triples:
        image[y, x] = value
    return image


def fixture_image(test_data_text: str) -> np.ndarray:
    """Build the masked ROI fixture from the shared shape2d morphology arrays."""
    intensity = to_image(parse_pixel_array(test_data_text, "shape2d_morphology_intensity"))
    mask = to_image(parse_pixel_array(test_data_text, "shape2d_morphology_mask")) > 0
    return intensity * mask


def rel_diff_pct(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(expected), np.finfo(float).eps) * 100.0


def main(argv: list[str]) -> int:
    repo_root = ensure_nyxus_repo(argv[1] if len(argv) > 1 else None)
    tests_dir = repo_root / "tests"
    test_data_text = (tests_dir / "test_data.h").read_text(encoding="utf-8")
    remaining_header_text = (tests_dir / "test_2d_remaining_features.h").read_text(encoding="utf-8")

    nyxus_zernike2d = parse_zernike_truth(remaining_header_text)
    image = fixture_image(test_data_text)

    out_path = OUT_DIR / "mahotas_remaining2d_zernike_results.csv"
    summary_path = OUT_DIR / "mahotas_remaining2d_zernike_summary.csv"

    tool_values = mh.features.zernike_moments(image, radius=4, degree=9)
    if len(tool_values) != len(nyxus_zernike2d):
        raise RuntimeError(
            f"Expected {len(nyxus_zernike2d)} mahotas coefficients, got {len(tool_values)}"
        )

    tolerance_pct = 5.0
    rows = []
    pass_count = 0
    for i, (nyxus_expected, tool_value) in enumerate(zip(nyxus_zernike2d, tool_values)):
        absdiff = abs(float(tool_value) - float(nyxus_expected))
        reldiff = rel_diff_pct(float(tool_value), float(nyxus_expected))
        status = "PASS" if reldiff <= tolerance_pct else "MISMATCH"
        pass_count += status == "PASS"
        rows.append(
            {
                "Family": "Remaining 2D Zernike",
                "Feature": f"ZERNIKE2D[{i}]",
                "NyxusExpected": f"{nyxus_expected:.17g}",
                "ToolValue": f"{tool_value:.17g}",
                "AbsDiff": f"{absdiff:.17g}",
                "RelDiffPct": f"{reldiff:.17g}",
                "TolerancePct": f"{tolerance_pct:.17g}",
                "Status": status,
                "Tool": f"mahotas {mh.__version__} features.zernike_moments(radius=4, degree=9)",
                "ImageSource": "tests/test_data.h shape2d_morphology_intensity*shape2d_morphology_mask",
                "Note": "Direct built-in Zernike moments, but normalization/order convention differs from current Nyxus ZERNIKE2D.",
            }
        )
        print(
            f"ZERNIKE2D[{i}] nyxus={nyxus_expected:.15g} tool={tool_value:.15g} "
            f"absdiff={absdiff:.15g} reldiff={reldiff:.6g}% tol={tolerance_pct:.6g}% {status}"
        )

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Feature",
                "Tool",
                "Coefficients",
                "PassCoefficients",
                "MismatchCoefficients",
                "AggregateStatus",
                "Note",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "Feature": "ZERNIKE2D",
                "Tool": f"mahotas {mh.__version__}",
                "Coefficients": len(rows),
                "PassCoefficients": pass_count,
                "MismatchCoefficients": len(rows) - pass_count,
                "AggregateStatus": "PASS" if pass_count == len(rows) else "MISMATCH",
                "Note": "Counted as comparable built-in coverage, not agreement, because at least one coefficient is outside 5%.",
            }
        )

    print(
        f"Remaining2D ZERNIKE2D mahotas rows: {pass_count} pass, "
        f"{len(rows) - pass_count} mismatch. Wrote {out_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
