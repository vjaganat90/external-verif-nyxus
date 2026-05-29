#!/usr/bin/env python3
"""Run PyRadiomics-based independent checks for Nyxus 2D feature tests."""

from __future__ import annotations

import csv
import ast
import operator
import math
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

# PyRadiomics 3.0.1 uses the removed alias np.NaN in a few places. Supplying
# it here keeps the oracle runnable on the Python 3.12 / NumPy 2.x nyxus_ci env.
if not hasattr(np, "NaN"):
    np.NaN = np.nan

import SimpleITK as sitk
import radiomics
from radiomics import firstorder, glcm, gldm, glrlm, glszm, ngtdm

radiomics.setVerbosity(logging.ERROR)


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
        for key, value in re.findall(r'\{"([^"]+)",\s*([^}]+)\}', match.group(1))
    }


def parse_first_order_truth(header_text: str) -> dict[str, float]:
    truth: dict[str, float] = {}
    pattern = r"Feature2D::([A-Z0-9_]+)\]\[0\],\s*([^,\)]+)(?:,\s*([0-9.]+))?\)\);"
    for feature, raw_value, _frac in re.findall(pattern, header_text):
        truth[feature] = parse_numeric_expr(raw_value)
    return truth


def parse_numeric_expr(raw_value: str) -> float:
    node = ast.parse(raw_value.strip(), mode="eval").body
    ops = {
        ast.Add: operator.add,
        ast.Sub: operator.sub,
        ast.Mult: operator.mul,
        ast.Div: operator.truediv,
        ast.USub: operator.neg,
        ast.UAdd: operator.pos,
    }

    def eval_node(item):
        if isinstance(item, ast.Constant) and isinstance(item.value, (int, float)):
            return float(item.value)
        if isinstance(item, ast.UnaryOp) and type(item.op) in ops:
            return ops[type(item.op)](eval_node(item.operand))
        if isinstance(item, ast.BinOp) and type(item.op) in ops:
            return ops[type(item.op)](eval_node(item.left), eval_node(item.right))
        raise ValueError(f"Unsupported numeric expression: {raw_value}")

    return float(eval_node(node))


def image_from_triples(
    intensity_triples: Iterable[tuple[int, int, float]],
    mask_triples: Iterable[tuple[int, int, float]] | None = None,
    *,
    include_aabb_background: bool,
) -> tuple[np.ndarray, np.ndarray]:
    intensities = {(x, y): value for x, y, value in intensity_triples}
    if mask_triples is None:
        mask_lookup = {xy: 1.0 for xy in intensities}
    else:
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

    if include_aabb_background:
        mask[:, :] = 1

    return image, mask


def nyxus_matlab_binning(image: np.ndarray, levels: int) -> np.ndarray:
    max_i = float(np.nanmax(image))
    if max_i <= 0:
        raise ValueError("Expected positive max intensity for texture binning")
    out = np.ones_like(image, dtype=float)
    positive = image > 0
    out[positive] = np.floor(levels / max_i * image[positive] + 1.0)
    out = np.clip(out, 1, levels)
    return out


def to_sitk(image2d: np.ndarray, mask2d: np.ndarray) -> tuple[sitk.Image, sitk.Image]:
    image3d = image2d[np.newaxis, :, :].astype(np.float32)
    mask3d = mask2d[np.newaxis, :, :].astype(np.uint32)
    image = sitk.GetImageFromArray(image3d)
    mask = sitk.GetImageFromArray(mask3d)
    image.SetSpacing((1.0, 1.0, 1.0))
    mask.CopyInformation(image)
    return image, mask


def execute_feature_class(
    cls,
    image2d: np.ndarray,
    mask2d: np.ndarray,
    *,
    settings: dict,
) -> dict[str, float]:
    image, mask = to_sitk(image2d, mask2d)
    calculator = cls(image, mask, **settings)
    calculator.enableAllFeatures()
    values = calculator.execute()
    return {key: float(value) for key, value in values.items() if not key.startswith("diagnostics_")}


def average_slice_feature(
    slices: Iterable[tuple[np.ndarray, np.ndarray]],
    cls,
    pyradiomics_feature: str,
    *,
    settings: dict,
) -> float:
    values = []
    for image2d, mask2d in slices:
        values.append(execute_feature_class(cls, image2d, mask2d, settings=settings)[pyradiomics_feature])
    return float(np.mean(values))


def load_ibsi_slices(test_data_text: str, *, levels: int) -> list[tuple[np.ndarray, np.ndarray]]:
    slices = []
    for z in range(1, 5):
        image, mask = image_from_triples(
            parse_pixel_array(test_data_text, f"ibsi_phantom_z{z}_intensity"),
            parse_pixel_array(test_data_text, f"ibsi_phantom_z{z}_mask"),
            include_aabb_background=True,
        )
        slices.append((nyxus_matlab_binning(image, levels), mask))
    return slices


def load_cat2500(test_data_text: str) -> tuple[np.ndarray, np.ndarray]:
    image, mask = image_from_triples(
        parse_pixel_array(test_data_text, "cat2500_int"),
        parse_pixel_array(test_data_text, "cat2500_seg"),
        include_aabb_background=False,
    )
    return image, mask


def load_first_order_image(test_data_text: str) -> tuple[np.ndarray, np.ndarray]:
    return image_from_triples(
        parse_pixel_array(test_data_text, "pixelIntensityFeaturesTestData"),
        None,
        include_aabb_background=False,
    )


FIRST_ORDER_CHECKS = {
    "ENTROPY": ("Entropy", {"binCount": 9, "voxelArrayShift": 0}),
    "UNIFORMITY": ("Uniformity", {"binCount": 256, "voxelArrayShift": 0}),
    "ROBUST_MEAN_ABSOLUTE_DEVIATION": (
        "RobustMeanAbsoluteDeviation",
        {"binWidth": 25, "voxelArrayShift": 0},
    ),
}


GLCM_MAP = {
    "GLCM_ACOR": "Autocorrelation",
    "GLCM_ASM": "JointEnergy",
    "GLCM_CLUPROM": "ClusterProminence",
    "GLCM_CLUSHADE": "ClusterShade",
    "GLCM_CLUTEND": "ClusterTendency",
    "GLCM_CONTRAST": "Contrast",
    "GLCM_CORRELATION": "Correlation",
    "GLCM_DIFAVE": "DifferenceAverage",
    "GLCM_DIFENTRO": "DifferenceEntropy",
    "GLCM_DIFVAR": "DifferenceVariance",
    "GLCM_DIS": "DifferenceAverage",
    "GLCM_ENERGY": "JointEnergy",
    "GLCM_ENTROPY": "JointEntropy",
    "GLCM_HOM1": "Id",
    "GLCM_HOM2": "Idm",
    "GLCM_ID": "Id",
    "GLCM_IDN": "Idn",
    "GLCM_IDM": "Idm",
    "GLCM_IDMN": "Idmn",
    "GLCM_INFOMEAS1": "Imc1",
    "GLCM_INFOMEAS2": "Imc2",
    "GLCM_IV": "InverseVariance",
    "GLCM_JAVE": "JointAverage",
    "GLCM_JE": "JointEntropy",
    "GLCM_JMAX": "MaximumProbability",
    "GLCM_JVAR": "SumSquares",
    "GLCM_SUMAVERAGE": "SumAverage",
    "GLCM_SUMENTROPY": "SumEntropy",
    "GLCM_SUMVARIANCE": "ClusterTendency",
    "GLCM_VARIANCE": "SumSquares",
}


GLCM_AVERAGED_ROWS = [
    "GLCM_ACOR",
    "GLCM_ASM",
    "GLCM_CLUPROM",
    "GLCM_CLUSHADE",
    "GLCM_CLUTEND",
    "GLCM_CONTRAST",
    "GLCM_CORRELATION",
    "GLCM_DIFAVE",
    "GLCM_DIFENTRO",
    "GLCM_DIFVAR",
    "GLCM_DIS",
    "GLCM_ENERGY",
    "GLCM_ENTROPY",
    "GLCM_HOM1",
    "GLCM_ID",
    "GLCM_IDN",
    "GLCM_IDM",
    "GLCM_IDMN",
    "GLCM_INFOMEAS1",
    "GLCM_INFOMEAS2",
    "GLCM_IV",
    "GLCM_JAVE",
    "GLCM_JE",
    "GLCM_JMAX",
    "GLCM_JVAR",
    "GLCM_SUMAVERAGE",
    "GLCM_SUMENTROPY",
    "GLCM_SUMVARIANCE",
    "GLCM_VARIANCE",
]


GLRLM_MAP = {
    "GLRLM_SRE": "ShortRunEmphasis",
    "GLRLM_LRE": "LongRunEmphasis",
    "GLRLM_LGLRE": "LowGrayLevelRunEmphasis",
    "GLRLM_HGLRE": "HighGrayLevelRunEmphasis",
    "GLRLM_SRLGLE": "ShortRunLowGrayLevelEmphasis",
    "GLRLM_SRHGLE": "ShortRunHighGrayLevelEmphasis",
    "GLRLM_LRLGLE": "LongRunLowGrayLevelEmphasis",
    "GLRLM_LRHGLE": "LongRunHighGrayLevelEmphasis",
    "GLRLM_GLN": "GrayLevelNonUniformity",
    "GLRLM_GLNN": "GrayLevelNonUniformityNormalized",
    "GLRLM_RLN": "RunLengthNonUniformity",
    "GLRLM_RLNN": "RunLengthNonUniformityNormalized",
    "GLRLM_RP": "RunPercentage",
    "GLRLM_GLV": "GrayLevelVariance",
    "GLRLM_RV": "RunVariance",
    "GLRLM_RE": "RunEntropy",
}


GLSZM_MAP = {
    "GLSZM_SAE": "SmallAreaEmphasis",
    "GLSZM_LAE": "LargeAreaEmphasis",
    "GLSZM_LGLZE": "LowGrayLevelZoneEmphasis",
    "GLSZM_HGLZE": "HighGrayLevelZoneEmphasis",
    "GLSZM_SALGLE": "SmallAreaLowGrayLevelEmphasis",
    "GLSZM_SAHGLE": "SmallAreaHighGrayLevelEmphasis",
    "GLSZM_LALGLE": "LargeAreaLowGrayLevelEmphasis",
    "GLSZM_LAHGLE": "LargeAreaHighGrayLevelEmphasis",
    "GLSZM_GLN": "GrayLevelNonUniformity",
    "GLSZM_GLNN": "GrayLevelNonUniformityNormalized",
    "GLSZM_SZN": "SizeZoneNonUniformity",
    "GLSZM_SZNN": "SizeZoneNonUniformityNormalized",
    "GLSZM_ZP": "ZonePercentage",
    "GLSZM_GLV": "GrayLevelVariance",
    "GLSZM_ZV": "ZoneVariance",
    "GLSZM_ZE": "ZoneEntropy",
}


GLDM_MAP = {
    "GLDM_SDE": "SmallDependenceEmphasis",
    "GLDM_LDE": "LargeDependenceEmphasis",
    "GLDM_LGLE": "LowGrayLevelEmphasis",
    "GLDM_HGLE": "HighGrayLevelEmphasis",
    "GLDM_SDLGLE": "SmallDependenceLowGrayLevelEmphasis",
    "GLDM_SDHGLE": "SmallDependenceHighGrayLevelEmphasis",
    "GLDM_LDLGLE": "LargeDependenceLowGrayLevelEmphasis",
    "GLDM_LDHGLE": "LargeDependenceHighGrayLevelEmphasis",
    "GLDM_GLN": "GrayLevelNonUniformity",
    "GLDM_DN": "DependenceNonUniformity",
    "GLDM_DNN": "DependenceNonUniformityNormalized",
    "GLDM_GLV": "GrayLevelVariance",
    "GLDM_DV": "DependenceVariance",
    "GLDM_DE": "DependenceEntropy",
}


NGTDM_MAP = {
    "NGTDM_COARSENESS": "Coarseness",
    "NGTDM_CONTRAST": "Contrast",
    "NGTDM_BUSYNESS": "Busyness",
    "NGTDM_COMPLEXITY": "Complexity",
    "NGTDM_STRENGTH": "Strength",
}


def add_rows_from_map(
    rows: list[Row],
    *,
    family: str,
    truth: dict[str, float],
    feature_map: dict[str, str],
    value_func: Callable[[str], float],
    oracle_prefix: str,
    image_source: str,
    settings: str,
    note: str,
    duplicate_ave_rows: bool = False,
) -> None:
    for feature, oracle_feature in feature_map.items():
        if feature in truth:
            rows.append(
                result_row(
                    family,
                    feature,
                    truth[feature],
                    value_func(oracle_feature),
                    f"{oracle_prefix}.{oracle_feature}",
                    image_source,
                    settings,
                    note,
                )
            )
        if duplicate_ave_rows:
            ave_feature = f"{feature}_AVE"
            if feature in truth:
                rows.append(
                    result_row(
                        family,
                        ave_feature,
                        truth[feature],
                        value_func(oracle_feature),
                        f"{oracle_prefix}.{oracle_feature}",
                        image_source,
                        settings,
                        "Nyxus _AVE row uses the same asserted test value as the base feature.",
                    )
                )


def run(repo_root: Path) -> list[Row]:
    tests_dir = repo_root / "tests"
    test_data_text = read_text(tests_dir / "test_data.h")
    rows: list[Row] = []

    first_truth = parse_first_order_truth(read_text(tests_dir / "test_pixel_intensity_features.h"))
    first_image, first_mask = load_first_order_image(test_data_text)
    for feature, (oracle_feature, settings) in FIRST_ORDER_CHECKS.items():
        values = execute_feature_class(firstorder.RadiomicsFirstOrder, first_image, first_mask, settings=settings)
        rows.append(
            result_row(
                "First-order intensity",
                feature,
                first_truth[feature],
                values[oracle_feature],
                f"RadiomicsFirstOrder.{oracle_feature}",
                "pixelIntensityFeaturesTestData ROI pixels",
                repr(settings),
                "PyRadiomics first-order feature with Nyxus-documented histogram setting where applicable.",
            )
        )

    texture_settings = {
        "binWidth": 1,
        "force2D": True,
        "force2Ddimension": 0,
        "distances": [1],
        "symmetricalGLCM": False,
    }
    texture_settings_text = repr(texture_settings)
    ibsi_100_slices = load_ibsi_slices(test_data_text, levels=100)

    glcm_truth = parse_unordered_map(read_text(tests_dir / "test_glcm.h"), "glcm_values")
    glcm_values = {
        feature: average_slice_feature(
            ibsi_100_slices,
            glcm.RadiomicsGLCM,
            feature,
            settings=texture_settings,
        )
        for feature in set(GLCM_MAP.values())
    }
    add_rows_from_map(
        rows,
        family="GLCM",
        truth=glcm_truth,
        feature_map=GLCM_MAP,
        value_func=lambda feature: glcm_values[feature],
        oracle_prefix="RadiomicsGLCM",
        image_source="IBSI 2D phantom slices, Nyxus MATLAB-style 100-level binning, AABB background included",
        settings=texture_settings_text,
        note="PyRadiomics GLCM value averaged over the four Nyxus 2D phantom slices.",
        duplicate_ave_rows=False,
    )
    for base_feature in GLCM_AVERAGED_ROWS:
        if base_feature in glcm_truth and base_feature in GLCM_MAP:
            oracle_feature = GLCM_MAP[base_feature]
            rows.append(
                result_row(
                    "GLCM",
                    f"{base_feature}_AVE",
                    glcm_truth[base_feature],
                    glcm_values[oracle_feature],
                    f"RadiomicsGLCM.{oracle_feature}",
                    "IBSI 2D phantom slices, Nyxus MATLAB-style 100-level binning, AABB background included",
                    texture_settings_text,
                    "Nyxus _AVE row uses the same asserted test value as the base feature.",
                )
            )

    glrlm_truth = parse_unordered_map(read_text(tests_dir / "test_glrlm.h"), "glrlm_values")
    glrlm_values = {
        feature: average_slice_feature(
            ibsi_100_slices,
            glrlm.RadiomicsGLRLM,
            feature,
            settings=texture_settings,
        )
        for feature in set(GLRLM_MAP.values())
    }
    add_rows_from_map(
        rows,
        family="GLRLM",
        truth=glrlm_truth,
        feature_map=GLRLM_MAP,
        value_func=lambda feature: glrlm_values[feature],
        oracle_prefix="RadiomicsGLRLM",
        image_source="IBSI 2D phantom slices, Nyxus MATLAB-style 100-level binning, AABB background included",
        settings=texture_settings_text,
        note="PyRadiomics GLRLM value averaged over the four Nyxus 2D phantom slices.",
        duplicate_ave_rows=True,
    )

    glszm_truth = parse_unordered_map(read_text(tests_dir / "test_glszm.h"), "glszm_values")
    glszm_values = {
        feature: average_slice_feature(
            ibsi_100_slices,
            glszm.RadiomicsGLSZM,
            feature,
            settings=texture_settings,
        )
        for feature in set(GLSZM_MAP.values())
    }
    add_rows_from_map(
        rows,
        family="GLSZM",
        truth=glszm_truth,
        feature_map=GLSZM_MAP,
        value_func=lambda feature: glszm_values[feature],
        oracle_prefix="RadiomicsGLSZM",
        image_source="IBSI 2D phantom slices, Nyxus MATLAB-style 100-level binning, AABB background included",
        settings=texture_settings_text,
        note="PyRadiomics GLSZM value averaged over the four Nyxus 2D phantom slices.",
    )

    gldm_truth = parse_unordered_map(read_text(tests_dir / "test_gldm.h"), "gldm_values")
    cat_image, cat_mask = load_cat2500(test_data_text)
    gldm_settings = {
        "binCount": 128,
        "force2D": True,
        "force2Ddimension": 0,
        "distances": [1],
    }
    gldm_values = execute_feature_class(gldm.RadiomicsGLDM, cat_image, cat_mask, settings=gldm_settings)
    add_rows_from_map(
        rows,
        family="GLDM",
        truth=gldm_truth,
        feature_map=GLDM_MAP,
        value_func=lambda feature: gldm_values[feature],
        oracle_prefix="RadiomicsGLDM",
        image_source="cat2500_int/cat2500_seg ROI",
        settings=repr(gldm_settings),
        note="PyRadiomics GLDM on the Nyxus cat2500 2D test ROI.",
    )

    ngtdm_truth = parse_unordered_map(read_text(tests_dir / "test_ngtdm.h"), "ngtdm_values")
    ngtdm_values = {
        feature: average_slice_feature(
            ibsi_100_slices,
            ngtdm.RadiomicsNGTDM,
            feature,
            settings=texture_settings,
        )
        for feature in set(NGTDM_MAP.values())
    }
    add_rows_from_map(
        rows,
        family="NGTDM",
        truth=ngtdm_truth,
        feature_map=NGTDM_MAP,
        value_func=lambda feature: ngtdm_values[feature],
        oracle_prefix="RadiomicsNGTDM",
        image_source="IBSI 2D phantom slices, Nyxus MATLAB-style 100-level binning, AABB background included",
        settings=texture_settings_text,
        note="PyRadiomics NGTDM value averaged over the four Nyxus 2D phantom slices.",
    )

    return rows


def write_results(rows: list[Row], out_dir: Path) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / "pyradiomics_independent_2d_results.csv"
    summary_path = out_dir / "pyradiomics_independent_2d_summary.csv"

    with results_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "Family",
                "Feature",
                "Expected",
                "PyRadiomicsValue",
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
                    f"{row.expected:.15g}",
                    f"{row.actual:.15g}",
                    f"{row.abs_diff:.15g}" if not math.isnan(row.abs_diff) else "",
                    f"{row.rel_diff_pct:.12g}" if math.isfinite(row.rel_diff_pct) else "Inf",
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
        writer.writerow(["Family", "Rows", "PASS", "MISMATCH"])
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


def print_rows(rows: list[Row]) -> None:
    print(f"PyRadiomics version: {radiomics.__version__}")
    for row in rows:
        actual = f"{row.actual:.15g}" if not math.isnan(row.actual) else "nan"
        rel = f"{row.rel_diff_pct:.6g}%" if math.isfinite(row.rel_diff_pct) else "Inf"
        print(
            f"{row.family:22s} {row.feature:42s} expected={row.expected:.15g} "
            f"pyradiomics={actual} absdiff={row.abs_diff:.15g} "
            f"reldiff={rel} tol={row.tolerance_pct:.2f}% status={row.status}"
        )


def main(argv: list[str]) -> int:
    repo_root = ensure_nyxus_repo(argv[1] if len(argv) > 1 else None)
    rows = run(repo_root)
    print_rows(rows)
    results_path, summary_path = write_results(rows, OUT_DIR)
    pass_count = sum(1 for row in rows if row.status == "PASS")
    mismatch_count = sum(1 for row in rows if row.status == "MISMATCH")
    print(
        f"PyRadiomics independent 2D rows: {pass_count} pass, "
        f"{mismatch_count} mismatch. Wrote {results_path} and {summary_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
