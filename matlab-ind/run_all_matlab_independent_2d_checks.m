function summaryTable = run_all_matlab_independent_2d_checks(repoRoot, outDir)
% ===========================================================================
%  ENTRY POINT  -  run_all_matlab_independent_2d_checks(repoRoot, outDir)
% ===========================================================================
% Single MATLAB built-in oracle verifier for Nyxus 2D checks.
%
% The comparisons below are arranged by feature family. Each comparable row
% uses a MATLAB built-in or MATLAB Image Processing Toolbox built-in as the
% oracle. Rows with no direct built-in are kept out of the comparison CSVs
% and written to matlab_incomparable_no_oracle.csv.
%
% Usage:
%   run_all_matlab_independent_2d_checks                 % find local Nyxus
%   run_all_matlab_independent_2d_checks(localRepoPath)  % use a local checkout
%
% If repoRoot is omitted (or empty), ensure_nyxus_repo() looks for an existing
% checkout via NYXUS_REPO_DIR, ../nyxus, the current directory, or
% <external-verif-root>/.cache/nyxus. It does not update or reset any checkout.
% Set NYXUS_ALLOW_CLONE=1 only if you want a one-time clone into .cache/nyxus.
%
% ---------------------------------------------------------------------------
%  FUNCTION MAP (this is a single-file program; sections are marked with %%):
%    ENTRY POINT          run_all_matlab_independent_2d_checks (this function)
%    FAMILY RUNNERS       run_first_order_family, run_shape2d_family,
%                         run_glcm_family, run_remaining_feret_family,
%                         write_matlab_incomparable_list
%    COMPARISON HELPERS   add_glcm_comparison, add_comparison,
%                         add_result_row, rows_to_result_table,
%                         summarize_family
%    HEADER/FIXTURE LOAD  load_feature_assertions, load_pixel_intensity_phantom,
%                         load_named_map, load_shape_fixture,
%                         load_shape2d_mask_image, load_ibsi_phantom_slices,
%                         load_nyxus_pixel_array, glcm_truth_key
%    IMAGE BUILDERS       fixture_to_images, build_masked_aux_image,
%                         bin_image_nyxus_matlab_style, mean_over_slices
%    MATH & I/O UTILS     relative_diff_pct, print_result_rows
% ===========================================================================

if nargin < 1
    repoRoot = '';
end
if nargin < 2 || strlength(string(outDir)) == 0
    outDir = fileparts(mfilename('fullpath'));
end

repoRoot = char(ensure_nyxus_repo(char(repoRoot)));
outDir = char(outDir);

fprintf('\nNyxus 2D MATLAB built-in independent verification\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Output dir: %s\n\n', outDir);

firstOrder = run_first_order_family(repoRoot, outDir);
shape2d = run_shape2d_family(repoRoot, outDir);
glcm = run_glcm_family(repoRoot, outDir);
remainingFeret = run_remaining_feret_family(repoRoot, outDir);
write_matlab_incomparable_list(outDir);

summaryRows = [
    summarize_family('First-order intensity', firstOrder)
    summarize_family('Shape and morphology subset', shape2d)
    summarize_family('GLCM MATLAB-exposed subset', glcm)
    summarize_family('Remaining 2D Feret subset', remainingFeret)
    ];

summaryTable = cell2table(summaryRows, 'VariableNames', { ...
    'Family', 'ComparedRows', 'PassRows', 'MismatchRows', 'Notes'});
summaryPath = fullfile(outDir, 'matlab_independent_2d_summary.csv');
writetable(summaryTable, summaryPath);

disp(summaryTable);
fprintf('\nWrote summary: %s\n', summaryPath);
end

%% ========================= FAMILY RUNNERS =========================

function results = run_first_order_family(repoRoot, outDir)
% First-order intensity feature family.
dataPath = fullfile(repoRoot, 'tests', 'test_data.h');
headerPath = fullfile(repoRoot, 'tests', 'test_pixel_intensity_features.h');

x = load_pixel_intensity_phantom(dataPath);
truth = load_feature_assertions(headerPath);

rows = {};
rows = add_comparison(rows, 'First-order intensity', truth, 'INTEGRATED_INTENSITY', sum(x), 'sum', 'Exact MATLAB sum of ROI intensities');
rows = add_comparison(rows, 'First-order intensity', truth, 'MIN', min(x), 'min', 'MATLAB minimum');
rows = add_comparison(rows, 'First-order intensity', truth, 'MAX', max(x), 'max', 'MATLAB maximum');
rows = add_comparison(rows, 'First-order intensity', truth, 'RANGE', range(x), 'range', 'MATLAB range');
rows = add_comparison(rows, 'First-order intensity', truth, 'MEAN', mean(x), 'mean', 'MATLAB arithmetic mean');
rows = add_comparison(rows, 'First-order intensity', truth, 'MEDIAN', median(x), 'median', 'MATLAB median');
rows = add_comparison(rows, 'First-order intensity', truth, 'MODE', mode(x), 'mode', 'MATLAB mode');
rows = add_comparison(rows, 'First-order intensity', truth, 'STANDARD_DEVIATION', std(x, 0), 'std(x,0)', 'Unbiased sample standard deviation');
rows = add_comparison(rows, 'First-order intensity', truth, 'STANDARD_DEVIATION_BIASED', std(x, 1), 'std(x,1)', 'Population standard deviation');
rows = add_comparison(rows, 'First-order intensity', truth, 'VARIANCE', var(x, 0), 'var(x,0)', 'Unbiased sample variance');
rows = add_comparison(rows, 'First-order intensity', truth, 'VARIANCE_BIASED', var(x, 1), 'var(x,1)', 'Population variance');
rows = add_comparison(rows, 'First-order intensity', truth, 'STANDARD_ERROR', std(x, 0) / sqrt(numel(x)), 'std/sqrt(n)', 'Composition of MATLAB built-ins');
rows = add_comparison(rows, 'First-order intensity', truth, 'ROOT_MEAN_SQUARED', sqrt(mean(x .^ 2)), 'sqrt(mean(x.^2))', 'RMS via MATLAB mean');
rows = add_comparison(rows, 'First-order intensity', truth, 'SKEWNESS', skewness(x, 1), 'skewness(x,1)', 'MATLAB population-normalized skewness');
rows = add_comparison(rows, 'First-order intensity', truth, 'KURTOSIS', kurtosis(x, 1), 'kurtosis(x,1)', 'MATLAB population-normalized kurtosis');
rows = add_comparison(rows, 'First-order intensity', truth, 'EXCESS_KURTOSIS', kurtosis(x, 1) - 3.0, 'kurtosis(x,1)-3', 'MATLAB excess kurtosis');
rows = add_comparison(rows, 'First-order intensity', truth, 'ENERGY', sum(x .^ 2), 'sum(x.^2)', 'Composition of MATLAB sum');
rows = add_comparison(rows, 'First-order intensity', truth, 'COV', std(x, 0) / mean(x), 'std/mean', 'Composition of MATLAB std and mean');
rows = add_comparison(rows, 'First-order intensity', truth, 'MEAN_ABSOLUTE_DEVIATION', mean(abs(x - mean(x))), 'mean(abs(x-mean(x)))', 'Composition of MATLAB mean');
rows = add_comparison(rows, 'First-order intensity', truth, 'MEDIAN_ABSOLUTE_DEVIATION', mad(x, 1), 'mad(x,1)', 'MATLAB built-in median absolute deviation = median(|x-median|). Nyxus MEDIAN_ABSOLUTE_DEVIATION is actually the mean of |x-median|, so a mismatch here legitimately flags that naming/convention divergence rather than a recomputation of the Nyxus formula');

pct = prctile(x, [1 10 25 75 90 99]);
rows = add_comparison(rows, 'First-order intensity', truth, 'P01', pct(1), 'prctile', 'Raw sample percentile');
rows = add_comparison(rows, 'First-order intensity', truth, 'P10', pct(2), 'prctile', 'Raw sample percentile');
rows = add_comparison(rows, 'First-order intensity', truth, 'P25', pct(3), 'prctile', 'Raw sample percentile');
rows = add_comparison(rows, 'First-order intensity', truth, 'P75', pct(4), 'prctile', 'Raw sample percentile');
rows = add_comparison(rows, 'First-order intensity', truth, 'P90', pct(5), 'prctile', 'Raw sample percentile');
rows = add_comparison(rows, 'First-order intensity', truth, 'P99', pct(6), 'prctile', 'Raw sample percentile');
rows = add_comparison(rows, 'First-order intensity', truth, 'INTERQUARTILE_RANGE', iqr(x), 'iqr', 'Raw sample IQR');
rows = add_comparison(rows, 'First-order intensity', truth, 'QCOD', (pct(4) - pct(3)) / (pct(4) + pct(3)), 'prctile P25/P75', 'Standard quartile coefficient of dispersion (P75-P25)/(P75+P25), composed transparently from MATLAB prctile percentiles. MATLAB exposes no single QCOD built-in; this reuses the same percentile oracle already verified by the P25/P75 rows');
rows = add_comparison(rows, 'First-order intensity', truth, 'ROBUST_MEAN', trimmean(x, 20), 'trimmean(x,20)', 'MATLAB built-in 20 percent trimmed mean; comparable to Nyxus P10-P90 robust mean');

results = rows_to_result_table(rows);
outPath = fullfile(outDir, 'first_order_builtin_independent_results.csv');
writetable(results, outPath);
print_result_rows(results, 'First-order');
fprintf('First-order MATLAB built-in comparisons: %d pass, %d mismatch. Wrote %s\n', ...
    sum(string(results.Status) == "PASS"), ...
    sum(string(results.Status) == "MISMATCH"), outPath);
end

function results = run_shape2d_family(repoRoot, outDir)
% Shape and morphology feature family.
dataPath = fullfile(repoRoot, 'tests', 'test_data.h');
headerPath = fullfile(repoRoot, 'tests', 'test_shape_morphology_2d.h');

truth = load_named_map(headerPath, 'shape2d_truth');
roi = load_shape_fixture(dataPath);
[bw, inten, minX, minY] = fixture_to_images(roi);

stats = regionprops(bw, inten, ...
    'Area', 'Centroid', 'WeightedCentroid', 'BoundingBox', ...
    'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', ...
    'Orientation', 'Extent', 'ConvexArea', 'Solidity', 'EulerNumber', ...
    'Extrema');

if isempty(stats)
    error('No object found in shape fixture.');
end
s = stats(1);

rows = {};
rows = add_comparison(rows, 'Shape and morphology', truth, 'AREA_PIXELS_COUNT', s.Area, 'regionprops Area', 'Binary object area in pixels');
rows = add_comparison(rows, 'Shape and morphology', truth, 'AREA_UM2', s.Area * 4.0, 'regionprops Area * pixel area', 'Uses Nyxus test PIXELSIZEUM=2');
rows = add_comparison(rows, 'Shape and morphology', truth, 'CENTROID_X', s.Centroid(1) + minX - 1, 'regionprops Centroid', 'Converted from MATLAB 1-based column coordinate to Nyxus x');
rows = add_comparison(rows, 'Shape and morphology', truth, 'CENTROID_Y', s.Centroid(2) + minY - 1, 'regionprops Centroid', 'Converted from MATLAB 1-based row coordinate to Nyxus y');
rows = add_comparison(rows, 'Shape and morphology', truth, 'WEIGHTED_CENTROID_X', s.WeightedCentroid(1) + minX, 'regionprops WeightedCentroid', 'Nyxus weighted centroid is 1-based in current 2D tests');
rows = add_comparison(rows, 'Shape and morphology', truth, 'WEIGHTED_CENTROID_Y', s.WeightedCentroid(2) + minY, 'regionprops WeightedCentroid', 'Nyxus weighted centroid is 1-based in current 2D tests');
rows = add_comparison(rows, 'Shape and morphology', truth, 'BBOX_XMIN', floor(s.BoundingBox(1) + minX - 0.5), 'regionprops BoundingBox', 'Converted from MATLAB half-pixel bounding box');
rows = add_comparison(rows, 'Shape and morphology', truth, 'BBOX_YMIN', floor(s.BoundingBox(2) + minY - 0.5), 'regionprops BoundingBox', 'Converted from MATLAB half-pixel bounding box');
rows = add_comparison(rows, 'Shape and morphology', truth, 'BBOX_WIDTH', s.BoundingBox(3), 'regionprops BoundingBox', 'MATLAB bounding box width');
rows = add_comparison(rows, 'Shape and morphology', truth, 'BBOX_HEIGHT', s.BoundingBox(4), 'regionprops BoundingBox', 'MATLAB bounding box height');
rows = add_comparison(rows, 'Shape and morphology', truth, 'EXTENT', s.Extent, 'regionprops Extent', 'Area divided by bounding-box area');
rows = add_comparison(rows, 'Shape and morphology', truth, 'ASPECT_RATIO', s.BoundingBox(3) / s.BoundingBox(4), 'regionprops BoundingBox composition', 'Width divided by height');
rows = add_comparison(rows, 'Shape and morphology', truth, 'MAJOR_AXIS_LENGTH', s.MajorAxisLength, 'regionprops MajorAxisLength', 'Equivalent ellipse major axis');
rows = add_comparison(rows, 'Shape and morphology', truth, 'MINOR_AXIS_LENGTH', s.MinorAxisLength, 'regionprops MinorAxisLength', 'Equivalent ellipse minor axis');
rows = add_comparison(rows, 'Shape and morphology', truth, 'ECCENTRICITY', s.Eccentricity, 'regionprops Eccentricity', 'Equivalent ellipse eccentricity');
rows = add_comparison(rows, 'Shape and morphology', truth, 'ELONGATION', s.MinorAxisLength / s.MajorAxisLength, 'regionprops MajorAxisLength/MinorAxisLength composition', 'MATLAB R2025b regionprops does not expose an Elongation property on this fixture; direct axis measurements are built-ins');
rows = add_comparison(rows, 'Shape and morphology', truth, 'ORIENTATION', s.Orientation, 'regionprops Orientation', 'MATLAB uses opposite sign convention from current Nyxus tests');
rows = add_comparison(rows, 'Shape and morphology', truth, 'EULER_NUMBER', bweuler(bw, 8), 'bweuler', 'MATLAB binary Euler number; may differ by object/background convention');
rows = add_comparison(rows, 'Shape and morphology', truth, 'CONVEX_HULL_AREA', s.ConvexArea, 'regionprops ConvexArea', 'Direct MATLAB filled-pixel convex area');
rows = add_comparison(rows, 'Shape and morphology', truth, 'SOLIDITY', s.Solidity, 'regionprops Solidity', 'Direct MATLAB solidity');
for i = 1:8
    extremaX = s.Extrema(i, 1) + minX - 0.5;
    extremaY = s.Extrema(i, 2) + minY - 0.5;
    rows = add_comparison(rows, 'Shape and morphology', truth, sprintf('EXTREMA_P%d_X', i), extremaX, 'regionprops Extrema', 'MATLAB boundary extrema converted from half-pixel to Nyxus zero-based pixel coordinates');
    rows = add_comparison(rows, 'Shape and morphology', truth, sprintf('EXTREMA_P%d_Y', i), extremaY, 'regionprops Extrema', 'MATLAB boundary extrema converted from half-pixel to Nyxus zero-based pixel coordinates');
end

results = rows_to_result_table(rows);
outPath = fullfile(outDir, 'shape2d_builtin_independent_results.csv');
writetable(results, outPath);
print_result_rows(results, 'Shape2D');
fprintf('Shape2D MATLAB built-in comparisons: %d pass, %d mismatch. Wrote %s\n', ...
    sum(string(results.Status) == "PASS"), ...
    sum(string(results.Status) == "MISMATCH"), outPath);
end

function results = run_glcm_family(repoRoot, outDir)
% GLCM feature family exposed by MATLAB graycomatrix/graycoprops.
dataPath = fullfile(repoRoot, 'tests', 'test_data.h');
headerPath = fullfile(repoRoot, 'tests', 'test_glcm.h');

slices = load_ibsi_phantom_slices(dataPath);
truth = load_named_map(headerPath, 'glcm_values');

offsets = [0 1; 1 1; 1 0; 1 -1];
sliceVals = cell(numel(slices), 1);
for i = 1:numel(slices)
    image = slices{i};
    discretized = bin_image_nyxus_matlab_style(image, 100);
    glcms = graycomatrix(discretized, 'Offset', offsets, 'NumLevels', 100, ...
        'GrayLimits', [1 100], 'Symmetric', false);
    sliceVals{i} = graycoprops(glcms, {'Contrast', 'Correlation', 'Energy', 'Homogeneity'});
end

rows = {};
rows = add_glcm_comparison(rows, truth, 'GLCM_CONTRAST', mean_over_slices(sliceVals, 'Contrast'), 'Nyxus binning + graycomatrix + graycoprops Contrast', 'MATLAB-defined GLCM contrast after Nyxus test discretization');
rows = add_glcm_comparison(rows, truth, 'GLCM_CORRELATION', mean_over_slices(sliceVals, 'Correlation'), 'Nyxus binning + graycomatrix + graycoprops Correlation', 'MATLAB-defined GLCM correlation after Nyxus test discretization');
rows = add_glcm_comparison(rows, truth, 'GLCM_ASM', mean_over_slices(sliceVals, 'Energy'), 'Nyxus binning + graycoprops Energy', 'MATLAB Energy is angular second moment after Nyxus test discretization');
rows = add_glcm_comparison(rows, truth, 'GLCM_ENERGY', mean_over_slices(sliceVals, 'Energy'), 'Nyxus binning + graycoprops Energy', 'Nyxus names this output ENERGY/ASM in 2D tests');
rows = add_glcm_comparison(rows, truth, 'GLCM_HOM1', mean_over_slices(sliceVals, 'Homogeneity'), 'Nyxus binning + graycoprops Homogeneity', 'MATLAB homogeneity after Nyxus test discretization');
rows = add_glcm_comparison(rows, truth, 'GLCM_CONTRAST_AVE', mean_over_slices(sliceVals, 'Contrast'), 'Nyxus binning + graycomatrix + graycoprops Contrast', 'MATLAB-defined averaged GLCM contrast after Nyxus test discretization');
rows = add_glcm_comparison(rows, truth, 'GLCM_CORRELATION_AVE', mean_over_slices(sliceVals, 'Correlation'), 'Nyxus binning + graycomatrix + graycoprops Correlation', 'MATLAB-defined averaged GLCM correlation after Nyxus test discretization');
rows = add_glcm_comparison(rows, truth, 'GLCM_ASM_AVE', mean_over_slices(sliceVals, 'Energy'), 'Nyxus binning + graycoprops Energy', 'MATLAB averaged Energy is angular second moment after Nyxus test discretization');
rows = add_glcm_comparison(rows, truth, 'GLCM_ENERGY_AVE', mean_over_slices(sliceVals, 'Energy'), 'Nyxus binning + graycoprops Energy', 'Nyxus names this averaged output ENERGY/ASM in 2D tests');
rows = add_glcm_comparison(rows, truth, 'GLCM_HOM1_AVE', mean_over_slices(sliceVals, 'Homogeneity'), 'Nyxus binning + graycoprops Homogeneity', 'MATLAB averaged homogeneity after Nyxus test discretization');

results = rows_to_result_table(rows);
outPath = fullfile(outDir, 'glcm_builtin_independent_results.csv');
writetable(results, outPath);
print_result_rows(results, 'GLCM');
fprintf('GLCM MATLAB built-in comparisons: %d pass, %d mismatch. Wrote %s\n', ...
    sum(string(results.Status) == "PASS"), ...
    sum(string(results.Status) == "MISMATCH"), outPath);
end

function results = run_remaining_feret_family(repoRoot, outDir)
% Remaining 2D Feret feature family.
% Fixture and expected constants are parsed from the Nyxus repo (the same
% shape2d_morphology_mask Nyxus itself feeds the caliper features), so this
% family cannot drift from the C++ test.
dataPath = fullfile(repoRoot, 'tests', 'test_data.h');
headerPath = fullfile(repoRoot, 'tests', 'test_2d_remaining_features.h');

BW = load_shape2d_mask_image(dataPath);
truth = load_named_map(headerPath, 'remaining2d_truth');

feret = bwferet(BW, 'all');
if height(feret) ~= 1
    error('Expected one connected object from the 2D morphology fixture, found %d.', height(feret));
end

rows = {};
rows = add_comparison(rows, 'Remaining 2D Feret', truth, 'MIN_FERET_ANGLE', feret.MinAngle, 'bwferet MinAngle', 'Direct MATLAB Feret angle; angle convention differs from current Nyxus test');
rows = add_comparison(rows, 'Remaining 2D Feret', truth, 'MAX_FERET_ANGLE', feret.MaxAngle, 'bwferet MaxAngle', 'Direct MATLAB Feret angle; angle convention differs from current Nyxus test');
rows = add_comparison(rows, 'Remaining 2D Feret', truth, 'STAT_FERET_DIAM_MIN', feret.MinDiameter, 'bwferet MinDiameter', 'Closest direct MATLAB Feret diameter primitive');
rows = add_comparison(rows, 'Remaining 2D Feret', truth, 'STAT_FERET_DIAM_MAX', feret.MaxDiameter, 'bwferet MaxDiameter', 'Closest direct MATLAB Feret diameter primitive');

results = rows_to_result_table(rows);
outPath = fullfile(outDir, 'remaining2d_feret_builtin_independent_results.csv');
writetable(results, outPath);
print_result_rows(results, 'Remaining2D Feret');
fprintf('Remaining2D Feret MATLAB built-in rows: %d pass, %d mismatch. Wrote %s\n', ...
    sum(string(results.Status) == "PASS"), sum(string(results.Status) == "MISMATCH"), outPath);
end

function write_matlab_incomparable_list(outDir)
% Plain reviewer-facing list of MATLAB rows/families with no direct oracle.
rows = {
    'MATLAB built-ins', 'First-order intensity', 'ENTROPY', 'Nyxus histogram entropy convention is not MATLAB entropy() on an image';
    'MATLAB built-ins', 'First-order intensity', 'UNIFORMITY', 'Nyxus histogram uniformity convention has no direct MATLAB first-order builtin';
    'MATLAB built-ins', 'First-order intensity', 'UNIFORMITY_PIU', 'Nyxus PIU convention has no direct MATLAB builtin';
    'MATLAB built-ins', 'First-order intensity', 'COVERED_IMAGE_INTENSITY_RANGE', 'Requires Nyxus slide min/max semantics';
    'MATLAB built-ins', 'First-order intensity', 'ROBUST_MEAN_ABSOLUTE_DEVIATION', 'No single MATLAB built-in feature; Nyxus derives this from the P10-P90 raw-sample percentile band';
    'MATLAB built-ins', 'First-order intensity', 'HYPERSKEWNESS', 'No MATLAB built-in feature oracle';
    'MATLAB built-ins', 'First-order intensity', 'HYPERFLATNESS', 'No MATLAB built-in feature oracle';
    'MATLAB built-ins', 'Shape and morphology', 'MASS_DISPLACEMENT', 'No single MATLAB morphology builtin; derived from Nyxus centroid conventions';
    'MATLAB built-ins', 'Shape and morphology', 'COMPACTNESS', 'Nyxus-specific compactness formula';
    'MATLAB built-ins', 'Shape and morphology', 'DIAMETER_EQUAL_AREA', 'Nyxus formula differs from MATLAB EquivDiameter convention';
    'MATLAB built-ins', 'Shape and morphology', 'ROUNDNESS', 'Nyxus formula depends on its major-axis convention';
    'MATLAB built-ins', 'Shape and morphology', 'PERIMETER', 'Nyxus ordered contour convention differs from regionprops Perimeter';
    'MATLAB built-ins', 'Shape and morphology', 'DIAMETER_EQUAL_PERIMETER', 'Depends on Nyxus perimeter convention';
    'MATLAB built-ins', 'Shape and morphology', 'EDGE_MEAN_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'MATLAB built-ins', 'Shape and morphology', 'EDGE_STDDEV_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'MATLAB built-ins', 'Shape and morphology', 'EDGE_MAX_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'MATLAB built-ins', 'Shape and morphology', 'EDGE_MIN_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'MATLAB built-ins', 'Shape and morphology', 'EDGE_INTEGRATED_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'MATLAB built-ins', 'Shape and morphology', 'CIRCULARITY', 'Depends on Nyxus perimeter convention';
    'MATLAB built-ins', 'Shape and morphology', 'FRACT_DIM_BOXCOUNT', 'No direct MATLAB built-in oracle';
    'MATLAB built-ins', 'Shape and morphology', 'FRACT_DIM_PERIMETER', 'No direct MATLAB built-in oracle';
    'MATLAB built-ins', 'Shape and morphology', 'DIAMETER_MIN_ENCLOSING_CIRCLE', 'No direct MATLAB built-in oracle in base/IPT';
    'MATLAB built-ins', 'Shape and morphology', 'DIAMETER_CIRCUMSCRIBING_CIRCLE', 'No direct MATLAB built-in oracle in base/IPT';
    'MATLAB built-ins', 'Shape and morphology', 'DIAMETER_INSCRIBING_CIRCLE', 'No direct MATLAB built-in oracle in base/IPT';
    'MATLAB built-ins', 'Shape and morphology', 'ROI_RADIUS_MEAN', 'Nyxus-specific contour distance convention';
    'MATLAB built-ins', 'Shape and morphology', 'ROI_RADIUS_MAX', 'Nyxus-specific contour distance convention';
    'MATLAB built-ins', 'Shape and morphology', 'ROI_RADIUS_MEDIAN', 'Nyxus-specific contour distance convention';
    'MATLAB built-ins', 'GLRLM', 'GLRLM_*', 'MATLAB Image Processing Toolbox has no full GLRLM matrix/feature family oracle';
    'MATLAB built-ins', 'GLDZM', 'GLDZM_*', 'MATLAB has no GLDZM feature-family oracle';
    'MATLAB built-ins', 'GLSZM', 'GLSZM_*', 'MATLAB can label connected zones, but has no GLSZM feature-family oracle';
    'MATLAB built-ins', 'GLDM', 'GLDM_*', 'MATLAB has no GLDM feature-family oracle';
    'MATLAB built-ins', 'NGLDM', 'NGLDM_*', 'MATLAB has no NGLDM feature-family oracle';
    'MATLAB built-ins', 'NGTDM', 'NGTDM_*', 'MATLAB has no NGTDM feature-family oracle';
    'MATLAB built-ins', 'Gabor score', 'GABOR', 'MATLAB has gabor/imgaborfilt but not Nyxus custom CPU score and threshold statistic';
    'MATLAB built-ins', 'Neighbor / touching', 'Neighbor / touching family', 'MATLAB has geometry primitives but no Nyxus multi-ROI neighbor/touching feature oracle';
    'MATLAB built-ins', 'Advanced shape long tail', 'Advanced shape long tail', 'Some primitives exist, but Feret/Martin/Nassenstein/chords/polygonality/geodesic/thickness need separate oracles';
    'MATLAB built-ins', 'Moments / Zernike / radial', 'Moments / Zernike / radial', 'Some moment primitives exist in MATLAB/File Exchange ecosystems, but no single trusted built-in oracle is wired here';
    };
glcmUnsupported = {
    'GLCM_ACOR'; 'GLCM_CLUPROM'; 'GLCM_CLUSHADE'; 'GLCM_CLUTEND'; 'GLCM_DIFAVE';
    'GLCM_DIFENTRO'; 'GLCM_DIFVAR'; 'GLCM_DIS'; 'GLCM_ENTROPY'; 'GLCM_HOM2';
    'GLCM_ID'; 'GLCM_IDN'; 'GLCM_IDM'; 'GLCM_IDMN'; 'GLCM_INFOMEAS1';
    'GLCM_INFOMEAS2'; 'GLCM_IV'; 'GLCM_JAVE'; 'GLCM_JE'; 'GLCM_JMAX';
    'GLCM_JVAR'; 'GLCM_SUMAVERAGE'; 'GLCM_SUMENTROPY'; 'GLCM_SUMVARIANCE'; 'GLCM_VARIANCE'};
for i = 1:numel(glcmUnsupported)
    feature = glcmUnsupported{i};
    rows(end + 1, :) = {'MATLAB built-ins', 'GLCM', feature, 'Not exposed by graycoprops'}; %#ok<AGROW>
    if ~strcmp(feature, 'GLCM_HOM2')
        rows(end + 1, :) = {'MATLAB built-ins', 'GLCM', [feature '_AVE'], 'Not exposed by graycoprops for the averaged Nyxus GLCM feature'}; %#ok<AGROW>
    end
end
results = cell2table(rows, 'VariableNames', {'Oracle', 'Family', 'Feature', 'Reason'});
outPath = fullfile(outDir, 'matlab_incomparable_no_oracle.csv');
writetable(results, outPath);
fprintf('MATLAB incomparable/no-oracle list wrote %s\n', outPath);
end

%% ========================= COMPARISON HELPERS =========================

function rows = add_glcm_comparison(rows, truth, feature, actual, builtinName, note)
truthFeature = glcm_truth_key(feature);
if ~isfield(truth, truthFeature)
    return;
end
rows = add_result_row(rows, 'GLCM', truth, truthFeature, actual, builtinName, note, feature);
end

function rows = add_comparison(rows, family, truth, feature, actual, builtinName, note)
rows = add_result_row(rows, family, truth, feature, actual, builtinName, note, feature);
end

function rows = add_result_row(rows, family, truth, feature, actual, builtinName, note, outputFeature)
if ~isfield(truth, feature)
    return;
end

expected = double(truth.(feature));
tolerancePct = 5.0;
absDiff = abs(actual - expected);
relDiffPct = relative_diff_pct(actual, expected);
if relDiffPct <= tolerancePct
    status = 'PASS';
else
    status = 'MISMATCH';
end
rows(end + 1, :) = {family, outputFeature, expected, actual, absDiff, relDiffPct, tolerancePct, status, builtinName, note};
end

function results = rows_to_result_table(rows)
results = cell2table(rows, 'VariableNames', { ...
    'Family', 'Feature', 'Expected', 'MatlabBuiltinValue', 'AbsDiff', ...
    'RelDiffPct', 'TolerancePct', 'Status', 'BuiltinUsed', 'Note'});
end

function row = summarize_family(name, t)
statuses = string(t.Status);
row = { ...
    name, ...
    height(t), ...
    sum(statuses == "PASS"), ...
    sum(statuses == "MISMATCH"), ...
    'MATLAB built-in oracle comparison'};
end

%% ===================== HEADER / FIXTURE PARSERS =====================

function truth = load_feature_assertions(testHeaderPath)
txt = fileread(testHeaderPath);
tokens = regexp(txt, 'Feature2D::([A-Z0-9_]+)\]\[0\],\s*([^,\)]+)(?:,\s*([0-9.]+))?\)\);', 'tokens');
truth = struct();
for i = 1:numel(tokens)
    name = tokens{i}{1};
    val = str2double(tokens{i}{2});
    if isnan(val)
        val = str2num(tokens{i}{2}); %#ok<ST2NM>
    end
    truth.(name) = double(val);
end
end

function x = load_pixel_intensity_phantom(testDataPath)
txt = fileread(testDataPath);
block = regexp(txt, 'pixelIntensityFeaturesTestData\[\]\s*=\s*\{(.*?)\};', 'tokens', 'once');
if isempty(block)
    error('Could not find pixelIntensityFeaturesTestData in %s', testDataPath);
end
tokens = regexp(block{1}, '\{\s*\d+\s*,\s*\d+\s*,\s*([0-9.+\-eE]+)\s*\}', 'tokens');
x = zeros(numel(tokens), 1);
for i = 1:numel(tokens)
    x(i) = str2double(tokens{i}{1});
end
end

function truth = load_named_map(headerPath, mapName)
% Parse a C++  std::unordered_map<std::string,double> NAME { {"K", v}, ... };
% into a MATLAB struct keyed by feature name. Used for shape2d_truth,
% glcm_values, and remaining2d_truth.
txt = fileread(headerPath);
block = regexp(txt, [mapName '\s*\{([\s\S]*?)\};'], 'tokens', 'once');
if isempty(block)
    error('Could not find %s in %s', mapName, headerPath);
end
tokens = regexp(block{1}, '\{"([^"]+)",\s*([^}]+)\}', 'tokens');
truth = struct();
for i = 1:numel(tokens)
    truth.(tokens{i}{1}) = str2double(tokens{i}{2});
end
end

function bw = load_shape2d_mask_image(dataPath)
% Build the full-extent binary ROI image from shape2d_morphology_mask. This is
% the same mask Nyxus feeds the 2D caliper/Feret features.
txt = fileread(dataPath);
mask = load_nyxus_pixel_array(txt, 'shape2d_morphology_mask');
width = max(mask(:, 1)) + 1;
height = max(mask(:, 2)) + 1;
bw = false(height, width);
for i = 1:size(mask, 1)
    if mask(i, 3) ~= 0
        bw(mask(i, 2) + 1, mask(i, 1) + 1) = true;
    end
end
end

function truthKey = glcm_truth_key(feature)
suffix = '_AVE';
if numel(feature) > numel(suffix) && strcmp(feature(end - numel(suffix) + 1:end), suffix)
    truthKey = feature(1:end - numel(suffix));
else
    truthKey = feature;
end
end

function roi = load_shape_fixture(dataPath)
txt = fileread(dataPath);
inten = load_nyxus_pixel_array(txt, 'shape2d_morphology_intensity');
mask = load_nyxus_pixel_array(txt, 'shape2d_morphology_mask');
roi = [];
for i = 1:size(mask, 1)
    if mask(i, 3) ~= 0
        idx = find(inten(:, 1) == mask(i, 1) & inten(:, 2) == mask(i, 2), 1);
        roi(end + 1, :) = [mask(i, 1), mask(i, 2), inten(idx, 3)]; %#ok<AGROW>
    end
end
end

function slices = load_ibsi_phantom_slices(dataPath)
txt = fileread(dataPath);
slices = {};
for z = 1:4
    inten = load_nyxus_pixel_array(txt, sprintf('ibsi_phantom_z%d_intensity', z));
    mask = load_nyxus_pixel_array(txt, sprintf('ibsi_phantom_z%d_mask', z));
    slices{end + 1} = build_masked_aux_image(inten, mask); %#ok<AGROW>
end
end

function triples = load_nyxus_pixel_array(txt, arrayName)
pattern = sprintf('const\\s+static\\s+NyxusPixel\\s+%s\\[\\]\\s*=\\s*\\{([\\s\\S]*?)\\};', arrayName);
block = regexp(txt, pattern, 'tokens', 'once');
if isempty(block)
    error('Could not find %s', arrayName);
end
tokens = regexp(block{1}, '\{\s*([0-9.+\-eE]+)\s*,\s*([0-9.+\-eE]+)\s*,\s*([0-9.+\-eE]+)\s*\}', 'tokens');
triples = zeros(numel(tokens), 3);
for i = 1:numel(tokens)
    triples(i, :) = [str2double(tokens{i}{1}), str2double(tokens{i}{2}), str2double(tokens{i}{3})];
end
end

%% ========================= IMAGE BUILDERS =========================

function [bw, inten, minX, minY] = fixture_to_images(roi)
x = roi(:, 1);
y = roi(:, 2);
minX = min(x);
minY = min(y);
maxX = max(x);
maxY = max(y);
bw = false(maxY - minY + 1, maxX - minX + 1);
inten = zeros(size(bw));
for i = 1:size(roi, 1)
    row = roi(i, 2) - minY + 1;
    col = roi(i, 1) - minX + 1;
    bw(row, col) = true;
    inten(row, col) = roi(i, 3);
end
end

function auxImage = build_masked_aux_image(intensityTriples, maskTriples)
x = maskTriples(:, 1);
y = maskTriples(:, 2);
minX = min(x);
minY = min(y);
maxX = max(x);
maxY = max(y);
auxImage = zeros(maxY - minY + 1, maxX - minX + 1);
for i = 1:size(maskTriples, 1)
    if maskTriples(i, 3) == 0
        continue;
    end
    idx = find(intensityTriples(:, 1) == maskTriples(i, 1) & intensityTriples(:, 2) == maskTriples(i, 2), 1);
    row = maskTriples(i, 2) - minY + 1;
    col = maskTriples(i, 1) - minX + 1;
    auxImage(row, col) = intensityTriples(idx, 3);
end
end

function discretized = bin_image_nyxus_matlab_style(auxImage, nLevels)
maxI = max(auxImage(:));
if maxI <= 0
    error('Expected positive max intensity for GLCM verifier.');
end

discretized = zeros(size(auxImage));
slope = double(nLevels) / double(maxI);
for idx = 1:numel(auxImage)
    x = auxImage(idx);
    if x == 0
        discretized(idx) = 1;
    else
        y = floor(slope * double(x) + 1.0);
        discretized(idx) = min(max(y, 1), nLevels);
    end
end
end

function val = mean_over_slices(sliceVals, propName)
vals = [];
for i = 1:numel(sliceVals)
    vals = [vals, sliceVals{i}.(propName)]; %#ok<AGROW>
end
val = mean(vals(:));
end

%% ===================== MATH & I/O UTILITIES =====================

function pct = relative_diff_pct(actual, expected)
if isnan(actual) || isnan(expected)
    pct = Inf;
elseif expected == 0
    if abs(actual) <= 1e-12
        pct = 0.0;
    else
        pct = Inf;
    end
else
    pct = 100.0 * abs(actual - expected) / abs(expected);
end
end

function print_result_rows(results, label)
for i = 1:height(results)
    fprintf('%s %-42s expected=%.15g matlab=%.15g absdiff=%.15g reldiff=%.6g%% tol=%.2f%% status=%s\n', ...
        label, results.Feature{i}, results.Expected(i), results.MatlabBuiltinValue(i), ...
        results.AbsDiff(i), results.RelDiffPct(i), results.TolerancePct(i), results.Status{i});
end
end
