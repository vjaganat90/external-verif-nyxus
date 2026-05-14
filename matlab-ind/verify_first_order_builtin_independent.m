function results = verify_first_order_builtin_independent(repoRoot, outDir)
% Compare Nyxus 2D first-order expectations to MATLAB built-in primitives.
%
% This script intentionally does not reimplement Nyxus histogram percentile
% logic. Histogram-specific Nyxus conventions are reported as mismatches or
% unsupported if MATLAB has no equivalent built-in definition.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = '/Users/jaganathv2/Code/Codex/nyxus';
end
if nargin < 2 || strlength(string(outDir)) == 0
    outDir = fileparts(mfilename('fullpath'));
end

dataPath = fullfile(repoRoot, 'tests', 'test_data.h');
headerPath = fullfile(repoRoot, 'tests', 'test_pixel_intensity_features.h');

x = load_pixel_intensity_phantom(dataPath);
truth = load_feature_assertions(headerPath);

rows = {};
rows = add_builtin(rows, truth, 'INTEGRATED_INTENSITY', sum(x), 'sum', 'Exact MATLAB sum of ROI intensities');
rows = add_builtin(rows, truth, 'MIN', min(x), 'min', 'MATLAB minimum');
rows = add_builtin(rows, truth, 'MAX', max(x), 'max', 'MATLAB maximum');
rows = add_builtin(rows, truth, 'RANGE', range(x), 'range', 'MATLAB range');
rows = add_builtin(rows, truth, 'MEAN', mean(x), 'mean', 'MATLAB arithmetic mean');
rows = add_builtin(rows, truth, 'MEDIAN', median(x), 'median', 'MATLAB median');
rows = add_builtin(rows, truth, 'MODE', mode(x), 'mode', 'MATLAB mode');
rows = add_builtin(rows, truth, 'STANDARD_DEVIATION', std(x, 0), 'std(x,0)', 'Unbiased sample standard deviation');
rows = add_builtin(rows, truth, 'STANDARD_DEVIATION_BIASED', std(x, 1), 'std(x,1)', 'Population standard deviation');
rows = add_builtin(rows, truth, 'VARIANCE', var(x, 0), 'var(x,0)', 'Unbiased sample variance');
rows = add_builtin(rows, truth, 'VARIANCE_BIASED', var(x, 1), 'var(x,1)', 'Population variance');
rows = add_builtin(rows, truth, 'STANDARD_ERROR', std(x, 0) / sqrt(numel(x)), 'std/sqrt(n)', 'Composition of MATLAB built-ins');
rows = add_builtin(rows, truth, 'ROOT_MEAN_SQUARED', sqrt(mean(x .^ 2)), 'sqrt(mean(x.^2))', 'RMS via MATLAB mean');
rows = add_builtin(rows, truth, 'SKEWNESS', skewness(x, 1), 'skewness(x,1)', 'MATLAB population-normalized skewness');
rows = add_builtin(rows, truth, 'KURTOSIS', kurtosis(x, 1), 'kurtosis(x,1)', 'MATLAB population-normalized kurtosis');
rows = add_builtin(rows, truth, 'EXCESS_KURTOSIS', kurtosis(x, 1) - 3.0, 'kurtosis(x,1)-3', 'MATLAB excess kurtosis');
rows = add_builtin(rows, truth, 'ENERGY', sum(x .^ 2), 'sum(x.^2)', 'Composition of MATLAB sum');
rows = add_builtin(rows, truth, 'COV', std(x, 0) / mean(x), 'std/mean', 'Composition of MATLAB std and mean');
rows = add_builtin(rows, truth, 'MEAN_ABSOLUTE_DEVIATION', mean(abs(x - mean(x))), 'mean(abs(x-mean(x)))', 'Composition of MATLAB mean');
rows = add_builtin(rows, truth, 'MEDIAN_ABSOLUTE_DEVIATION', mean(abs(x - median(x))), 'mean(abs(x-median(x)))', 'Nyxus uses mean absolute deviation around median, not MATLAB mad default');

% Built-in raw-sample percentiles are compared directly because Nyxus first
% order percentile outputs now use MATLAB-compatible exact sample interpolation.
pct = prctile(x, [1 10 25 75 90 99]);
rows = add_builtin(rows, truth, 'P01', pct(1), 'prctile', 'Raw sample percentile');
rows = add_builtin(rows, truth, 'P10', pct(2), 'prctile', 'Raw sample percentile');
rows = add_builtin(rows, truth, 'P25', pct(3), 'prctile', 'Raw sample percentile');
rows = add_builtin(rows, truth, 'P75', pct(4), 'prctile', 'Raw sample percentile');
rows = add_builtin(rows, truth, 'P90', pct(5), 'prctile', 'Raw sample percentile');
rows = add_builtin(rows, truth, 'P99', pct(6), 'prctile', 'Raw sample percentile');
rows = add_builtin(rows, truth, 'INTERQUARTILE_RANGE', iqr(x), 'iqr', 'Raw sample IQR');
rows = add_builtin(rows, truth, 'QCOD', iqr(x) / (pct(4) + pct(3)), 'iqr/prctile', 'Raw sample QCOD');

unsupported = {
    'ENTROPY', 'NO_BUILTIN', 'Nyxus histogram entropy convention is not MATLAB entropy() on an image';
    'UNIFORMITY', 'NO_BUILTIN', 'Nyxus histogram uniformity convention has no direct MATLAB first-order builtin';
    'UNIFORMITY_PIU', 'NO_BUILTIN', 'Nyxus PIU convention has no direct MATLAB builtin';
    'COVERED_IMAGE_INTENSITY_RANGE', 'NO_BUILTIN', 'Requires Nyxus slide min/max semantics';
    'ROBUST_MEAN', 'NO_BUILTIN', 'No single MATLAB built-in feature; Nyxus derives this from the P10-P90 raw-sample percentile band';
    'ROBUST_MEAN_ABSOLUTE_DEVIATION', 'NO_BUILTIN', 'No single MATLAB built-in feature; Nyxus derives this from the P10-P90 raw-sample percentile band';
    'HYPERSKEWNESS', 'NO_BUILTIN', 'No MATLAB built-in feature oracle';
    'HYPERFLATNESS', 'NO_BUILTIN', 'No MATLAB built-in feature oracle';
    };
for i = 1:size(unsupported, 1)
    rows = add_status(rows, truth, unsupported{i, 1}, NaN, unsupported{i, 2}, '', unsupported{i, 3});
end

results = rows_to_table(rows);
outPath = fullfile(outDir, 'first_order_builtin_independent_results.csv');
writetable(results, outPath);
print_result_rows(results, 'First-order');
fprintf('First-order MATLAB built-in independent rows: %d pass, %d mismatch, %d no-builtin. Wrote %s\n', ...
    sum(string(results.Status) == "PASS"), ...
    sum(string(results.Status) == "MISMATCH"), ...
    sum(string(results.Status) == "NO_BUILTIN"), outPath);
end

function rows = add_builtin(rows, truth, feature, actual, builtinName, note)
rows = add_status(rows, truth, feature, actual, '', builtinName, note);
end

function rows = add_status(rows, truth, feature, actual, forcedStatus, builtinName, note)
if ~isfield(truth, feature)
    return;
end
expected = truth.(feature).value;
tolerancePct = 5.0;
if strlength(string(forcedStatus)) > 0
    status = forcedStatus;
    absDiff = NaN;
    relDiffPct = NaN;
else
    absDiff = abs(actual - expected);
    relDiffPct = relative_diff_pct(actual, expected);
    if relDiffPct <= tolerancePct
        status = 'PASS';
    else
        status = 'MISMATCH';
    end
end
rows(end + 1, :) = {'First-order intensity', feature, expected, actual, absDiff, relDiffPct, tolerancePct, status, builtinName, note}; %#ok<AGROW>
end

function results = rows_to_table(rows)
results = cell2table(rows, 'VariableNames', { ...
    'Family', 'Feature', 'Expected', 'MatlabBuiltinValue', 'AbsDiff', 'RelDiffPct', 'TolerancePct', 'Status', 'BuiltinUsed', 'Note'});
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
    frac = 1000.0;
    if numel(tokens{i}) >= 3 && ~isempty(tokens{i}{3})
        frac = str2double(tokens{i}{3});
    end
    truth.(name) = struct('value', double(val), 'fracTolerance', double(frac));
end
end

function pct = relative_diff_pct(fval, gt)
if isnan(fval) || isnan(gt)
    pct = Inf;
elseif gt == 0
    if abs(fval) <= 1e-12
        pct = 0.0;
    else
        pct = Inf;
    end
else
    pct = 100.0 * abs(fval - gt) / abs(gt);
end
end

function print_result_rows(results, label)
for i = 1:height(results)
    fprintf('%s %-42s expected=%.15g matlab=%.15g absdiff=%.15g reldiff=%.6g%% tol=%.2f%% status=%s\n', ...
        label, results.Feature{i}, results.Expected(i), results.MatlabBuiltinValue(i), ...
        results.AbsDiff(i), results.RelDiffPct(i), results.TolerancePct(i), results.Status{i});
end
end
