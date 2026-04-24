function verify_nyxus_first_order_new_tests(repoRoot)
% Cross-verify the full Nyxus 2D first-order intensity test suite.
%
% Function name is kept for backward compatibility even though the script
% now verifies the full current intensity surface covered in
% tests/test_pixel_intensity_features.h.
%
% This script:
% 1. Reads the existing phantom from Nyxus tests/test_data.h
% 2. Reads the asserted constants from tests/test_pixel_intensity_features.h
% 3. Recomputes the same features in MATLAB
% 4. Reports whether each value would pass Nyxus's agrees_gt() tolerance
% 5. Shows where raw MATLAB defaults differ from Nyxus-specific logic
% 6. Clones vjaganat90/nyxus branch unittest_add locally on demand when
%    repoRoot is not provided

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
testHeaderPath = fullfile(repoRoot, "tests", "test_pixel_intensity_features.h");

featureCases = build_feature_cases();
featureNames = cellfun(@(s) s.feature, featureCases, 'UniformOutput', false);

slideMin = 0.0;
slideMax = 65535.0;

intensities = load_pixel_intensity_phantom(testDataPath);
[assertions, missingFeatures] = load_expected_assertions(testHeaderPath, featureNames);

actualGrey24 = compute_nyxus_style_values(intensities, slideMin, slideMax, 24);
actualGrey20 = compute_nyxus_style_values(intensities, slideMin, slideMax, 20);

fprintf('\nNyxus first-order intensity verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Phantom source: %s\n', testDataPath);
fprintf('Test constants source: %s\n', testHeaderPath);
fprintf('Pixels in phantom: %d\n', numel(intensities));
fprintf('Covered-image slide range used by test: [%.1f, %.1f]\n', slideMin, slideMax);
fprintf('Known first-order features in verifier: %d\n', numel(featureCases));
fprintf('Features with direct assertions in this Nyxus checkout: %d\n\n', numel(fieldnames(assertions)));

fprintf('%-34s %18s %18s %12s %8s %8s\n', 'Feature', 'Expected', 'MATLAB', 'AbsDiff', 'Grey', 'Pass');
fprintf('%s\n', repmat('-', 1, 108));

checkedCount = 0;
passCount = 0;
failedFeatures = {};
for i = 1:numel(featureCases)
    featureCase = featureCases{i};
    featureName = featureCase.feature;
    if ~isfield(assertions, featureName)
        continue;
    end
    meta = assertions.(featureName);
    checkedCount = checkedCount + 1;

    if featureCase.greyBins == 20
        actual = actualGrey20;
    else
        actual = actualGrey24;
    end

    val = actual.(featureName);
    gt = meta.value;
    diffVal = abs(val - gt);
    pass = agrees_gt_matlab(val, gt, meta.fracTolerance);
    if pass
        passCount = passCount + 1;
    else
        failedFeatures{end + 1} = featureName; %#ok<AGROW>
    end

    fprintf('%-34s %18.12g %18.12g %12.6g %8d %8s\n', ...
        featureName, gt, val, diffVal, featureCase.greyBins, pass_fail(pass));
end

fprintf('\n');
if isempty(failedFeatures)
    fprintf('All %d checked first-order intensity features pass Nyxus''s C++ test tolerances.\n', checkedCount);
else
    fprintf('%d/%d checked first-order intensity features pass Nyxus''s C++ test tolerances.\n', passCount, checkedCount);
    fprintf('Mismatches in this checkout: %s\n', strjoin(failedFeatures, ', '));
end

if ~isempty(missingFeatures)
    fprintf('Verifier-supported features not directly asserted in this Nyxus checkout: %s\n', strjoin(missingFeatures, ', '));
end

rawPct = prctile(intensities, [1 10 25 75 90 99]);
nyxPct = [actualGrey24.P01 actualGrey24.P10 actualGrey24.P25 actualGrey24.P75 actualGrey24.P90 actualGrey24.P99];
labels = {'P01', 'P10', 'P25', 'P75', 'P90', 'P99'};

fprintf('\nRaw MATLAB prctile() sanity check\n');
fprintf('%-8s %18s %18s %18s\n', 'Pct', 'prctile()', 'Nyxus-style', 'AbsDiff');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:numel(labels)
    fprintf('%-8s %18.12g %18.12g %18.6g\n', ...
        labels{i}, rawPct(i), nyxPct(i), abs(rawPct(i) - nyxPct(i)));
end

fprintf('\nEntropy / uniformity grey-depth sanity check\n');
fprintf('%-14s %18s %18s\n', 'Metric', 'GREYDEPTH=20', 'GREYDEPTH=24');
fprintf('%s\n', repmat('-', 1, 54));
fprintf('%-14s %18.12g %18.12g\n', 'ENTROPY', actualGrey20.ENTROPY, actualGrey24.ENTROPY);
fprintf('%-14s %18.12g %18.12g\n', 'UNIFORMITY', actualGrey20.UNIFORMITY, actualGrey24.UNIFORMITY);

fprintf('\nNotes:\n');
fprintf('1. This verifies the direct first-order intensity assertions present in the selected Nyxus checkout over the shared inline phantom.\n');
fprintf('2. Percentiles, IQR, QCOD, robust mean, and robust MAD follow Nyxus''s histogram-based implementation, not raw sample percentiles.\n');
fprintf('3. ENTROPY and UNIFORMITY use GREYDEPTH=20 in the C++ tests, while most other histogram-dependent features use Nyxus''s default path.\n');
fprintf('4. This is a MATLAB cross-check of the current Nyxus formulas and asserted constants; it is not an external IBSI oracle.\n');
fprintf('5. Different Nyxus checkouts may expose different intensity assertions, so missing features are reported separately from mismatches.\n\n');

end

function repoRoot = resolve_default_nyxus_repo()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
cacheRoot = fullfile(projectRoot, ".cache");
repoRoot = fullfile(cacheRoot, "nyxus_vjaganat90_unittest_add");

if ~exist(cacheRoot, 'dir')
    mkdir(cacheRoot);
end

repoRoot = ensure_nyxus_clone(repoRoot);
end

function repoRoot = ensure_nyxus_clone(repoRoot)
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
testHeaderPath = fullfile(repoRoot, "tests", "test_pixel_intensity_features.h");

if exist(testDataPath, 'file') && exist(testHeaderPath, 'file')
    return;
end

if exist(repoRoot, 'dir')
    error(['Expected Nyxus checkout at "%s", but required test files were not found there. ' ...
        'Delete the directory or pass a valid repoRoot explicitly.'], repoRoot);
end

gitUrl = 'https://github.com/vjaganat90/nyxus.git';
gitBranch = 'unittest_add';
fprintf('Nyxus checkout not found. Cloning %s (branch %s) into %s\n', gitUrl, gitBranch, repoRoot);
cmd = sprintf('git clone --depth 1 --branch %s --single-branch %s "%s"', gitBranch, gitUrl, repoRoot);
[status, cmdout] = system(cmd);
if status ~= 0
    error('Failed to clone Nyxus.\nCommand: %s\nOutput:\n%s', cmd, cmdout);
end

if ~(exist(testDataPath, 'file') && exist(testHeaderPath, 'file'))
    error('Clone completed, but expected Nyxus test files are missing under %s', repoRoot);
end
end

function featureCases = build_feature_cases()
featureCases = {
    make_case('INTEGRATED_INTENSITY', 24)
    make_case('MIN', 24)
    make_case('MAX', 24)
    make_case('RANGE', 24)
    make_case('MEAN', 24)
    make_case('MEDIAN', 24)
    make_case('MODE', 24)
    make_case('STANDARD_DEVIATION', 24)
    make_case('SKEWNESS', 24)
    make_case('KURTOSIS', 24)
    make_case('EXCESS_KURTOSIS', 24)
    make_case('HYPERSKEWNESS', 24)
    make_case('HYPERFLATNESS', 24)
    make_case('MEAN_ABSOLUTE_DEVIATION', 24)
    make_case('ROBUST_MEAN_ABSOLUTE_DEVIATION', 24)
    make_case('STANDARD_ERROR', 24)
    make_case('ROOT_MEAN_SQUARED', 24)
    make_case('ENTROPY', 20)
    make_case('ENERGY', 24)
    make_case('UNIFORMITY', 20)
    make_case('UNIFORMITY_PIU', 24)
    make_case('P01', 24)
    make_case('P10', 24)
    make_case('P25', 24)
    make_case('P75', 24)
    make_case('P90', 24)
    make_case('P99', 24)
    make_case('INTERQUARTILE_RANGE', 24)
    make_case('COV', 24)
    make_case('COVERED_IMAGE_INTENSITY_RANGE', 24)
    make_case('MEDIAN_ABSOLUTE_DEVIATION', 24)
    make_case('QCOD', 24)
    make_case('ROBUST_MEAN', 24)
    make_case('STANDARD_DEVIATION_BIASED', 24)
    make_case('VARIANCE', 24)
    make_case('VARIANCE_BIASED', 24)
    };
end

function item = make_case(featureName, greyBins)
item = struct('feature', featureName, 'greyBins', greyBins);
end

function intensities = load_pixel_intensity_phantom(testDataPath)
txt = fileread(testDataPath);
block = regexp(txt, 'pixelIntensityFeaturesTestData\[\]\s*=\s*\{(.*?)\};', 'tokens', 'once');
if isempty(block)
    error('Could not find pixelIntensityFeaturesTestData in %s', testDataPath);
end

tokens = regexp(block{1}, '\{\s*\d+\s*,\s*\d+\s*,\s*(\d+)\s*\}', 'tokens');
if isempty(tokens)
    error('Could not parse any phantom intensities from %s', testDataPath);
end

intensities = zeros(numel(tokens), 1);
for i = 1:numel(tokens)
    intensities(i) = str2double(tokens{i}{1});
end
intensities = double(intensities(:));
end

function [assertions, missingFeatures] = load_expected_assertions(testHeaderPath, featureNames)
txt = fileread(testHeaderPath);
assertions = struct();
missingFeatures = {};

for i = 1:numel(featureNames)
    featureName = featureNames{i};
    pattern = sprintf('Feature2D::%s\\]\\[0\\],\\s*(.+?)(?:,\\s*([0-9.]+))?\\)\\);', featureName);
    token = regexp(txt, pattern, 'tokens', 'once');
    if isempty(token)
        missingFeatures{end + 1} = featureName; %#ok<AGROW>
        continue;
    end

    expr = strtrim(token{1});
    value = str2num(expr); %#ok<ST2NM>
    if isempty(value)
        error('Could not evaluate asserted constant "%s" for feature %s', expr, featureName);
    end

    fracTolerance = 1000.0;
    if numel(token) > 1 && ~isempty(token{2})
        fracTolerance = str2double(token{2});
    end

    assertions.(featureName) = struct( ...
        'value', double(value), ...
        'fracTolerance', double(fracTolerance), ...
        'expression', expr);
end
end

function values = compute_nyxus_style_values(x, slideMin, slideMax, nGreyBins)
x = double(x(:));
n = numel(x);
if n == 0
    error('Phantom is empty.');
end

minVal = min(x);
maxVal = max(x);
roiRange = maxVal - minVal;
meanVal = mean(x);
medianVal = median(x);
energyVal = sum(x .^ 2);
centered = x - meanVal;
sumSq = sum(centered .^ 2);
sum3 = sum(centered .^ 3);
sum4 = sum(centered .^ 4);
sum5 = sum(centered .^ 5);
sum6 = sum(centered .^ 6);

[pct, rmad, entropyVal, uniformityVal] = nyxus_histogram_stats(x, minVal, maxVal, nGreyBins);

values = struct();
values.INTEGRATED_INTENSITY = sum(x);
values.MIN = minVal;
values.MAX = maxVal;
values.RANGE = roiRange;
values.MEAN = meanVal;
values.MEDIAN = medianVal;
values.MODE = nyxus_mode(x);
values.STANDARD_DEVIATION = sqrt(sumSq / (n - 1));
values.STANDARD_DEVIATION_BIASED = sqrt(sumSq / n);
values.VARIANCE = sumSq / (n - 1);
values.VARIANCE_BIASED = sumSq / n;
values.COV = values.STANDARD_DEVIATION / meanVal;
values.STANDARD_ERROR = values.STANDARD_DEVIATION / sqrt(n);
values.SKEWNESS = sqrt(n) * sum3 / (sumSq ^ 1.5);
values.KURTOSIS = n * sum4 / (sumSq * sumSq);
values.EXCESS_KURTOSIS = values.KURTOSIS - 3.0;
values.HYPERSKEWNESS = sum5 / (n * (values.STANDARD_DEVIATION ^ 5));
values.HYPERFLATNESS = sum6 / (n * (values.STANDARD_DEVIATION ^ 6));
values.MEAN_ABSOLUTE_DEVIATION = mean(abs(x - meanVal));
values.MEDIAN_ABSOLUTE_DEVIATION = mean(abs(x - medianVal));
values.ENERGY = energyVal;
values.ROOT_MEAN_SQUARED = sqrt(energyVal / n);
values.ENTROPY = entropyVal;
values.UNIFORMITY = uniformityVal;
values.UNIFORMITY_PIU = (1.0 - roiRange / (maxVal + minVal)) * 100.0;
values.P01 = pct.P01;
values.P10 = pct.P10;
values.P25 = pct.P25;
values.P75 = pct.P75;
values.P90 = pct.P90;
values.P99 = pct.P99;
values.INTERQUARTILE_RANGE = pct.P75 - pct.P25;
values.QCOD = (pct.P75 - pct.P25) / (pct.P75 + pct.P25);

robustMask = x >= pct.P10 & x <= pct.P90;
values.ROBUST_MEAN = mean(x(robustMask));
values.ROBUST_MEAN_ABSOLUTE_DEVIATION = rmad;
values.COVERED_IMAGE_INTENSITY_RANGE = roiRange / (slideMax - slideMin);
end

function [pct, rmad, entropyVal, uniformityVal] = nyxus_histogram_stats(x, minVal, maxVal, nGreyBins)
roiRange = maxVal - minVal;
if roiRange == 0
    pct = struct('P01', minVal, 'P10', minVal, 'P25', minVal, ...
        'P75', minVal, 'P90', minVal, 'P99', minVal);
    rmad = 0.0;
    entropyVal = 0.0;
    uniformityVal = 1.0;
    return;
end

binW100 = roiRange / 100.0;
bins100 = zeros(101, 1);
for i = 1:numel(x)
    realIdx = (x(i) - minVal) / binW100;
    idx = floor(realIdx);
    if isnan(realIdx)
        idx = 0;
    end
    bins100(idx + 1) = bins100(idx + 1) + 1;
end
bins100(100) = bins100(100) + bins100(101);
bins100(101) = 0;

targets = struct( ...
    'P01', numel(x) * 0.01, ...
    'P10', numel(x) * 0.10, ...
    'P25', numel(x) * 0.25, ...
    'P75', numel(x) * 0.75, ...
    'P90', numel(x) * 0.90, ...
    'P99', numel(x) * 0.99);

percentileNames = fieldnames(targets);
percentileValues = struct('P01', 0.0, 'P10', 0.0, 'P25', 0.0, 'P75', 0.0, 'P90', 0.0, 'P99', 0.0);

runSum = 0;
for i = 0:99
    binCount = bins100(i + 1);
    if binCount == 0
        runSum = runSum + binCount;
        continue;
    end

    for j = 1:numel(percentileNames)
        name = percentileNames{j};
        target = targets.(name);
        if runSum <= target && target <= runSum + binCount
            percentileValues.(name) = (target - runSum) * binW100 / binCount + minVal + binW100 * i;
        end
    end

    runSum = runSum + binCount;
end

bandMask = x >= percentileValues.P10 & x <= percentileValues.P90;
bandValues = x(bandMask);
bandMean = mean(bandValues);
rmad = mean(abs(bandValues - bandMean));

binsCust = zeros(nGreyBins + 1, 1);
for i = 1:numel(x)
    gray = floor(((x(i) - minVal) / roiRange) * nGreyBins);
    binsCust(gray + 1) = binsCust(gray + 1) + 1;
end
binsCust(nGreyBins) = binsCust(nGreyBins) + binsCust(nGreyBins + 1);
binsCust(nGreyBins + 1) = 0;

entropyVal = 0.0;
uniformityVal = 0.0;
for i = 1:numel(binsCust)
    p = binsCust(i) / numel(x);
    entropyVal = entropyVal + p * log2(p + 2.2e-16);
    uniformityVal = uniformityVal + p * p;
end
entropyVal = -entropyVal;

pct = percentileValues;
end

function m = nyxus_mode(x)
uniqueVals = unique(x);
counts = zeros(numel(uniqueVals), 1);
for i = 1:numel(uniqueVals)
    counts(i) = sum(x == uniqueVals(i));
end
[~, idx] = max(counts);
m = uniqueVals(idx);
end

function tf = agrees_gt_matlab(fval, gt, fracTolerance)
tolerance = gt / fracTolerance;
tf = abs(fval - gt) <= abs(tolerance);
end

function label = pass_fail(tf)
if tf
    label = 'PASS';
else
    label = 'FAIL';
end
end
