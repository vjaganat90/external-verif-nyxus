function verify_nyxus_glcm_tests(repoRoot)
% Cross-verify Nyxus GLCM unit-test expectations in MATLAB.
%
% This script:
% 1. Clones the Nyxus branch under test if repoRoot is not provided
% 2. Reads the IBSI phantom slices from tests/test_data.h
% 3. Reads the asserted constants from tests/test_glcm.h
% 4. Recomputes all non-IBSI 2D GLCM features in MATLAB
% 5. Reports whether each value would pass Nyxus's agrees_gt() tolerance

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
glcmHeaderPath = fullfile(repoRoot, "tests", "test_glcm.h");

slices = load_ibsi_phantom_slices(testDataPath);
glcmTruth = load_named_constants(glcmHeaderPath);
glcmCases = build_glcm_cases();

glcmResults = cell(numel(slices), 1);
for i = 1:numel(slices)
    glcmResults{i} = compute_glcm_slice_features(slices{i}.matrix);
end

fprintf('\nNyxus GLCM verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Phantom source: %s\n', testDataPath);
fprintf('GLCM test constants source: %s\n', glcmHeaderPath);
fprintf('Slices in phantom: %d\n', numel(slices));
fprintf('GLCM configuration: MATLAB-style binning with 100 levels, offset 1, asymmetric matrix, angles {0,45,90,135}\n\n');

[checkedCount, passCount, missing] = report_glcm(glcmCases, glcmTruth, glcmResults);

fprintf('\nSummary\n');
fprintf('GLCM: %d/%d features pass Nyxus''s C++ test tolerances.\n', passCount, checkedCount);
if ~isempty(missing)
    fprintf('GLCM features not directly asserted in this checkout: %s\n', strjoin(missing, ', '));
end
fprintf('\nNotes:\n');
fprintf('1. GLCM _AVE features are cross-checked the same way Nyxus tests them: they use the same ground-truth values as the base metrics and average over angles per slice.\n');
fprintf('2. The MATLAB code intentionally reproduces Nyxus''s current non-IBSI semantics, including asymmetric GLCMs, 100-level MATLAB-style binning, raw-count entropy, and the current JVAR/HOM2 conventions.\n');
fprintf('3. This is an independent reimplementation of the matrix construction and feature formulas over the same phantom, not a call back into Nyxus.\n\n');

end

function repoRoot = resolve_default_nyxus_repo()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
cacheRoot = fullfile(projectRoot, ".cache");
repoRoot = fullfile(cacheRoot, "nyxus_vjaganat90_main");

if ~exist(cacheRoot, 'dir')
    mkdir(cacheRoot);
end

repoRoot = ensure_nyxus_clone(repoRoot);
end

function repoRoot = ensure_nyxus_clone(repoRoot)
requiredFiles = {
    fullfile(repoRoot, "tests", "test_data.h")
    fullfile(repoRoot, "tests", "test_glcm.h")
    };

if all(cellfun(@(p) exist(p, 'file') ~= 0, requiredFiles))
    return;
end

if exist(repoRoot, 'dir')
    error(['Expected Nyxus checkout at "%s", but required test files were not found there. ' ...
        'Delete the directory or pass a valid repoRoot explicitly.'], repoRoot);
end

gitUrl = 'https://github.com/vjaganat90/nyxus.git';
gitBranch = 'main';
fprintf('Nyxus checkout not found. Cloning %s (branch %s) into %s\n', gitUrl, gitBranch, repoRoot);
cmd = sprintf('git clone --depth 1 --branch %s --single-branch %s "%s"', gitBranch, gitUrl, repoRoot);
[status, cmdout] = system(cmd);
if status ~= 0
    error('Failed to clone Nyxus.\nCommand: %s\nOutput:\n%s', cmd, cmdout);
end
end

function slices = load_ibsi_phantom_slices(testDataPath)
intensityNames = {
    'ibsi_phantom_z1_intensity'
    'ibsi_phantom_z2_intensity'
    'ibsi_phantom_z3_intensity'
    'ibsi_phantom_z4_intensity'
    };
maskNames = {
    'ibsi_phantom_z1_mask'
    'ibsi_phantom_z2_mask'
    'ibsi_phantom_z3_mask'
    'ibsi_phantom_z4_mask'
    };

slices = cell(numel(intensityNames), 1);
for i = 1:numel(intensityNames)
    intenTriples = load_nyxus_pixel_array(testDataPath, intensityNames{i});
    maskTriples = load_nyxus_pixel_array(testDataPath, maskNames{i});
    slices{i} = struct( ...
        'name', intensityNames{i}, ...
        'matrix', build_masked_aux_image(intenTriples, maskTriples));
end
end

function triples = load_nyxus_pixel_array(headerPath, arrayName)
txt = fileread(headerPath);
pattern = sprintf('const\\s+static\\s+NyxusPixel\\s+%s\\[\\]\\s*=\\s*\\{([\\s\\S]*?)\\};', arrayName);
block = regexp(txt, pattern, 'tokens', 'once');
if isempty(block)
    error('Could not find array %s in %s', arrayName, headerPath);
end

tokens = regexp(block{1}, '\{\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*\}', 'tokens');
if isempty(tokens)
    error('Could not parse array %s in %s', arrayName, headerPath);
end

triples = zeros(numel(tokens), 3);
for i = 1:numel(tokens)
    triples(i, 1) = str2double(tokens{i}{1});
    triples(i, 2) = str2double(tokens{i}{2});
    triples(i, 3) = str2double(tokens{i}{3});
end
end

function auxImage = build_masked_aux_image(intensityTriples, maskTriples)
roiMask = maskTriples(:, 3) ~= 0;
roiPixels = intensityTriples(roiMask, :);
if isempty(roiPixels)
    error('ROI is empty.');
end

minX = min(roiPixels(:, 1));
maxX = max(roiPixels(:, 1));
minY = min(roiPixels(:, 2));
maxY = max(roiPixels(:, 2));

auxImage = zeros(maxY - minY + 1, maxX - minX + 1);
for i = 1:size(roiPixels, 1)
    x = roiPixels(i, 1) - minX + 1;
    y = roiPixels(i, 2) - minY + 1;
    auxImage(y, x) = roiPixels(i, 3);
end
end

function cases = build_glcm_cases()
baseFeatures = {
    'GLCM_ACOR'
    'GLCM_ASM'
    'GLCM_CLUPROM'
    'GLCM_CLUSHADE'
    'GLCM_CLUTEND'
    'GLCM_CONTRAST'
    'GLCM_CORRELATION'
    'GLCM_DIFAVE'
    'GLCM_DIFENTRO'
    'GLCM_DIFVAR'
    'GLCM_DIS'
    'GLCM_ENERGY'
    'GLCM_ENTROPY'
    'GLCM_HOM1'
    'GLCM_HOM2'
    'GLCM_ID'
    'GLCM_IDN'
    'GLCM_IDM'
    'GLCM_IDMN'
    'GLCM_INFOMEAS1'
    'GLCM_INFOMEAS2'
    'GLCM_IV'
    'GLCM_JAVE'
    'GLCM_JE'
    'GLCM_JMAX'
    'GLCM_JVAR'
    'GLCM_SUMAVERAGE'
    'GLCM_SUMENTROPY'
    'GLCM_SUMVARIANCE'
    'GLCM_VARIANCE'
    };

aveFeatures = {
    'GLCM_ASM_AVE'
    'GLCM_ACOR_AVE'
    'GLCM_CLUPROM_AVE'
    'GLCM_CLUSHADE_AVE'
    'GLCM_CLUTEND_AVE'
    'GLCM_CONTRAST_AVE'
    'GLCM_CORRELATION_AVE'
    'GLCM_DIFAVE_AVE'
    'GLCM_DIFENTRO_AVE'
    'GLCM_DIFVAR_AVE'
    'GLCM_DIS_AVE'
    'GLCM_ENERGY_AVE'
    'GLCM_ENTROPY_AVE'
    'GLCM_HOM1_AVE'
    'GLCM_ID_AVE'
    'GLCM_IDN_AVE'
    'GLCM_IDM_AVE'
    'GLCM_IDMN_AVE'
    'GLCM_IV_AVE'
    'GLCM_JAVE_AVE'
    'GLCM_JE_AVE'
    'GLCM_INFOMEAS1_AVE'
    'GLCM_INFOMEAS2_AVE'
    'GLCM_VARIANCE_AVE'
    'GLCM_JMAX_AVE'
    'GLCM_JVAR_AVE'
    'GLCM_SUMAVERAGE_AVE'
    'GLCM_SUMENTROPY_AVE'
    'GLCM_SUMVARIANCE_AVE'
    };

cases = cell(numel(baseFeatures) + numel(aveFeatures), 1);
for i = 1:numel(baseFeatures)
    cases{i} = struct('feature', baseFeatures{i}, 'truthKey', baseFeatures{i});
end
for i = 1:numel(aveFeatures)
    truthKey = extractBefore(aveFeatures{i}, "_AVE");
    cases{numel(baseFeatures) + i} = struct('feature', aveFeatures{i}, 'truthKey', char(truthKey));
end
end

function truth = load_named_constants(headerPath)
txt = fileread(headerPath);
tokens = regexp(txt, '\{\s*"([A-Z0-9_]+)"\s*,\s*([-+0-9.eE]+)\s*\}', 'tokens');
if isempty(tokens)
    error('Could not parse named constants in %s', headerPath);
end

truth = struct();
for i = 1:numel(tokens)
    name = tokens{i}{1};
    value = str2double(tokens{i}{2});
    if ~isfield(truth, name)
        truth.(name) = value;
    end
end
end

function [checkedCount, passCount, missing] = report_glcm(cases, truth, sliceResults)
checkedCount = 0;
passCount = 0;
fracTolerance = 100.0;
missing = {};

fprintf('GLCM\n');
fprintf('%-22s %18s %18s %12s %8s\n', 'Feature', 'Expected', 'MATLAB', 'AbsDiff', 'Pass');
fprintf('%s\n', repmat('-', 1, 86));

for i = 1:numel(cases)
    featureName = cases{i}.feature;
    truthKey = cases{i}.truthKey;
    if ~isfield(truth, truthKey)
        missing{end + 1} = featureName; %#ok<AGROW>
        continue;
    end

    expected = truth.(truthKey);
    actual = aggregate_feature(sliceResults, featureName);
    pass = agrees_gt_matlab(actual, expected, fracTolerance);

    checkedCount = checkedCount + 1;
    if pass
        passCount = passCount + 1;
    end

    fprintf('%-22s %18.12g %18.12g %12.6g %8s\n', ...
        featureName, expected, actual, abs(actual - expected), pass_fail(pass));
end
end

function actual = aggregate_feature(sliceResults, featureName)
vals = zeros(numel(sliceResults), 1);
for i = 1:numel(sliceResults)
    vals(i) = sliceResults{i}.(featureName);
end
actual = mean(vals);
end

function f = compute_glcm_slice_features(auxImage)
nLevels = 100;
angles = [1 0; 1 1; 0 1; -1 1];
baseNames = {
    'GLCM_ACOR'
    'GLCM_ASM'
    'GLCM_CLUPROM'
    'GLCM_CLUSHADE'
    'GLCM_CLUTEND'
    'GLCM_CONTRAST'
    'GLCM_CORRELATION'
    'GLCM_DIFAVE'
    'GLCM_DIFENTRO'
    'GLCM_DIFVAR'
    'GLCM_DIS'
    'GLCM_ENERGY'
    'GLCM_ENTROPY'
    'GLCM_HOM1'
    'GLCM_HOM2'
    'GLCM_ID'
    'GLCM_IDN'
    'GLCM_IDM'
    'GLCM_IDMN'
    'GLCM_INFOMEAS1'
    'GLCM_INFOMEAS2'
    'GLCM_IV'
    'GLCM_JAVE'
    'GLCM_JE'
    'GLCM_JMAX'
    'GLCM_JVAR'
    'GLCM_SUMAVERAGE'
    'GLCM_SUMENTROPY'
    'GLCM_SUMVARIANCE'
    'GLCM_VARIANCE'
    };

angled = struct();
for i = 1:numel(baseNames)
    angled.(baseNames{i}) = zeros(1, size(angles, 1));
end

binned = bin_image_matlab(auxImage, nLevels);
for ai = 1:size(angles, 1)
    P = build_glcm_matrix(binned, angles(ai, 1), angles(ai, 2), nLevels);
    vals = compute_glcm_angle_values(P, nLevels);
    for i = 1:numel(baseNames)
        angled.(baseNames{i})(ai) = vals.(baseNames{i});
    end
end

f = struct();
for i = 1:numel(baseNames)
    baseName = baseNames{i};
    f.(baseName) = mean(angled.(baseName));
end

aveNames = {
    'GLCM_ASM_AVE'
    'GLCM_ACOR_AVE'
    'GLCM_CLUPROM_AVE'
    'GLCM_CLUSHADE_AVE'
    'GLCM_CLUTEND_AVE'
    'GLCM_CONTRAST_AVE'
    'GLCM_CORRELATION_AVE'
    'GLCM_DIFAVE_AVE'
    'GLCM_DIFENTRO_AVE'
    'GLCM_DIFVAR_AVE'
    'GLCM_DIS_AVE'
    'GLCM_ENERGY_AVE'
    'GLCM_ENTROPY_AVE'
    'GLCM_HOM1_AVE'
    'GLCM_ID_AVE'
    'GLCM_IDN_AVE'
    'GLCM_IDM_AVE'
    'GLCM_IDMN_AVE'
    'GLCM_IV_AVE'
    'GLCM_JAVE_AVE'
    'GLCM_JE_AVE'
    'GLCM_INFOMEAS1_AVE'
    'GLCM_INFOMEAS2_AVE'
    'GLCM_VARIANCE_AVE'
    'GLCM_JMAX_AVE'
    'GLCM_JVAR_AVE'
    'GLCM_SUMAVERAGE_AVE'
    'GLCM_SUMENTROPY_AVE'
    'GLCM_SUMVARIANCE_AVE'
    };

for i = 1:numel(aveNames)
    baseName = char(extractBefore(aveNames{i}, "_AVE"));
    f.(aveNames{i}) = f.(baseName);
end
end

function D = bin_image_matlab(auxImage, nLevels)
maxI = max(auxImage(:));
if maxI <= 0
    error('Expected positive max intensity for GLCM verifier.');
end

slope = double(nLevels) / double(maxI);
intercept = 1.0;
D = zeros(size(auxImage));
for idx = 1:numel(auxImage)
    x = auxImage(idx);
    if x == 0
        y = 1;
    else
        y = floor(slope * double(x) + intercept);
        if y > nLevels
            y = nLevels;
        end
        if y < 1
            y = 1;
        end
    end
    D(idx) = y;
end
end

function P = build_glcm_matrix(D, dx, dy, nLevels)
[h, w] = size(D);
P = zeros(nLevels, nLevels);
for row = 1:h
    for col = 1:w
        rr = row + dy;
        cc = col + dx;
        if rr < 1 || rr > h || cc < 1 || cc > w
            continue;
        end
        a = D(rr, cc);
        b = D(row, col);
        P(a, b) = P(a, b) + 1;
    end
end
end

function vals = compute_glcm_angle_values(P, nLevels)
sumP = sum(P(:));
I = 1:nLevels;
EPS = 2.2e-16;

Pxpy = zeros(1, 2 * nLevels);
Pxmy = zeros(1, nLevels);
kValuesSum = zeros(1, 2 * nLevels);
kValuesDiff = zeros(1, nLevels);
for x = 1:nLevels
    for y = 1:nLevels
        pNorm = P(x, y) / sumP;
        Pxpy(x + y - 1) = Pxpy(x + y - 1) + pNorm;
        Pxmy(abs(x - y) + 1) = Pxmy(abs(x - y) + 1) + pNorm;
        kValuesSum(x + y - 1) = I(x) + I(y);
        kValuesDiff(abs(x - y) + 1) = abs(I(x) - I(y));
    end
end

px = zeros(1, nLevels);
for j = 1:nLevels
    for i = 1:nLevels
        px(i) = px(i) + P(i, j) / sumP;
    end
end
byRowMean = sum(px .* I);

vals = struct();
vals.GLCM_ASM = sum((P(:) / sumP) .^ 2);
vals.GLCM_ENERGY = vals.GLCM_ASM;
vals.GLCM_CONTRAST = sum(sum(P .* ((I' - I) .^ 2))) / sumP;

mr = 0.0;
mc = 0.0;
for c = 1:nLevels
    for r = 1:nLevels
        mr = mr + P(r, c) * I(r);
        mc = mc + P(r, c) * I(c);
    end
end
mr = mr / sumP;
mc = mc / sumP;

s2r = 0.0;
s2c = 0.0;
tmp1 = 0.0;
for c = 1:nLevels
    for r = 1:nLevels
        s2r = s2r + P(r, c) / sumP * (I(r) - mr) ^ 2;
        s2c = s2c + P(r, c) / sumP * (I(c) - mc) ^ 2;
        tmp1 = tmp1 + (I(r) - mr) * (I(c) - mc) * P(r, c) / sumP;
    end
end
vals.GLCM_CORRELATION = tmp1 / (sqrt(s2r) * sqrt(s2c));

meanVal = 0.0;
for r = 1:nLevels
    meanVal = meanVal + sum(P(r, :)) * I(r);
end
meanVal = meanVal / sumP;
vals.GLCM_VARIANCE = 0.0;
for r = 1:nLevels
    d = I(r) - meanVal;
    vals.GLCM_VARIANCE = vals.GLCM_VARIANCE + d * d * sum(P(r, :));
end
vals.GLCM_VARIANCE = vals.GLCM_VARIANCE / sumP;

vals.GLCM_IDM = sum(Pxmy ./ (1 + (0:nLevels - 1) .^ 2));
vals.GLCM_SUMAVERAGE = sum(kValuesSum .* Pxpy);
vals.GLCM_SUMENTROPY = -sum(Pxpy .* log2(Pxpy + EPS));
vals.GLCM_ENTROPY = -sum(P(:) .* log2(P(:) + EPS));
vals.GLCM_DIFAVE = sum(kValuesDiff .* Pxmy);

n = numel(Pxmy);
varBuf = zeros(1, n);
for x = 1:n
    for k = 1:n
        varBuf(k) = varBuf(k) + ((k - 1) - vals.GLCM_DIFAVE) ^ 2 * Pxmy(k);
    end
end
vals.GLCM_DIFVAR = sum(varBuf) / n;
vals.GLCM_DIFENTRO = -sum(Pxmy(Pxmy ~= 0) .* log2(Pxmy(Pxmy ~= 0) + EPS));

px = zeros(1, nLevels);
py = zeros(1, nLevels);
for i = 1:nLevels
    for j = 1:nLevels
        px(i) = px(i) + P(i, j) / sumP;
        py(j) = py(j) + P(i, j) / sumP;
    end
end
HXY = 0.0;
HXY1 = 0.0;
HXY2 = 0.0;
HX = 0.0;
for i = 1:nLevels
    HX = HX + px(i) * log2(px(i) + EPS);
    for j = 1:nLevels
        HXY = HXY + P(i, j) / sumP * log2(P(i, j) / sumP + EPS);
        HXY1 = HXY1 + P(i, j) / sumP * log2(px(i) * py(j) + EPS);
        HXY2 = HXY2 + px(i) * py(j) * log2(px(i) * py(j) + EPS);
    end
end
vals.GLCM_INFOMEAS1 = (HXY - HXY1) / HX;
vals.GLCM_INFOMEAS2 = sqrt(abs(1 - exp(-2 * (-HXY2 + HXY))));

vals.GLCM_ACOR = sum(sum(P .* (I' * I))) / sumP;
vals.GLCM_CLUPROM = 0.0;
vals.GLCM_CLUSHADE = 0.0;
vals.GLCM_CLUTEND = 0.0;
vals.GLCM_DIS = 0.0;
vals.GLCM_HOM2 = 0.0;
for r = 1:nLevels
    for c = 1:nLevels
        m = I(r) + I(c) - 2 * byRowMean;
        vals.GLCM_CLUPROM = vals.GLCM_CLUPROM + m ^ 4 * P(r, c) / sumP;
        vals.GLCM_CLUSHADE = vals.GLCM_CLUSHADE + m ^ 3 * P(r, c) / sumP;
        vals.GLCM_CLUTEND = vals.GLCM_CLUTEND + m ^ 2 * P(r, c) / sumP;
        vals.GLCM_DIS = vals.GLCM_DIS + abs(r - c) * P(r, c) / sumP;
        vals.GLCM_HOM2 = vals.GLCM_HOM2 + P(r, c) / (1 + abs(r - c) * abs(r - c));
    end
end
vals.GLCM_SUMVARIANCE = vals.GLCM_CLUTEND;
vals.GLCM_HOM1 = sum(sum(P / sumP ./ (1 + abs((1:nLevels)' - (1:nLevels)))));
vals.GLCM_ID = sum(Pxmy ./ (1 + (0:nLevels - 1)));
vals.GLCM_IDN = sum(Pxmy ./ (1 + (0:nLevels - 1) / nLevels));
vals.GLCM_IDMN = sum(Pxmy ./ (1 + ((0:nLevels - 1) .^ 2) / (nLevels * nLevels)));

vals.GLCM_IV = 0.0;
for k = 2:nLevels
    kval = kValuesDiff(k);
    if kval ~= 0
        vals.GLCM_IV = vals.GLCM_IV + Pxmy(k) / (kval * kval);
    end
end

vals.GLCM_JAVE = sum(sum(P .* I')) / sumP;
vals.GLCM_JE = -sum(sum((P / sumP) .* log2(P / sumP + EPS)));
vals.GLCM_JMAX = max(P(:) / sumP);
vals.GLCM_JVAR = 0.0;
for x = 1:nLevels
    for y = 1:nLevels
        vals.GLCM_JVAR = vals.GLCM_JVAR + ((x) - vals.GLCM_JAVE) ^ 2 * P(x, y) / sumP;
    end
end
end

function pass = agrees_gt_matlab(fval, groundTruth, fracTolerance)
diffVal = fval - groundTruth;
tolerance = groundTruth / fracTolerance;
pass = abs(diffVal) <= abs(tolerance);
end

function txt = pass_fail(pass)
if pass
    txt = 'PASS';
else
    txt = 'FAIL';
end
end
