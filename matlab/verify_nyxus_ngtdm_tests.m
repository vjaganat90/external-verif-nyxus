function verify_nyxus_ngtdm_tests(repoRoot)
% Cross-verify Nyxus 2D NGTDM unit-test expectations in MATLAB.
%
% This script:
% 1. Clones the Nyxus repo under test if repoRoot is not provided
% 2. Reads the IBSI digital phantom slices from tests/test_data.h
% 3. Reads active non-IBSI and IBSI NGTDM constants from C++ test headers
% 4. Recomputes the 2D NGTDM features in MATLAB
% 5. Reports whether each value would pass Nyxus's C++ test tolerances

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
nonIbsiHeaderPath = fullfile(repoRoot, "tests", "test_ngtdm.h");
ibsiHeaderPath = fullfile(repoRoot, "tests", "test_ibsi_ngtdm.h");

slices = load_ibsi_phantom_slices(testDataPath);
nonIbsiTruth = load_named_constants_from_map(nonIbsiHeaderPath, 'ngtdm_values');
ibsiTruth = load_named_constants_from_map(ibsiHeaderPath, 'IBSI_ngtdm_values');
features = build_ngtdm_features();

nonIbsiResults = cell(numel(slices), 1);
ibsiResults = cell(numel(slices), 1);
for i = 1:numel(slices)
    nonIbsiResults{i} = compute_ngtdm_slice_features(slices{i}.matrix, false, 100);
    ibsiResults{i} = compute_ngtdm_slice_features(slices{i}.matrix, true, 0);
end

fprintf('\nNyxus 2D NGTDM verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Phantom source: %s\n', testDataPath);
fprintf('Non-IBSI constants source: %s\n', nonIbsiHeaderPath);
fprintf('IBSI constants source: %s\n', ibsiHeaderPath);
fprintf('Slices in phantom: %d\n', numel(slices));
fprintf('Configuration: active 2D NGTDM C++ regression surfaces over the IBSI digital phantom slices\n\n');

[checkedNonIbsi, passedNonIbsi] = report_ngtdm('2D NGTDM non-IBSI', features, nonIbsiTruth, nonIbsiResults);
fprintf('\n');
[checkedIbsi, passedIbsi] = report_ngtdm('2D NGTDM IBSI', features, ibsiTruth, ibsiResults);

fprintf('\nSummary\n');
fprintf('NGTDM non-IBSI: %d/%d features pass Nyxus''s C++ test tolerances.\n', passedNonIbsi, checkedNonIbsi);
fprintf('NGTDM IBSI:     %d/%d features pass Nyxus''s C++ test tolerances.\n', passedIbsi, checkedIbsi);
fprintf('NGTDM total:    %d/%d checks pass.\n', passedNonIbsi + passedIbsi, checkedNonIbsi + checkedIbsi);
fprintf('\nNotes:\n');
fprintf('1. Non-IBSI mode mirrors Nyxus''s current 100-level MATLAB-style binning, including zero-valued AABB background mapping to gray level 1.\n');
fprintf('2. IBSI mode uses unbinned phantom gray levels and excludes zero-valued background from neighborhood accumulation.\n');
fprintf('3. This is an independent MATLAB implementation over the same shared fixture, not a call back into Nyxus.\n\n');

if passedNonIbsi + passedIbsi ~= checkedNonIbsi + checkedIbsi
    error('NGTDM cross-verification failed.');
end
end

function repoRoot = resolve_default_nyxus_repo()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
cacheRoot = fullfile(projectRoot, ".cache");
repoRoot = fullfile(cacheRoot, "nyxus_vjaganat90_main");
gitUrl = 'https://github.com/vjaganat90/nyxus.git';
gitBranch = 'main';

if ~exist(cacheRoot, 'dir')
    mkdir(cacheRoot);
end

repoRoot = ensure_nyxus_clone(repoRoot, gitUrl, gitBranch);
end

function repoRoot = ensure_nyxus_clone(repoRoot, gitUrl, gitBranch)
requiredFiles = {
    fullfile(repoRoot, "tests", "test_data.h")
    fullfile(repoRoot, "tests", "test_ngtdm.h")
    fullfile(repoRoot, "tests", "test_ibsi_ngtdm.h")
    };

if all(cellfun(@(p) exist(p, 'file') ~= 0, requiredFiles))
    refresh_existing_clone(repoRoot, gitBranch);
    return;
end

if exist(repoRoot, 'dir')
    error(['Expected Nyxus checkout at "%s", but required test files were not found there. ' ...
        'Delete the directory or pass a valid repoRoot explicitly.'], repoRoot);
end

fprintf('Nyxus checkout not found. Cloning %s (branch %s) into %s\n', gitUrl, gitBranch, repoRoot);
cmd = sprintf('git clone --depth 1 --branch %s --single-branch %s "%s"', gitBranch, gitUrl, repoRoot);
[status, cmdout] = system(cmd);
if status ~= 0
    error('Failed to clone Nyxus.\nCommand: %s\nOutput:\n%s', cmd, cmdout);
end
end

function refresh_existing_clone(repoRoot, gitBranch)
fprintf('Refreshing cached Nyxus checkout at %s to origin/%s\n', repoRoot, gitBranch);
cmd = sprintf('git -C "%s" fetch --depth 1 origin %s && git -C "%s" checkout -q %s && git -C "%s" reset --hard -q origin/%s', ...
    repoRoot, gitBranch, repoRoot, gitBranch, repoRoot, gitBranch);
[status, cmdout] = system(cmd);
if status ~= 0
    error('Failed to refresh Nyxus clone.\nOutput:\n%s', cmdout);
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

function truth = load_named_constants_from_map(headerPath, mapName)
txt = fileread(headerPath);
pattern = sprintf('static\\s+std::unordered_map<std::string,\\s*double>\\s+%s\\s*\\{([\\s\\S]*?)\\};', mapName);
block = regexp(txt, pattern, 'tokens', 'once');
if isempty(block)
    error('Could not find map %s in %s', mapName, headerPath);
end

tokens = regexp(block{1}, '\{\s*"([A-Z0-9_]+)"\s*,\s*([-+0-9.eE]+)\s*\}', 'tokens');
if isempty(tokens)
    error('Could not parse named constants in %s', headerPath);
end

truth = struct();
for i = 1:numel(tokens)
    truth.(tokens{i}{1}) = str2double(tokens{i}{2});
end
end

function features = build_ngtdm_features()
features = {
    'NGTDM_COARSENESS'
    'NGTDM_CONTRAST'
    'NGTDM_BUSYNESS'
    'NGTDM_COMPLEXITY'
    'NGTDM_STRENGTH'
    };
end

function [checkedCount, passCount] = report_ngtdm(titleText, features, truth, sliceResults)
checkedCount = 0;
passCount = 0;
fracTolerance = 100.0;

fprintf('%s\n', titleText);
fprintf('%-24s %18s %18s %12s %8s\n', 'Feature', 'Expected', 'MATLAB', 'AbsDiff', 'Pass');
fprintf('%s\n', repmat('-', 1, 88));

for i = 1:numel(features)
    featureName = features{i};
    expected = truth.(featureName);
    actual = aggregate_feature(sliceResults, featureName);
    pass = agrees_gt_matlab(actual, expected, fracTolerance);

    checkedCount = checkedCount + 1;
    if pass
        passCount = passCount + 1;
    end

    fprintf('%-24s %18.12g %18.12g %12.6g %8s\n', ...
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

function f = compute_ngtdm_slice_features(auxImage, ibsiMode, nLevels)
if ibsiMode
    D = auxImage;
else
    D = bin_image_matlab(auxImage, max(auxImage(:)), nLevels);
end

uniqueLevels = unique(D(:));
uniqueLevels(uniqueLevels == 0) = [];
if ibsiMode
    levels = 0:max(uniqueLevels);
else
    levels = sort(uniqueLevels(:))';
end

if numel(levels) < 2
    error('NGTDM verifier expected at least two gray levels.');
end

N = zeros(1, numel(levels) + 1);
S = zeros(1, numel(levels) + 1);
Nvp = 0;

for row = 1:size(D, 1)
    for col = 1:size(D, 2)
        piVal = D(row, col);
        if piVal == 0
            continue;
        end

        [neighborSum, nd] = gather_ngtdm_neighbors(D, row, col);
        if nd == 0
            continue;
        end

        aveNeighbor = neighborSum / nd;
        if ibsiMode
            levelIdx = piVal + 1;
        else
            levelIdx = find(levels == piVal, 1);
        end

        N(levelIdx) = N(levelIdx) + 1;
        S(levelIdx) = S(levelIdx) + abs(levels(levelIdx) - aveNeighbor);
        if aveNeighbor > 0.0
            Nvp = Nvp + 1;
        end
    end
end

Ng = numel(levels);
Ngp = numel(uniqueLevels);
Nvc = sum(N);
P = N / Nvc;

f = struct();
f.NGTDM_COARSENESS = calc_coarseness(P, S, Ng);
f.NGTDM_CONTRAST = calc_contrast(P, S, levels, Ng, Ngp, Nvc);
f.NGTDM_BUSYNESS = calc_busyness(P, S, levels, Ng, Ngp);
f.NGTDM_COMPLEXITY = calc_complexity(P, S, levels, Ng, Nvp);
f.NGTDM_STRENGTH = calc_strength(P, S, levels, Ng);
end

function [neighborSum, nd] = gather_ngtdm_neighbors(D, row, col)
offsets = [-1 0; -1 1; 0 1; 1 1; 1 0; 1 -1; 0 -1; -1 -1];
neighborSum = 0.0;
nd = 0;
for k = 1:size(offsets, 1)
    rr = row + offsets(k, 1);
    cc = col + offsets(k, 2);
    if rr < 1 || rr > size(D, 1) || cc < 1 || cc > size(D, 2)
        continue;
    end
    if D(rr, cc) ~= 0
        neighborSum = neighborSum + D(rr, cc);
        nd = nd + 1;
    end
end
end

function D = bin_image_matlab(auxImage, maxI, nLevels)
if maxI <= 0
    error('Expected positive max intensity for NGTDM verifier.');
end

slope = double(nLevels) / double(maxI);
intercept = 1.0;
D = zeros(size(auxImage));
for idx = 1:numel(auxImage)
    if auxImage(idx) == 0
        D(idx) = 1;
    else
        y = floor(slope * double(auxImage(idx)) + intercept);
        if y > nLevels
            y = nLevels;
        end
        if y < 1
            y = 1;
        end
        D(idx) = y;
    end
end
end

function retval = calc_coarseness(P, S, Ng)
sumVal = 0.0;
for i = 1:Ng
    sumVal = sumVal + P(i) * S(i);
end
retval = 1.0 / sumVal;
end

function retval = calc_contrast(P, S, levels, Ng, Ngp, Nvc)
sumVal = 0.0;
for i = 1:Ng
    for j = 1:Ng
        sumVal = sumVal + P(i) * P(j) * (levels(i) - levels(j)) ^ 2;
    end
end
NgpP2 = Ngp;
if Ngp > 1
    NgpP2 = Ngp * (Ngp - 1);
end
term1 = sumVal / double(NgpP2);
term2 = sum(S(1:Ng)) / Nvc;
retval = term1 * term2;
end

function retval = calc_busyness(P, S, levels, Ng, Ngp)
if Ngp == 1
    retval = 0.0;
    return;
end

sum1 = 0.0;
for i = 1:Ng
    sum1 = sum1 + P(i) * S(i);
end

sum2 = 0.0;
for i = 1:Ng
    for j = 1:Ng
        if P(i) ~= 0 && P(j) ~= 0
            sum2 = sum2 + abs(P(i) * levels(i) - P(j) * levels(j));
        end
    end
end

if sum2 == 0
    retval = 0.0;
else
    retval = sum1 / sum2;
end
end

function retval = calc_complexity(P, S, levels, Ng, Nvp)
sumVal = 0.0;
for i = 1:Ng
    for j = 1:Ng
        if P(i) ~= 0 && P(j) ~= 0
            sumVal = sumVal + abs(levels(i) - levels(j)) * (P(i) * S(i) + P(j) * S(j)) / (P(i) + P(j));
        end
    end
end
retval = sumVal / double(Nvp);
end

function retval = calc_strength(P, S, levels, Ng)
sum1 = 0.0;
for i = 1:Ng
    for j = 1:Ng
        if P(i) ~= 0 && P(j) ~= 0
            sum1 = sum1 + (P(i) + P(j)) * (levels(i) - levels(j)) ^ 2;
        end
    end
end
sum2 = sum(S(1:Ng));
retval = sum1 / sum2;
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
