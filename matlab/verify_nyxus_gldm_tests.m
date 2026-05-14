function verify_nyxus_gldm_tests(repoRoot)
% Cross-verify Nyxus 2D GLDM unit-test expectations in MATLAB.
%
% This script:
% 1. Clones the Nyxus repo under test if repoRoot is not provided
% 2. Reads the IBSI phantom slices from tests/test_data.h
% 3. Reads the asserted constants from tests/test_ibsi_gldm.h
% 4. Recomputes all 2D IBSI GLDM features in MATLAB
% 5. Reports whether each value would pass Nyxus's C++ test tolerances

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
headerPath = fullfile(repoRoot, "tests", "test_ibsi_gldm.h");

slices = load_ibsi_phantom_slices(testDataPath);
truth = load_named_constants_from_map(headerPath, 'IBSI_gldm_values');
features = build_gldm_features();

results = cell(numel(slices), 1);
for i = 1:numel(slices)
    results{i} = compute_gldm_slice_features(slices{i}.matrix);
end

fprintf('\nNyxus GLDM verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Phantom source: %s\n', testDataPath);
fprintf('Test constants source: %s\n', headerPath);
fprintf('Slices in phantom: %d\n', numel(slices));
fprintf('Configuration: 2D IBSI GLDM over the shared digital phantom slices, averaged across slices\n\n');

[checkedCount, passCount] = report_gldm(features, truth, results);

fprintf('\nSummary\n');
fprintf('GLDM: %d/%d features pass Nyxus''s C++ test tolerances.\n', passCount, checkedCount);
fprintf('\nNotes:\n');
fprintf('1. This verifier targets the active 2D IBSI GLDM constants in tests/test_ibsi_gldm.h.\n');
fprintf('2. The MATLAB code independently rebuilds the gray-level dependence matrix by counting 8-neighbor same-intensity dependencies for each nonzero ROI pixel.\n');
fprintf('3. This is a reimplementation of the matrix and feature formulas over the same phantom, not a call back into Nyxus.\n\n');

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
    fullfile(repoRoot, "tests", "test_ibsi_gldm.h")
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

function features = build_gldm_features()
features = {
    'GLDM_SDE'
    'GLDM_LDE'
    'GLDM_LGLE'
    'GLDM_HGLE'
    'GLDM_SDLGLE'
    'GLDM_SDHGLE'
    'GLDM_LDLGLE'
    'GLDM_LDHGLE'
    'GLDM_GLN'
    'GLDM_DN'
    'GLDM_DNN'
    'GLDM_GLV'
    'GLDM_DV'
    'GLDM_DE'
    };
end

function [checkedCount, passCount] = report_gldm(features, truth, sliceResults)
checkedCount = 0;
passCount = 0;
fracTolerance = 100.0;

fprintf('GLDM\n');
fprintf('%-22s %18s %18s %12s %8s\n', 'Feature', 'Expected', 'MATLAB', 'AbsDiff', 'Pass');
fprintf('%s\n', repmat('-', 1, 86));

for i = 1:numel(features)
    featureName = features{i};
    expected = truth.(featureName);
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

function f = compute_gldm_slice_features(auxImage)
D = auxImage;
maxLevel = max(D(:));
levels = (1:maxLevel)';
deps = 1:9;
P = zeros(numel(levels), numel(deps));

for row = 1:size(D, 1)
    for col = 1:size(D, 2)
        grey = D(row, col);
        if grey == 0
            continue;
        end

        depCount = 1;
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0
                    continue;
                end
                rr = row + dr;
                cc = col + dc;
                if rr < 1 || rr > size(D, 1) || cc < 1 || cc > size(D, 2)
                    continue;
                end
                if D(rr, cc) == grey
                    depCount = depCount + 1;
                end
            end
        end

        P(grey, depCount) = P(grey, depCount) + 1;
    end
end

Nz = sum(P(:));
if Nz == 0
    error('GLDM matrix is empty.');
end

greyWeights = repmat(levels, 1, numel(deps));
depWeights = repmat(deps, numel(levels), 1);
prob = P / Nz;
greyMarginal = sum(P, 2);
depMarginal = sum(P, 1);

greyMean = sum(sum(prob .* greyWeights));
depMean = sum(sum(prob .* depWeights));

f = struct();
f.GLDM_SDE = sum(sum(P ./ (depWeights .^ 2))) / Nz;
f.GLDM_LDE = sum(sum(P .* (depWeights .^ 2))) / Nz;
f.GLDM_LGLE = sum(greyMarginal ./ (levels .^ 2)) / Nz;
f.GLDM_HGLE = sum(greyMarginal .* (levels .^ 2)) / Nz;
f.GLDM_SDLGLE = sum(sum(P ./ ((greyWeights .^ 2) .* (depWeights .^ 2)))) / Nz;
f.GLDM_SDHGLE = sum(sum(P .* (greyWeights .^ 2) ./ (depWeights .^ 2))) / Nz;
f.GLDM_LDLGLE = sum(sum(P .* (depWeights .^ 2) ./ (greyWeights .^ 2))) / Nz;
f.GLDM_LDHGLE = sum(sum(P .* (greyWeights .^ 2) .* (depWeights .^ 2))) / Nz;
f.GLDM_GLN = sum(greyMarginal .^ 2) / Nz;
f.GLDM_DN = sum(depMarginal .^ 2) / Nz;
f.GLDM_DNN = sum(depMarginal .^ 2) / (Nz * Nz);
f.GLDM_GLV = sum(sum(prob .* ((greyWeights - greyMean) .^ 2)));
f.GLDM_DV = sum(sum(prob .* ((depWeights - depMean) .^ 2)));
f.GLDM_DE = -sum(sum(prob .* log2(prob + 2.2e-16)));
end

function pass = agrees_gt_matlab(actual, expected, fracTolerance)
tolerance = abs(expected / fracTolerance);
pass = abs(actual - expected) <= tolerance;
end

function text = pass_fail(pass)
if pass
    text = 'PASS';
else
    text = 'FAIL';
end
end
