function verify_nyxus_ngldm_tests(repoRoot)
% Cross-verify Nyxus 2D NGLDM unit-test expectations in MATLAB.
%
% This script:
% 1. Clones the Nyxus repo under test if repoRoot is not provided
% 2. Reads the IBSI phantom slices from tests/test_data.h
% 3. Reads the asserted constants from tests/test_ibsi_ngldm.h
% 4. Recomputes all active 2D NGLDM regression features in MATLAB
% 5. Reports whether each value would pass Nyxus's C++ test tolerances

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
headerPath = fullfile(repoRoot, "tests", "test_ibsi_ngldm.h");

slices = load_ibsi_phantom_slices(testDataPath);
truth = load_named_constants_from_map(headerPath, 'ngldm_regression_gtruth');
features = build_ngldm_features();

results = cell(numel(slices), 1);
for i = 1:numel(slices)
    results{i} = compute_ngldm_slice_features(slices{i}.matrix);
end

fprintf('\nNyxus NGLDM verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Phantom source: %s\n', testDataPath);
fprintf('Test constants source: %s\n', headerPath);
fprintf('Slices in phantom: %d\n', numel(slices));
fprintf('Configuration: 2D NGLDM regression surface over the IBSI phantom slices using Nyxus''s current implementation semantics\n\n');

[checkedCount, passCount] = report_ngldm(features, truth, results);

fprintf('\nSummary\n');
fprintf('NGLDM: %d/%d features pass Nyxus''s C++ test tolerances.\n', passCount, checkedCount);
fprintf('\nNotes:\n');
fprintf('1. This verifier targets the active 2D NGLDM regression constants in tests/test_ibsi_ngldm.h.\n');
fprintf('2. The MATLAB code intentionally mirrors Nyxus''s current behavior, including counting zero-valued background pixels inside the ROI AABB when constructing the dependence matrix.\n');
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
    fullfile(repoRoot, "tests", "test_ibsi_ngldm.h")
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

function features = build_ngldm_features()
features = {
    'NGLDM_LDE'
    'NGLDM_HDE'
    'NGLDM_LGLCE'
    'NGLDM_HGLCE'
    'NGLDM_LDLGLE'
    'NGLDM_LDHGLE'
    'NGLDM_HDLGLE'
    'NGLDM_HDHGLE'
    'NGLDM_GLNU'
    'NGLDM_GLNUN'
    'NGLDM_DCNU'
    'NGLDM_DCNUN'
    'NGLDM_DCP'
    'NGLDM_GLM'
    'NGLDM_GLV'
    'NGLDM_DCM'
    'NGLDM_DCV'
    'NGLDM_DCENT'
    'NGLDM_DCENE'
    };
end

function [checkedCount, passCount] = report_ngldm(features, truth, sliceResults)
checkedCount = 0;
passCount = 0;
fracTolerance = 2.0;

fprintf('NGLDM\n');
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

function f = compute_ngldm_slice_features(auxImage)
maxInten = max(auxImage(:));
levels = unique(auxImage(:));
levels = sort(levels(:))';

N = zeros(numel(levels), 9);
maxDep = 0;
for row = 2:size(auxImage, 1) - 1
    for col = 2:size(auxImage, 2) - 1
        cpi = auxImage(row, col);
        rowIdx = find(levels == cpi, 1);
        matches = 0;
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0
                    continue;
                end
                if auxImage(row + dr, col + dc) == cpi
                    matches = matches + 1;
                end
            end
        end
        N(rowIdx, matches + 1) = N(rowIdx, matches + 1) + 1;
        maxDep = max(maxDep, matches);
    end
end

Nr = maxDep + 1;
N = N(:, 1:Nr);
Ng = size(N, 1);
Ns = sum(N(:));

f = struct();
f.NGLDM_LDE = 0.0;
f.NGLDM_HDE = 0.0;
f.NGLDM_LGLCE = 0.0;
f.NGLDM_HGLCE = 0.0;
f.NGLDM_LDLGLE = 0.0;
f.NGLDM_LDHGLE = 0.0;
f.NGLDM_HDLGLE = 0.0;
f.NGLDM_HDHGLE = 0.0;
f.NGLDM_GLNU = 0.0;
f.NGLDM_GLNUN = 0.0;
f.NGLDM_DCNU = 0.0;
f.NGLDM_DCNUN = 0.0;
f.NGLDM_DCP = 1.0;
f.NGLDM_GLM = 0.0;
f.NGLDM_GLV = 0.0;
f.NGLDM_DCM = 0.0;
f.NGLDM_DCV = 0.0;
f.NGLDM_DCENT = 0.0;
f.NGLDM_DCENE = 0.0;

for i = 1:Ng
    sj = 0.0;
    for j = 2:Nr
        grey = levels(i);
        sij = N(i, j);
        k = j;
        depIdx = j - 1;
        pij = sij / Ns;

        f.NGLDM_LDE = f.NGLDM_LDE + sij / (depIdx * depIdx);
        f.NGLDM_HDE = f.NGLDM_HDE + sij * depIdx * depIdx;
        if grey ~= 0
            f.NGLDM_LGLCE = f.NGLDM_LGLCE + sij / (grey * grey);
        end
        f.NGLDM_HGLCE = f.NGLDM_HGLCE + sij * grey * grey;
        if grey ~= 0
            f.NGLDM_LDLGLE = f.NGLDM_LDLGLE + sij / (depIdx * depIdx * grey * grey);
            f.NGLDM_HDLGLE = f.NGLDM_HDLGLE + sij * k * k / (grey * grey);
        end
        f.NGLDM_LDHGLE = f.NGLDM_LDHGLE + sij * grey * grey / (k * k);
        f.NGLDM_HDHGLE = f.NGLDM_HDHGLE + sij * k * k * grey * grey;
        f.NGLDM_GLM = f.NGLDM_GLM + grey * pij;
        f.NGLDM_DCM = f.NGLDM_DCM + k * pij;
        if pij > 0
            f.NGLDM_DCENT = f.NGLDM_DCENT - pij * log2(pij);
        end
        f.NGLDM_DCENE = f.NGLDM_DCENE + pij * pij;
        sj = sj + sij;
    end
    f.NGLDM_GLNU = f.NGLDM_GLNU + sj * sj;
    f.NGLDM_GLNUN = f.NGLDM_GLNUN + sj * sj;
end

for i = 1:Ng
    si = 0.0;
    for j = 2:Nr
        si = si + N(i, j);
    end
    f.NGLDM_DCNU = f.NGLDM_DCNU + si * si;
    f.NGLDM_DCNUN = f.NGLDM_DCNUN + si * si;
end

for i = 1:Ng
    for j = 2:Nr
        sij = N(i, j);
        k = j;
        pij = sij / Ns;
        rowIndex1Based = i;
        f.NGLDM_GLV = f.NGLDM_GLV + (rowIndex1Based - f.NGLDM_GLM) ^ 2 * pij;
        f.NGLDM_DCV = f.NGLDM_DCV + (k - f.NGLDM_DCM) ^ 2 * pij;
    end
end

f.NGLDM_LDE = f.NGLDM_LDE / Ns;
f.NGLDM_HDE = f.NGLDM_HDE / Ns;
f.NGLDM_LGLCE = f.NGLDM_LGLCE / Ns;
f.NGLDM_HGLCE = f.NGLDM_HGLCE / Ns;
f.NGLDM_LDLGLE = f.NGLDM_LDLGLE / Ns;
f.NGLDM_LDHGLE = f.NGLDM_LDHGLE / Ns;
f.NGLDM_HDLGLE = f.NGLDM_HDLGLE / Ns;
f.NGLDM_HDHGLE = f.NGLDM_HDHGLE / Ns;
f.NGLDM_GLNU = f.NGLDM_GLNU / Ns;
f.NGLDM_GLNUN = f.NGLDM_GLNUN / (Ns * Ns);
f.NGLDM_DCNU = f.NGLDM_DCNU / Ns;
f.NGLDM_DCNUN = f.NGLDM_DCNUN / (Ns * Ns);
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
