function verify_nyxus_glszm_tests(repoRoot)
% Cross-verify Nyxus 2D GLSZM unit-test expectations in MATLAB.
%
% This script:
% 1. Clones the Nyxus repo under test if repoRoot is not provided
% 2. Reads the IBSI phantom slices from tests/test_data.h
% 3. Reads the asserted constants from tests/test_glszm.h
% 4. Recomputes all 2D GLSZM features in MATLAB
% 5. Reports whether each value would pass Nyxus's C++ test tolerances

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
headerPath = fullfile(repoRoot, "tests", "test_glszm.h");

slices = load_ibsi_phantom_slices(testDataPath);
truth = load_named_constants_from_map(headerPath, 'glszm_values');
features = build_glszm_features();

results = cell(numel(slices), 1);
for i = 1:numel(slices)
    results{i} = compute_glszm_slice_features(slices{i}.matrix);
end

fprintf('\nNyxus GLSZM verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Phantom source: %s\n', testDataPath);
fprintf('Test constants source: %s\n', headerPath);
fprintf('Slices in phantom: %d\n', numel(slices));
fprintf('Configuration: 2D non-IBSI GLSZM with 64-level MATLAB-style binning over the shared digital phantom slices\n\n');

[checkedCount, passCount] = report_glszm(features, truth, results);

fprintf('\nSummary\n');
fprintf('GLSZM: %d/%d features pass Nyxus''s C++ test tolerances.\n', passCount, checkedCount);
fprintf('\nNotes:\n');
fprintf('1. This verifier targets the active 2D GLSZM constants in tests/test_glszm.h.\n');
fprintf('2. The MATLAB code independently rebuilds 8-connected same-intensity zones after Nyxus-style 64-level MATLAB binning.\n');
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
    fullfile(repoRoot, "tests", "test_glszm.h")
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

function features = build_glszm_features()
features = {
    'GLSZM_SAE'
    'GLSZM_LAE'
    'GLSZM_LGLZE'
    'GLSZM_HGLZE'
    'GLSZM_SALGLE'
    'GLSZM_SAHGLE'
    'GLSZM_LALGLE'
    'GLSZM_LAHGLE'
    'GLSZM_GLN'
    'GLSZM_GLNN'
    'GLSZM_SZN'
    'GLSZM_SZNN'
    'GLSZM_ZP'
    'GLSZM_GLV'
    'GLSZM_ZV'
    'GLSZM_ZE'
    };
end

function [checkedCount, passCount] = report_glszm(features, truth, sliceResults)
checkedCount = 0;
passCount = 0;
fracTolerance = 100.0;

fprintf('GLSZM\n');
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

function f = compute_glszm_slice_features(auxImage)
nLevels = 64;
D = bin_image_matlab(auxImage, max(auxImage(:)), nLevels);
levels = unique(D(:));
levels(levels == 0) = [];
levels = sort(levels(:))';
zones = scan_zones_forward(D);
P = zeros(numel(levels), numel(D));
for i = 1:size(zones, 1)
    rowIdx = find(levels == zones(i, 1), 1);
    colIdx = zones(i, 2);
    P(rowIdx, colIdx) = P(rowIdx, colIdx) + 1;
end

Ng = numel(levels);
Ns = size(D, 1) * size(D, 2);
sumP = sum(P(:));
Np = nnz(D);

si = zeros(1, Ng);
sj = zeros(1, Ns);
muGLV = 0.0;
muZV = 0.0;
sumSALGLE = 0.0;
sumSAHGLE = 0.0;
sumLALGLE = 0.0;
sumLAHGLE = 0.0;
zoneEntropy = 0.0;

for i = 1:Ng
    inten = levels(i);
    for j = 1:Ns
        pRaw = P(i, j);
        si(i) = si(i) + pRaw;
        sj(j) = sj(j) + pRaw;
        pNorm = pRaw / sumP;
        i2 = inten * inten;
        j2 = j * j;
        sumLAHGLE = sumLAHGLE + pRaw * i2 * j2;
        sumLALGLE = sumLALGLE + pRaw * j2 / i2;
        sumSAHGLE = sumSAHGLE + pRaw * i2 / j2;
        sumSALGLE = sumSALGLE + pRaw / (i2 * j2);
        zoneEntropy = zoneEntropy + pNorm * log2(pNorm + 2.2e-16);
        muZV = muZV + pNorm * j;
        muGLV = muGLV + pNorm * inten;
    end
end

f = struct();
f.GLSZM_SAE = sum(sj ./ ((1:Ns) .^ 2)) / sumP;
f.GLSZM_LAE = sum(sj .* ((1:Ns) .^ 2)) / sumP;
f.GLSZM_GLN = sum(si .^ 2) / sumP;
f.GLSZM_GLNN = sum(si .^ 2) / (sumP * sumP);
f.GLSZM_SZN = sum(sj .^ 2) / sumP;
f.GLSZM_SZNN = sum(sj .^ 2) / (sumP * sumP);
f.GLSZM_ZP = sumP / double(Np);
f.GLSZM_GLV = 0.0;
f.GLSZM_ZV = 0.0;
f.GLSZM_ZE = -zoneEntropy;
f.GLSZM_LGLZE = sum(si ./ (levels .^ 2)) / sumP;
f.GLSZM_HGLZE = sum(si .* (levels .^ 2)) / sumP;
f.GLSZM_SALGLE = sumSALGLE / sumP;
f.GLSZM_SAHGLE = sumSAHGLE / sumP;
f.GLSZM_LALGLE = sumLALGLE / sumP;
f.GLSZM_LAHGLE = sumLAHGLE / sumP;

for i = 1:Ng
    inten = levels(i);
    for j = 1:Ns
        pNorm = P(i, j) / sumP;
        f.GLSZM_GLV = f.GLSZM_GLV + pNorm * (inten - muGLV) ^ 2;
        f.GLSZM_ZV = f.GLSZM_ZV + pNorm * (j - muZV) ^ 2;
    end
end
end

function D = bin_image_matlab(auxImage, maxI, nLevels)
if maxI <= 0
    error('Expected positive max intensity for GLSZM verifier.');
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

function zones = scan_zones_forward(Din)
% Mirror Nyxus's current non-IBSI GLSZM traversal:
% 1. Matlab binning maps zero-valued background pixels to gray level 1.
% 2. The zone walk explores E, SE, S, SW, then backtracks through history.
VISITED = -1;
D = Din;
zones = zeros(0, 2);

for row = 1:size(D, 1)
    for col = 1:size(D, 2)
        inten = D(row, col);
        if inten == 0 || inten == VISITED
            continue;
        end

        x = col;
        y = row;
        zoneArea = 1;
        history = zeros(numel(D), 2);
        historySize = 0;
        D(y, x) = VISITED;

        while true
            if is_same_unvisited(D, y, x + 1, inten, VISITED)
                historySize = historySize + 1;
                history(historySize, :) = [x, y];
                x = x + 1;
                D(y, x) = VISITED;
                zoneArea = zoneArea + 1;
                continue;
            end
            if is_same_unvisited(D, y + 1, x + 1, inten, VISITED)
                historySize = historySize + 1;
                history(historySize, :) = [x, y];
                x = x + 1;
                y = y + 1;
                D(y, x) = VISITED;
                zoneArea = zoneArea + 1;
                continue;
            end
            if is_same_unvisited(D, y + 1, x, inten, VISITED)
                historySize = historySize + 1;
                history(historySize, :) = [x, y];
                y = y + 1;
                D(y, x) = VISITED;
                zoneArea = zoneArea + 1;
                continue;
            end
            if is_same_unvisited(D, y + 1, x - 1, inten, VISITED)
                historySize = historySize + 1;
                history(historySize, :) = [x, y];
                x = x - 1;
                y = y + 1;
                D(y, x) = VISITED;
                zoneArea = zoneArea + 1;
                continue;
            end
            if historySize > 0
                x = history(historySize, 1);
                y = history(historySize, 2);
                historySize = historySize - 1;
                continue;
            end
            break;
        end

        zones(end + 1, :) = [inten, zoneArea]; %#ok<AGROW>
    end
end
end

function ok = is_same_unvisited(D, row, col, inten, visitedValue)
ok = row >= 1 && row <= size(D, 1) && ...
     col >= 1 && col <= size(D, 2) && ...
     D(row, col) ~= visitedValue && ...
     D(row, col) == inten;
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
