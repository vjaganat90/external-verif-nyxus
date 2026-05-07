function verify_nyxus_gldzm_tests(repoRoot)
% Cross-verify Nyxus 2D GLDZM unit-test expectations in MATLAB.
%
% This script:
% 1. Clones the Nyxus repo under test if repoRoot is not provided
% 2. Reads the IBSI phantom slices from tests/test_data.h
% 3. Reads the asserted constants from tests/test_ibsi_gldzm.h
% 4. Recomputes all 2D GLDZM features in MATLAB
% 5. Reports whether each value would pass Nyxus's C++ test tolerances

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
headerPath = fullfile(repoRoot, "tests", "test_ibsi_gldzm.h");

slices = load_ibsi_phantom_slices(testDataPath);
truth = load_named_constants_from_map(headerPath, 'ibsi_gldzm_gtruth');
features = build_gldzm_features();

results = cell(numel(slices), 1);
for i = 1:numel(slices)
    results{i} = compute_gldzm_slice_features(slices{i}.matrix);
end

fprintf('\nNyxus GLDZM verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Phantom source: %s\n', testDataPath);
fprintf('Test constants source: %s\n', headerPath);
fprintf('Slices in phantom: %d\n', numel(slices));
fprintf('Configuration: 2D IBSI-style GLDZM over the shared digital phantom slices\n\n');

[checkedCount, passCount] = report_gldzm(features, truth, results);

fprintf('\nSummary\n');
fprintf('GLDZM: %d/%d features pass Nyxus''s C++ test tolerances.\n', passCount, checkedCount);
fprintf('\nNotes:\n');
fprintf('1. This verifier targets the active 2D GLDZM constants in tests/test_ibsi_gldzm.h.\n');
fprintf('2. The MATLAB code independently reconstructs same-intensity zones and 4-connected border distances over the masked ROI matrices.\n');
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
    fullfile(repoRoot, "tests", "test_ibsi_gldzm.h")
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

function features = build_gldzm_features()
features = {
    'GLDZM_SDE'
    'GLDZM_LDE'
    'GLDZM_LGLZE'
    'GLDZM_HGLZE'
    'GLDZM_SDLGLE'
    'GLDZM_SDHGLE'
    'GLDZM_LDLGLE'
    'GLDZM_LDHGLE'
    'GLDZM_GLNU'
    'GLDZM_GLNUN'
    'GLDZM_ZDNU'
    'GLDZM_ZDNUN'
    'GLDZM_ZP'
    'GLDZM_GLM'
    'GLDZM_GLV'
    'GLDZM_ZDM'
    'GLDZM_ZDV'
    'GLDZM_ZDE'
    };
end

function [checkedCount, passCount] = report_gldzm(features, truth, sliceResults)
checkedCount = 0;
passCount = 0;
fracTolerance = 2.0;

fprintf('GLDZM\n');
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

function f = compute_gldzm_slice_features(auxImage)
D = auxImage;
maxLevel = max(D(:));
levels = 1:maxLevel;
visited = false(size(D));
zones = zeros(0, 3);

for row = 1:size(D, 1)
    for col = 1:size(D, 2)
        inten = D(row, col);
        if inten == 0 || visited(row, col)
            continue;
        end

        [coords, visited] = flood_zone(D, visited, row, col, inten);
        zoneMetric = inf;
        for k = 1:size(coords, 1)
            zoneMetric = min(zoneMetric, dist2border(D, coords(k, 2), coords(k, 1)));
        end
        zones(end + 1, :) = [inten, zoneMetric, size(coords, 1)]; %#ok<AGROW>
    end
end

Nd = max(zones(:, 2));
P = zeros(numel(levels), Nd);
for i = 1:size(zones, 1)
    P(zones(i, 1), zones(i, 2)) = P(zones(i, 1), zones(i, 2)) + 1;
end

Mx = sum(P, 2);
Md = sum(P, 1);
Ns = sum(P(:));
Nv = nnz(D);

f = struct();
f.GLDZM_SDE = 0.0;
f.GLDZM_LDE = 0.0;
f.GLDZM_LGLZE = 0.0;
f.GLDZM_HGLZE = 0.0;
f.GLDZM_SDLGLE = 0.0;
f.GLDZM_SDHGLE = 0.0;
f.GLDZM_LDLGLE = 0.0;
f.GLDZM_LDHGLE = 0.0;
f.GLDZM_GLNU = 0.0;
f.GLDZM_GLNUN = 0.0;
f.GLDZM_ZDNU = 0.0;
f.GLDZM_ZDNUN = 0.0;
f.GLDZM_ZP = Ns / Nv;
f.GLDZM_GLM = 0.0;
f.GLDZM_GLV = 0.0;
f.GLDZM_ZDM = 0.0;
f.GLDZM_ZDV = 0.0;
f.GLDZM_ZDE = 0.0;

for d = 1:Nd
    f.GLDZM_SDE = f.GLDZM_SDE + Md(d) / (d * d);
    f.GLDZM_LDE = f.GLDZM_LDE + Md(d) * (d * d);
    f.GLDZM_ZDNU = f.GLDZM_ZDNU + Md(d) * Md(d);
end
f.GLDZM_SDE = f.GLDZM_SDE / Ns;
f.GLDZM_LDE = f.GLDZM_LDE / Ns;
f.GLDZM_ZDNU = f.GLDZM_ZDNU / Ns;
f.GLDZM_ZDNUN = f.GLDZM_ZDNU / Ns;

for g = 1:numel(levels)
    grey = levels(g);
    f.GLDZM_LGLZE = f.GLDZM_LGLZE + Mx(g) / (grey * grey);
    f.GLDZM_HGLZE = f.GLDZM_HGLZE + (grey * grey) * Mx(g);
    f.GLDZM_GLNU = f.GLDZM_GLNU + Mx(g) * Mx(g);
end
f.GLDZM_LGLZE = f.GLDZM_LGLZE / Ns;
f.GLDZM_HGLZE = f.GLDZM_HGLZE / Ns;
f.GLDZM_GLNU = f.GLDZM_GLNU / Ns;
f.GLDZM_GLNUN = f.GLDZM_GLNU / Ns;

for g = 1:numel(levels)
    grey = levels(g);
    for d = 1:Nd
        pRaw = P(g, d);
        f.GLDZM_SDLGLE = f.GLDZM_SDLGLE + pRaw / (grey * grey * d * d);
        f.GLDZM_SDHGLE = f.GLDZM_SDHGLE + (grey * grey) * pRaw / (d * d);
        f.GLDZM_LDLGLE = f.GLDZM_LDLGLE + (d * d) * pRaw / (grey * grey);
        f.GLDZM_LDHGLE = f.GLDZM_LDHGLE + (grey * grey) * (d * d) * pRaw;
        f.GLDZM_GLM = f.GLDZM_GLM + grey * pRaw;
        f.GLDZM_ZDM = f.GLDZM_ZDM + d * pRaw;
        if pRaw ~= 0
            pNorm = pRaw / Ns;
            f.GLDZM_ZDE = f.GLDZM_ZDE + pNorm * log2(pNorm + 2.2e-16);
        end
    end
end
f.GLDZM_SDLGLE = f.GLDZM_SDLGLE / Ns;
f.GLDZM_SDHGLE = f.GLDZM_SDHGLE / Ns;
f.GLDZM_LDLGLE = f.GLDZM_LDLGLE / Ns;
f.GLDZM_LDHGLE = f.GLDZM_LDHGLE / Ns;
f.GLDZM_GLM = f.GLDZM_GLM / Ns;
f.GLDZM_ZDM = f.GLDZM_ZDM / Ns;
f.GLDZM_ZDE = -f.GLDZM_ZDE;

for g = 1:numel(levels)
    grey = levels(g);
    for d = 1:Nd
        pNorm = P(g, d) / Ns;
        f.GLDZM_GLV = f.GLDZM_GLV + pNorm * (grey - f.GLDZM_GLM) ^ 2;
        f.GLDZM_ZDV = f.GLDZM_ZDV + pNorm * (d - f.GLDZM_ZDM) ^ 2;
    end
end
end

function [coords, visited] = flood_zone(D, visited, row0, col0, inten)
queue = zeros(numel(D), 2);
head = 1;
tail = 1;
queue(tail, :) = [row0, col0];
visited(row0, col0) = true;
coords = zeros(0, 2);

while head <= tail
    row = queue(head, 1);
    col = queue(head, 2);
    head = head + 1;
    coords(end + 1, :) = [row, col]; %#ok<AGROW>

    neighbors = [row, col + 1; row + 1, col; row, col - 1; row - 1, col];
    for i = 1:4
        r = neighbors(i, 1);
        c = neighbors(i, 2);
        if r < 1 || r > size(D, 1) || c < 1 || c > size(D, 2)
            continue;
        end
        if visited(r, c) || D(r, c) ~= inten
            continue;
        end
        tail = tail + 1;
        queue(tail, :) = [r, c];
        visited(r, c) = true;
    end
end
end

function retval = dist2border(I, x, y)
[h, w] = size(I);

dist2l = 0;
for x0 = x - 1:-1:1
    if I(y, x0) == 0 || x0 == 1
        dist2l = x - x0;
        break;
    end
end

dist2r = 0;
for x0 = x + 1:w
    if I(y, x0) == 0 || x0 == w
        dist2r = x0 - x;
        break;
    end
end

dist2t = 0;
for y0 = y - 1:-1:1
    if I(y0, x) == 0 || y0 == 1
        dist2t = y - y0;
        break;
    end
end

dist2b = 0;
for y0 = y + 1:h
    if I(y0, x) == 0 || y0 == h
        dist2b = y0 - y;
        break;
    end
end

retval = min([dist2l + 1, dist2r + 1, dist2t + 1, dist2b + 1]);
if retval == 0
    retval = 1;
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
