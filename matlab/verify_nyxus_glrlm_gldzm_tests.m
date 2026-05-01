function verify_nyxus_glrlm_gldzm_tests(repoRoot)
% Cross-verify Nyxus GLRLM and GLDZM unit-test expectations in MATLAB.
%
% This script:
% 1. Clones the Nyxus branch under test if repoRoot is not provided
% 2. Reads the IBSI phantom slices from tests/test_data.h
% 3. Reads the asserted constants from tests/test_glrlm.h and tests/test_ibsi_gldzm.h
% 4. Recomputes all GLRLM and GLDZM features in MATLAB
% 5. Reports whether each value would pass Nyxus's agrees_gt() tolerance

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
glrlmHeaderPath = fullfile(repoRoot, "tests", "test_glrlm.h");
gldzmHeaderPath = fullfile(repoRoot, "tests", "test_ibsi_gldzm.h");

slices = load_ibsi_phantom_slices(testDataPath);
glrlmTruth = load_named_constants(glrlmHeaderPath);
gldzmTruth = load_named_constants(gldzmHeaderPath);

glrlmCases = build_glrlm_cases();
gldzmFeatures = build_gldzm_features();

glrlmResults = cell(numel(slices), 1);
gldzmResults = cell(numel(slices), 1);
for i = 1:numel(slices)
    glrlmResults{i} = compute_glrlm_slice_features(slices{i}.matrix);
    gldzmResults{i} = compute_gldzm_slice_features(slices{i}.matrix);
end

fprintf('\nNyxus GLRLM/GLDZM verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Phantom source: %s\n', testDataPath);
fprintf('GLRLM test constants source: %s\n', glrlmHeaderPath);
fprintf('GLDZM test constants source: %s\n', gldzmHeaderPath);
fprintf('Slices in phantom: %d\n', numel(slices));
fprintf('GLRLM grey-level configuration: MATLAB-style binning with 100 levels\n');
fprintf('GLDZM grey-level configuration: IBSI mode (no binning)\n\n');

[glrlmChecked, glrlmPassed] = report_glrlm(glrlmCases, glrlmTruth, glrlmResults);
fprintf('\n');
[gldzmChecked, gldzmPassed] = report_gldzm(gldzmFeatures, gldzmTruth, gldzmResults);

fprintf('\nSummary\n');
fprintf('GLRLM: %d/%d features pass Nyxus''s C++ test tolerances.\n', glrlmPassed, glrlmChecked);
fprintf('GLDZM: %d/%d features pass Nyxus''s C++ test tolerances.\n', gldzmPassed, gldzmChecked);
fprintf('\nNotes:\n');
fprintf('1. GLRLM _AVE features are cross-checked the same way Nyxus tests them: they share the same ground-truth values as the base 4-angle metrics and average over angles per slice.\n');
fprintf('2. The MATLAB code intentionally reproduces Nyxus quirks that matter for test parity, including MATLAB-style 100-level binning for GLRLM and 4-connected zone distances for GLDZM.\n');
fprintf('3. This is an independent reimplementation of the matrices and feature formulas over the same test phantom, not a call back into Nyxus.\n\n');

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
requiredFiles = {
    fullfile(repoRoot, "tests", "test_data.h")
    fullfile(repoRoot, "tests", "test_glrlm.h")
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
gitBranch = 'unittest_add';
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

function cases = build_glrlm_cases()
baseFeatures = {
    'GLRLM_SRE'
    'GLRLM_LRE'
    'GLRLM_LGLRE'
    'GLRLM_HGLRE'
    'GLRLM_SRLGLE'
    'GLRLM_SRHGLE'
    'GLRLM_LRLGLE'
    'GLRLM_LRHGLE'
    'GLRLM_GLN'
    'GLRLM_GLNN'
    'GLRLM_RLN'
    'GLRLM_RLNN'
    'GLRLM_RP'
    'GLRLM_GLV'
    'GLRLM_RV'
    'GLRLM_RE'
    };

cases = cell(numel(baseFeatures) * 2, 1);
for i = 1:numel(baseFeatures)
    cases{i} = struct('feature', baseFeatures{i}, 'truthKey', baseFeatures{i});
    cases{numel(baseFeatures) + i} = struct( ...
        'feature', [baseFeatures{i} '_AVE'], ...
        'truthKey', baseFeatures{i});
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

function truth = load_named_constants(headerPath)
txt = fileread(headerPath);
tokens = regexp(txt, '\{\s*"([A-Z0-9_]+)"\s*,\s*([-+0-9.eE]+)\s*\}', 'tokens');
if isempty(tokens)
    error('Could not parse named constants in %s', headerPath);
end

truth = struct();
for i = 1:numel(tokens)
    truth.(tokens{i}{1}) = str2double(tokens{i}{2});
end
end

function [checkedCount, passCount] = report_glrlm(cases, truth, sliceResults)
checkedCount = 0;
passCount = 0;
fracTolerance = 100.0;

fprintf('GLRLM\n');
fprintf('%-22s %18s %18s %12s %8s\n', 'Feature', 'Expected', 'MATLAB', 'AbsDiff', 'Pass');
fprintf('%s\n', repmat('-', 1, 86));

for i = 1:numel(cases)
    featureName = cases{i}.feature;
    truthKey = cases{i}.truthKey;
    if ~isfield(truth, truthKey)
        error('Missing GLRLM ground truth for %s in tests/test_glrlm.h', truthKey);
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

function [checkedCount, passCount] = report_gldzm(features, truth, sliceResults)
checkedCount = 0;
passCount = 0;
fracTolerance = 2.0;

fprintf('GLDZM\n');
fprintf('%-22s %18s %18s %12s %8s\n', 'Feature', 'Expected', 'MATLAB', 'AbsDiff', 'Pass');
fprintf('%s\n', repmat('-', 1, 86));

for i = 1:numel(features)
    featureName = features{i};
    if ~isfield(truth, featureName)
        error('Missing GLDZM ground truth for %s in tests/test_ibsi_gldzm.h', featureName);
    end

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

function f = compute_glrlm_slice_features(auxImage)
nLevels = 100;
D = bin_image_matlab(auxImage, max(auxImage(:)), nLevels);
countNonZero = nnz(auxImage ~= 0);

baseNames = {
    'GLRLM_SRE'
    'GLRLM_LRE'
    'GLRLM_LGLRE'
    'GLRLM_HGLRE'
    'GLRLM_SRLGLE'
    'GLRLM_SRHGLE'
    'GLRLM_LRLGLE'
    'GLRLM_LRHGLE'
    'GLRLM_GLN'
    'GLRLM_GLNN'
    'GLRLM_RLN'
    'GLRLM_RLNN'
    'GLRLM_RP'
    'GLRLM_GLV'
    'GLRLM_RV'
    'GLRLM_RE'
    };

angled = struct();
for i = 1:numel(baseNames)
    angled.(baseNames{i}) = zeros(1, 4);
end

for angleIdx = 1:4
    [P, levels] = build_glrlm_matrix(D, angleIdx);
    sumP = sum(P(:));
    if sumP == 0
        continue;
    end

    runTotals = sum(P, 1);
    greyTotals = sum(P, 2);
    nr = size(P, 2);
    ng = size(P, 1);
    j = 1:nr;
    iLevels = levels(:);

    angled.GLRLM_SRE(angleIdx) = sum(runTotals ./ (j .^ 2)) / sumP;
    angled.GLRLM_LRE(angleIdx) = sum(runTotals .* (j .^ 2)) / sumP;
    angled.GLRLM_GLN(angleIdx) = sum(greyTotals .^ 2) / sumP;
    angled.GLRLM_GLNN(angleIdx) = sum(greyTotals .^ 2) / (sumP * sumP);
    angled.GLRLM_RLN(angleIdx) = sum(runTotals .^ 2) / sumP;
    angled.GLRLM_RLNN(angleIdx) = sum(runTotals .^ 2) / (sumP * sumP);
    angled.GLRLM_RP(angleIdx) = sumP / countNonZero;

    muGrey = 0.0;
    muRun = 0.0;
    entropyVal = 0.0;
    lglre = 0.0;
    hglre = 0.0;
    srlgle = 0.0;
    srhgle = 0.0;
    lrlgle = 0.0;
    lrhgle = 0.0;
    for r = 1:ng
        inten = iLevels(r);
        for c = 1:nr
            p = P(r, c) / sumP;
            muGrey = muGrey + p * inten;
            muRun = muRun + p * c;
            entropyVal = entropyVal + p * log2(p + 2.2e-16);
            lglre = lglre + P(r, c) / (inten * inten);
            hglre = hglre + P(r, c) * (inten * inten);
            srlgle = srlgle + P(r, c) / (inten * inten * c * c);
            srhgle = srhgle + P(r, c) * (inten * inten) / (c * c);
            lrlgle = lrlgle + P(r, c) * (c * c) / (inten * inten);
            lrhgle = lrhgle + P(r, c) * (inten * inten * c * c);
        end
    end

    glv = 0.0;
    rv = 0.0;
    for r = 1:ng
        inten = iLevels(r);
        for c = 1:nr
            p = P(r, c) / sumP;
            glv = glv + p * (inten - muGrey) ^ 2;
            rv = rv + p * (c - muRun) ^ 2;
        end
    end

    angled.GLRLM_GLV(angleIdx) = glv;
    angled.GLRLM_RV(angleIdx) = rv;
    angled.GLRLM_RE(angleIdx) = -entropyVal;
    angled.GLRLM_LGLRE(angleIdx) = lglre / sumP;
    angled.GLRLM_HGLRE(angleIdx) = hglre / sumP;
    angled.GLRLM_SRLGLE(angleIdx) = srlgle / sumP;
    angled.GLRLM_SRHGLE(angleIdx) = srhgle / sumP;
    angled.GLRLM_LRLGLE(angleIdx) = lrlgle / sumP;
    angled.GLRLM_LRHGLE(angleIdx) = lrhgle / sumP;
end

f = struct();
for i = 1:numel(baseNames)
    baseName = baseNames{i};
    f.(baseName) = mean(angled.(baseName));
    f.([baseName '_AVE']) = mean(angled.(baseName));
end
end

function D = bin_image_matlab(auxImage, maxI, nLevels)
if maxI <= 0
    error('Expected positive max intensity for GLRLM verifier.');
end

slope = double(nLevels) / double(maxI);
intercept = 1.0;
D = zeros(size(auxImage));
for idx = 1:numel(auxImage)
    D(idx) = bin_pixel_matlab(auxImage(idx), slope, intercept, nLevels);
end
end

function y = bin_pixel_matlab(x, slope, intercept, nLevels)
if x == 0
    y = 1;
    return;
end

y = floor(slope * double(x) + intercept);
if y > nLevels
    y = nLevels;
end
if y < 1
    y = 1;
end
end

function [P, levels] = build_glrlm_matrix(D, angleIdx)
levels = unique(D(:));
levels(levels == 0) = [];
levels = sort(levels(:))';

[h, w] = size(D);
visited = false(h, w);
runs = zeros(0, 2);
maxRun = 0;

for row = 1:h
    for col = 1:w
        pi = D(row, col);
        if pi == 0 || visited(row, col)
            continue;
        end

        visited(row, col) = true;
        runLen = 1;
        x = col;
        y = row;

        while true
            [yn, xn] = glrlm_next_pixel(y, x, angleIdx);
            if yn < 1 || yn > h || xn < 1 || xn > w
                break;
            end
            if visited(yn, xn) || D(yn, xn) ~= pi
                break;
            end

            visited(yn, xn) = true;
            runLen = runLen + 1;
            y = yn;
            x = xn;
        end

        runs(end + 1, :) = [pi, runLen]; %#ok<AGROW>
        maxRun = max(maxRun, runLen);
    end
end

P = zeros(numel(levels), maxRun);
for i = 1:size(runs, 1)
    rowIdx = find(levels == runs(i, 1), 1);
    colIdx = runs(i, 2);
    P(rowIdx, colIdx) = P(rowIdx, colIdx) + 1;
end
end

function [yn, xn] = glrlm_next_pixel(y, x, angleIdx)
switch angleIdx
    case 1
        yn = y;
        xn = x + 1;
    case 2
        yn = y + 1;
        xn = x + 1;
    case 3
        yn = y + 1;
        xn = x;
    case 4
        yn = y + 1;
        xn = x - 1;
    otherwise
        error('Unexpected GLRLM angle index: %d', angleIdx);
end
end

function f = compute_gldzm_slice_features(auxImage)
D = auxImage;
maxLevel = max(D(:));
levels = 1:maxLevel;
[h, w] = size(D);
visited = false(h, w);
zones = zeros(0, 3);

for row = 1:h
    for col = 1:w
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
    rowIdx = zones(i, 1);
    colIdx = zones(i, 2);
    P(rowIdx, colIdx) = P(rowIdx, colIdx) + 1;
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
f.GLDZM_ZP = 0.0;
f.GLDZM_GLM = 0.0;
f.GLDZM_GLV = 0.0;
f.GLDZM_ZDM = 0.0;
f.GLDZM_ZDV = 0.0;
f.GLDZM_ZDE = 0.0;

for d = 1:Nd
    m = Md(d);
    f.GLDZM_SDE = f.GLDZM_SDE + m / (d * d);
    f.GLDZM_LDE = f.GLDZM_LDE + (d * d) * m;
    f.GLDZM_ZDNU = f.GLDZM_ZDNU + m * m;
end
f.GLDZM_SDE = f.GLDZM_SDE / Ns;
f.GLDZM_LDE = f.GLDZM_LDE / Ns;
f.GLDZM_ZDNU = f.GLDZM_ZDNU / Ns;
f.GLDZM_ZDNUN = f.GLDZM_ZDNU / Ns;

for g = 1:numel(levels)
    grey = levels(g);
    x = Mx(g);
    f.GLDZM_LGLZE = f.GLDZM_LGLZE + x / (grey * grey);
    f.GLDZM_HGLZE = f.GLDZM_HGLZE + (grey * grey) * x;
    f.GLDZM_GLNU = f.GLDZM_GLNU + x * x;
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
f.GLDZM_ZP = Ns / Nv;

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
[h, w] = size(D);
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
        if r < 1 || r > h || c < 1 || c > w
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

dist2l = dist2l + 1;
dist2r = dist2r + 1;
dist2t = dist2t + 1;
dist2b = dist2b + 1;

retval = min([dist2l, dist2r, dist2t, dist2b]);
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
