function verify_nyxus_gabor_tests(repoRoot)
% Cross-verify Nyxus 2D Gabor unit-test expectations in MATLAB.
%
% This script:
% 1. Clones the Nyxus repo under test if repoRoot is not provided
% 2. Reads the DSB2018 ROI fixtures from tests/test_dsb2018_data.h
% 3. Reads the asserted Gabor vector constants from tests/test_gabor_truth.h
% 4. Recomputes Nyxus's CPU Gabor score in MATLAB
% 5. Reports whether each value would pass Nyxus's C++ test tolerances

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
dataPath = fullfile(repoRoot, "tests", "test_dsb2018_data.h");
truthPath = fullfile(repoRoot, "tests", "test_gabor_truth.h");

images = load_dsb_images(dataPath);
truth = load_gabor_truth(truthPath);

fprintf('\nNyxus 2D Gabor verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('ROI fixture source: %s\n', dataPath);
fprintf('Test constants source: %s\n', truthPath);
fprintf('ROIs in fixture: %d\n', numel(images));
fprintf('Configuration: CPU Gabor defaults from src/nyx/features/gabor.cpp, including current f0/theta pair unpacking semantics\n\n');

[checkedCount, passCount] = report_gabor(images, truth);

fprintf('\nSummary\n');
fprintf('Gabor: %d/%d vector values pass Nyxus''s C++ test tolerances.\n', passCount, checkedCount);
fprintf('\nNotes:\n');
fprintf('1. The MATLAB code independently reconstructs the complex Gabor kernels, full convolution, central crop, unsigned-integer energy image, and threshold score used by Nyxus CPU tests.\n');
fprintf('2. The active C++ implementation stores pairs as documented theta/frequency values but currently unpacks them as f0/theta during calculation; this verifier mirrors the executable behavior under test.\n');
fprintf('3. This is an independent MATLAB implementation over the same DSB ROI fixtures, not a call back into Nyxus.\n\n');

if passCount ~= checkedCount
    error('Gabor cross-verification failed.');
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
    fullfile(repoRoot, "tests", "test_dsb2018_data.h")
    fullfile(repoRoot, "tests", "test_gabor_truth.h")
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

function images = load_dsb_images(dataPath)
txt = fileread(dataPath);
block = regexp(txt, 'const\s+static\s+std::vector<ImageData>\s+dsb_data\s*\{([\s\S]*?)\n\s*\};', 'tokens', 'once');
if isempty(block)
    error('Could not find dsb_data in %s', dataPath);
end

tokens = regexp(block{1}, '\{\s*(\d+)\s*,\s*(\d+)\s*,\s*\{([^{}]*)\}\s*\}', 'tokens');
if isempty(tokens)
    error('Could not parse dsb_data in %s', dataPath);
end

images = cell(numel(tokens), 1);
for i = 1:numel(tokens)
    width = str2double(tokens{i}{1});
    height = str2double(tokens{i}{2});
    vals = regexp(tokens{i}{3}, '-?\d+', 'match');
    pixels = zeros(numel(vals), 1);
    for j = 1:numel(vals)
        pixels(j) = str2double(vals{j});
    end
    if numel(pixels) ~= width * height
        error('DSB ROI %d has %d pixels, expected %d.', i, numel(pixels), width * height);
    end
    images{i} = reshape(pixels, width, height)';
end
end

function truth = load_gabor_truth(truthPath)
txt = fileread(truthPath);
block = regexp(txt, 'gabor_truth\s*=\s*\{([\s\S]*?)\n\s*\};', 'tokens', 'once');
if isempty(block)
    error('Could not find gabor_truth in %s', truthPath);
end

rows = regexp(block{1}, '\{([^{}]*)\}', 'tokens');
if isempty(rows)
    error('Could not parse gabor_truth in %s', truthPath);
end

truth = cell(numel(rows), 1);
for i = 1:numel(rows)
    vals = regexp(rows{i}{1}, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
    truth{i} = zeros(numel(vals), 1);
    for j = 1:numel(vals)
        truth{i}(j) = str2double(vals{j});
    end
end
end

function [checkedCount, passCount] = report_gabor(images, truth)
if numel(images) ~= numel(truth)
    error('Fixture count %d does not match truth count %d.', numel(images), numel(truth));
end

checkedCount = 0;
passCount = 0;
fracTolerance = 1000.0;

fprintf('Gabor\n');
fprintf('%-8s %-8s %18s %18s %12s %8s\n', 'ROI', 'Filter', 'Expected', 'MATLAB', 'AbsDiff', 'Pass');
fprintf('%s\n', repmat('-', 1, 88));

for roiIdx = 1:numel(images)
    actual = compute_gabor_features(images{roiIdx});
    expected = truth{roiIdx};
    if numel(actual) ~= numel(expected)
        error('ROI %d produced %d values, expected %d.', roiIdx, numel(actual), numel(expected));
    end
    for filterIdx = 1:numel(expected)
        pass = agrees_gt_matlab(actual(filterIdx), expected(filterIdx), fracTolerance);
        checkedCount = checkedCount + 1;
        if pass
            passCount = passCount + 1;
        end
        fprintf('%-8d %-8d %18.12g %18.12g %12.6g %8s\n', ...
            roiIdx, filterIdx, expected(filterIdx), actual(filterIdx), abs(actual(filterIdx) - expected(filterIdx)), pass_fail(pass));
    end
end
end

function values = compute_gabor_features(image)
kernelSide = 16;
sig2lam = 0.8;
gamma = 0.1;
f0LP = 0.1;
grayThreshold = 0.025;

baseline = gabor_energy(image, f0LP, sig2lam, gamma, pi / 2.0, kernelSide);
maxVal = max(baseline(:));
cmpVal = min(baseline(:));

if maxVal == cmpVal
    values = nan(4, 1);
    return;
end

baselineScore = nnz(baseline > cmpVal);
if baselineScore == 0
    error('Gabor baseline score is zero.');
end

% These are the exact pairs in GaborFeature::f0_theta_pairs, unpacked by
% calculate() as f0=first, theta=second.
pairs = [
    0.0,        4.0
    pi / 4.0,   16.0
    pi / 2.0,   32.0
    3.0*pi/4.0, 64.0
    ];

values = zeros(size(pairs, 1), 1);
for i = 1:size(pairs, 1)
    response = gabor_energy(image, pairs(i, 1), sig2lam, gamma, pairs(i, 2), kernelSide);
    afterGaborScore = nnz((double(response) ./ double(maxVal)) > grayThreshold);
    values(i) = double(afterGaborScore) / double(baselineScore);
end
end

function out = gabor_energy(image, f0, sig2lam, gamma, theta, kernelSide)
G = gabor_kernel(f0, sig2lam, gamma, theta, 0.0, kernelSide);
C = conv_dud_matlab(image, G);
startIdx = ceil(double(kernelSide) / 2.0) + 1;
crop = C(startIdx:startIdx + size(image, 1) - 1, startIdx:startIdx + size(image, 2) - 1);
out = floor(abs(crop));
end

function G = gabor_kernel(f0, sig2lam, gamma, theta, fi, kernelSide)
if f0 == 0
    lambda = Inf;
else
    lambda = 2.0 * pi / f0;
end
sig = sig2lam * lambda;

if mod(kernelSide, 2) > 0
    tx = -((kernelSide - 1) / 2):((kernelSide - 1) / 2);
    ty = tx;
else
    tx = -(kernelSide / 2):(kernelSide / 2 - 1);
    ty = tx;
end

G = complex(zeros(kernelSide, kernelSide));
sumAbs = 0.0;
cosTheta = cos(theta);
sinTheta = sin(theta);

for y = 1:kernelSide
    for x = 1:kernelSide
        xte = tx(x) * cosTheta + ty(y) * sinTheta;
        yte = ty(y) * cosTheta - tx(x) * sinTheta;
        rte = xte * xte + gamma * gamma * yte * yte;
        if isinf(sig)
            ge = 1.0;
        else
            ge = exp(-rte / (2.0 * sig * sig));
        end
        argm = xte * f0 + fi;
        G(y, x) = complex(ge * cos(argm), ge * sin(argm));
        sumAbs = sumAbs + abs(G(y, x));
    end
end

G = G ./ sumAbs;
end

function C = conv_dud_matlab(A, B)
% Mirror GaborFeature::conv_dud, including full convolution dimensions.
[ma, na] = size(A);
[mb, nb] = size(B);
C = complex(zeros(ma + mb - 1, na + nb - 1));

for j = 1:mb
    for i = 1:nb
        C(j:j + ma - 1, i:i + na - 1) = C(j:j + ma - 1, i:i + na - 1) + double(A) .* B(j, i);
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
