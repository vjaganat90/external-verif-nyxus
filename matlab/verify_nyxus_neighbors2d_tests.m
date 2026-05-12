function verify_nyxus_neighbors2d_tests(repoRoot)
% Cross-verify Nyxus 2D neighborhood/touching unit-test expectations in MATLAB.
%
% This script:
% 1. Clones Nyxus if repoRoot is not provided
% 2. Reads the shared synthetic multi-label scene from tests/test_data.h
% 3. Reads the asserted constants from tests/test_neighbors_2d.h
% 4. Recomputes the active 2D neighborhood/touching regression surface in MATLAB
% 5. Reports whether each value would pass Nyxus's C++ test tolerances

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
headerPath = fullfile(repoRoot, "tests", "test_neighbors_2d.h");

truth = load_neighbor_truth(headerPath);
rois = load_neighbor_fixture(testDataPath);
results = compute_neighbor_features(rois);
featureList = build_neighbor_feature_list();

fprintf('\nNyxus 2D neighborhood/touching verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Fixture source: %s\n', testDataPath);
fprintf('Test constants source: %s\n', headerPath);
fprintf('Configuration: synthetic 5-ROI 2D scene, PIXELDISTANCE=1, current Nyxus 2D neighborhood/touching semantics\n\n');

[checkedCount, passCount] = report_neighbor_features(featureList, truth, results);

fprintf('\nSummary\n');
fprintf('2D neighborhood/touching: %d/%d values pass Nyxus''s C++ test tolerances.\n', passCount, checkedCount);
fprintf('\nNotes:\n');
fprintf('1. This verifier mirrors the current branch semantics of the neighborhood reducer, including radius-aware broad-phase pairing, Euclidean centroid distances, and direction angles in [0, 360).\n');
fprintf('2. Touching percentages are measured as the percentage of contour pixels within the configured neighborhood radius of another ROI.\n');
fprintf('3. This is an independent MATLAB implementation over the same shared scene fixture, not a call back into Nyxus.\n\n');

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
    fullfile(repoRoot, "tests", "test_neighbors_2d.h")
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
cmd = sprintf('git -C \"%s\" fetch origin %s --depth 1 && git -C \"%s\" checkout %s && git -C \"%s\" reset --hard origin/%s', ...
    repoRoot, gitBranch, repoRoot, gitBranch, repoRoot, gitBranch);
[status, cmdout] = system(cmd);
if status ~= 0
    error('Failed to refresh Nyxus clone.\nOutput:\n%s', cmdout);
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

function truth = load_neighbor_truth(headerPath)
txt = fileread(headerPath);
labelBlocks = regexp(txt, '\{(\d+),\s*\{([\s\S]*?)\}\s*\}', 'tokens');
if isempty(labelBlocks)
    error('Could not parse label-indexed neighbor truth in %s', headerPath);
end

truth = struct();
for i = 1:numel(labelBlocks)
    label = labelBlocks{i}{1};
    entries = regexp(labelBlocks{i}{2}, '\{\s*"([A-Z0-9_]+)"\s*,\s*([-+0-9.eE]+)\s*\}', 'tokens');
    row = struct();
    for j = 1:numel(entries)
        row.(entries{j}{1}) = str2double(entries{j}{2});
    end
    truth.(sprintf('L%s', label)) = row;
end
end

function rois = load_neighbor_fixture(testDataPath)
triples = load_nyxus_pixel_array(testDataPath, 'neighborhood2d_scene_labels');
labels = unique(triples(:, 3));
labels = sort(labels(:))';

rois = struct();
for label = labels
    idx = triples(:, 3) == label;
    pts = triples(idx, 1:2);
    minX = min(pts(:, 1));
    maxX = max(pts(:, 1));
    minY = min(pts(:, 2));
    maxY = max(pts(:, 2));

    contourPts = nyxus_ordered_contour_2d(pts, ones(size(pts, 1), 1));

    key = sprintf('L%d', label);
    rois.(key) = struct( ...
        'label', label, ...
        'points', pts, ...
        'bbox', [minX, maxX, minY, maxY], ...
        'centroid', [mean(pts(:, 1)), mean(pts(:, 2))], ...
        'contour', contourPts, ...
        'neighborLabels', zeros(0, 1), ...
        'touchingCount', 0.0, ...
        'results', struct());
end
end

function featureList = build_neighbor_feature_list()
featureList = {
    'NUM_NEIGHBORS'
    'PERCENT_TOUCHING'
    'CLOSEST_NEIGHBOR1_DIST'
    'CLOSEST_NEIGHBOR1_ANG'
    'CLOSEST_NEIGHBOR2_DIST'
    'CLOSEST_NEIGHBOR2_ANG'
    'ANG_BW_NEIGHBORS_MEAN'
    'ANG_BW_NEIGHBORS_STDDEV'
    'ANG_BW_NEIGHBORS_MODE'
    };
end

function [checkedCount, passCount] = report_neighbor_features(featureList, truth, results)
checkedCount = 0;
passCount = 0;
fracTolerance = 1000.0;
labels = {'L1', 'L2', 'L3', 'L4', 'L5'};

fprintf('2D Neighborhood/Touching\n');
fprintf('%-10s %-26s %18s %18s %12s %8s\n', 'ROI', 'Feature', 'Expected', 'MATLAB', 'AbsDiff', 'Pass');
fprintf('%s\n', repmat('-', 1, 102));

for i = 1:numel(labels)
    labelKey = labels{i};
    for j = 1:numel(featureList)
        featureName = featureList{j};
        expected = truth.(labelKey).(featureName);
        actual = results.(labelKey).(featureName);
        pass = agrees_gt_matlab(actual, expected, fracTolerance);

        checkedCount = checkedCount + 1;
        if pass
            passCount = passCount + 1;
        end

        fprintf('%-10s %-26s %18.12g %18.12g %12.6g %8s\n', ...
            labelKey, featureName, expected, actual, abs(actual - expected), pass_fail(pass));
    end
end
end

function results = compute_neighbor_features(rois)
radius = 1;
radius2 = radius * radius;
labels = sort(cellfun(@(k) rois.(k).label, fieldnames(rois)));

for i = 1:numel(labels)
    l1 = labels(i);
    r1 = rois.(sprintf('L%d', l1));
    for j = i+1:numel(labels)
        l2 = labels(j);
        r2 = rois.(sprintf('L%d', l2));
        if aabb_no_overlap(r1.bbox, r2.bbox, radius)
            continue;
        end

        [minD12, touch1] = contour_relation(r1.contour, r2.contour, radius2);
        [~, touch2] = contour_relation(r2.contour, r1.contour, radius2);
        if minD12 > radius2
            continue;
        end

        r1.touchingCount = r1.touchingCount + touch1;
        r2.touchingCount = r2.touchingCount + touch2;
        r1.neighborLabels(end + 1, 1) = l2; %#ok<AGROW>
        r2.neighborLabels(end + 1, 1) = l1; %#ok<AGROW>
        rois.(sprintf('L%d', l1)) = r1;
        rois.(sprintf('L%d', l2)) = r2;
    end
end

for i = 1:numel(labels)
    label = labels(i);
    key = sprintf('L%d', label);
    roi = rois.(key);
    nNeighbors = numel(roi.neighborLabels);
    contourLen = size(roi.contour, 1);

    roi.results.NUM_NEIGHBORS = nNeighbors;
    if contourLen > 0
        roi.results.PERCENT_TOUCHING = 100.0 * roi.touchingCount / contourLen;
    else
        roi.results.PERCENT_TOUCHING = 0.0;
    end

    roi.results.CLOSEST_NEIGHBOR1_DIST = 0.0;
    roi.results.CLOSEST_NEIGHBOR1_ANG = 0.0;
    roi.results.CLOSEST_NEIGHBOR2_DIST = 0.0;
    roi.results.CLOSEST_NEIGHBOR2_ANG = 0.0;
    roi.results.ANG_BW_NEIGHBORS_MEAN = 0.0;
    roi.results.ANG_BW_NEIGHBORS_STDDEV = 0.0;
    roi.results.ANG_BW_NEIGHBORS_MODE = 0.0;

    if nNeighbors > 0
        cen = roi.centroid;
        dists = zeros(nNeighbors, 1);
        dirs = zeros(nNeighbors, 1);
        for k = 1:nNeighbors
            neigh = rois.(sprintf('L%d', roi.neighborLabels(k)));
            neighCen = neigh.centroid;
            dists(k) = hypot(cen(1) - neighCen(1), cen(2) - neighCen(2));
            dirs(k) = direction_angle_deg(cen(1), cen(2), neighCen(1), neighCen(2));
        end

        [~, idx1] = min(dists);
        roi.results.CLOSEST_NEIGHBOR1_DIST = dists(idx1);
        roi.results.CLOSEST_NEIGHBOR1_ANG = dirs(idx1);

        if nNeighbors > 1
            d2 = dists;
            d2(idx1) = inf;
            [~, idx2] = min(d2);
            roi.results.CLOSEST_NEIGHBOR2_DIST = dists(idx2);
            roi.results.CLOSEST_NEIGHBOR2_ANG = dirs(idx2);
        end

        roi.results.ANG_BW_NEIGHBORS_MEAN = mean(dirs);
        roi.results.ANG_BW_NEIGHBORS_STDDEV = sample_std(dirs);
        rounded = sort(round(dirs));
        roi.results.ANG_BW_NEIGHBORS_MODE = nyx_mode(rounded);
    end

    rois.(key) = roi;
end

results = struct();
for i = 1:numel(labels)
    key = sprintf('L%d', labels(i));
    results.(key) = rois.(key).results;
end
end

function [minD, touchingCount] = contour_relation(contourA, contourB, radius2)
touchingCount = 0;
minD = inf;
for i = 1:size(contourA, 1)
    localMin = nyxus_min_sqdist_ordered_2d(contourA(i, :), contourB);
    minD = min(minD, localMin);
    if localMin <= radius2
        touchingCount = touchingCount + 1;
    end
end
end

function tf = aabb_no_overlap(b1, b2, radius)
tf = (b2(1) - radius > b1(2) + radius) || ...
     (b2(2) + radius < b1(1) - radius) || ...
     (b2(3) - radius > b1(4) + radius) || ...
     (b2(4) + radius < b1(3) - radius);
end

function ang = direction_angle_deg(x1, y1, x2, y2)
ang = atan2d(y2 - y1, x2 - x1);
if ang < 0.0
    ang = ang + 360.0;
end
end

function s = sample_std(v)
if numel(v) > 2
    s = std(v, 0);
else
    s = 0.0;
end
end

function m = nyx_mode(v)
prev = v(end);
m = 0;
maxcount = 0;
currcount = 0;
for i = 1:numel(v)
    n = v(i);
    if n == prev
        currcount = currcount + 1;
        if currcount > maxcount
            maxcount = currcount;
            m = n;
        end
    else
        currcount = 1;
    end
    prev = n;
end
end

function pass = agrees_gt_matlab(actual, expected, fracTolerance)
tolerance = abs(expected) / fracTolerance;
pass = abs(actual - expected) <= tolerance;
end

function s = pass_fail(pass)
if pass
    s = 'PASS';
else
    s = 'FAIL';
end
end
