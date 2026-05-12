function verify_nyxus_shape2d_morphology_tests(repoRoot)
% Cross-verify Nyxus 2D shape/morphology unit-test expectations in MATLAB.
%
% This script:
% 1. Clones Nyxus if repoRoot is not provided
% 2. Reads the shared synthetic 2D shape fixture from tests/test_data.h
% 3. Reads the asserted constants from tests/test_shape_morphology_2d.h
% 4. Recomputes the active 2D shape/morphology regression surface in MATLAB
% 5. Reports whether each value would pass Nyxus's C++ test tolerances

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = resolve_default_nyxus_repo();
end

repoRoot = char(repoRoot);
testDataPath = fullfile(repoRoot, "tests", "test_data.h");
headerPath = fullfile(repoRoot, "tests", "test_shape_morphology_2d.h");

truth = load_named_constants(headerPath);
roiPixels = load_shape_fixture(testDataPath);
results = compute_shape_features(roiPixels);
features = build_shape_feature_list();

fprintf('\nNyxus 2D shape/morphology verifier\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Fixture source: %s\n', testDataPath);
fprintf('Test constants source: %s\n', headerPath);
fprintf('Configuration: synthetic irregular 2D ROI, PIXELSIZEUM=2, XYRES=1, current Nyxus 2D shape/morphology semantics\n\n');

[checkedCount, passCount] = report_shape_features(features, truth, results);

fprintf('\nSummary\n');
fprintf('2D shape/morphology: %d/%d features pass Nyxus''s C++ test tolerances.\n', passCount, checkedCount);
fprintf('\nNotes:\n');
fprintf('1. This verifier mirrors Nyxus''s current implementation semantics, including the current diameter-equal-area, compactness, ROI-radius, Euler, and convex-hull conventions.\n');
fprintf('2. Contour-driven features use an independent MATLAB implementation of Nyxus''s padded ordered-contour convention, then recompute the feature formulas over that contour.\n');
fprintf('3. This is an independent MATLAB implementation over the same shared fixture, not a call back into Nyxus.\n\n');

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
    fullfile(repoRoot, "tests", "test_shape_morphology_2d.h")
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

function truth = load_named_constants(headerPath)
txt = fileread(headerPath);
tokens = regexp(txt, '\{\s*"([A-Z0-9_]+)"\s*,\s*([-+0-9.eE]+)\s*\}', 'tokens');
if isempty(tokens)
    error('Could not parse named constants in %s', headerPath);
end

truth = struct();
for i = 1:numel(tokens)
    name = tokens{i}{1};
    if ~isfield(truth, name)
        truth.(name) = str2double(tokens{i}{2});
    end
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

function roiPixels = load_shape_fixture(testDataPath)
intensityTriples = load_nyxus_pixel_array(testDataPath, 'shape2d_morphology_intensity');
maskTriples = load_nyxus_pixel_array(testDataPath, 'shape2d_morphology_mask');

roiMask = maskTriples(:, 3) ~= 0;
roiPixels = intensityTriples(roiMask, :);
end

function features = build_shape_feature_list()
features = {
    'AREA_PIXELS_COUNT'
    'AREA_UM2'
    'CENTROID_X'
    'CENTROID_Y'
    'WEIGHTED_CENTROID_X'
    'WEIGHTED_CENTROID_Y'
    'MASS_DISPLACEMENT'
    'COMPACTNESS'
    'BBOX_XMIN'
    'BBOX_YMIN'
    'BBOX_WIDTH'
    'BBOX_HEIGHT'
    'DIAMETER_EQUAL_AREA'
    'EXTENT'
    'ASPECT_RATIO'
    'MAJOR_AXIS_LENGTH'
    'MINOR_AXIS_LENGTH'
    'ELONGATION'
    'ECCENTRICITY'
    'ORIENTATION'
    'ROUNDNESS'
    'PERIMETER'
    'DIAMETER_EQUAL_PERIMETER'
    'EDGE_MEAN_INTENSITY'
    'EDGE_STDDEV_INTENSITY'
    'EDGE_MAX_INTENSITY'
    'EDGE_MIN_INTENSITY'
    'EDGE_INTEGRATED_INTENSITY'
    'CIRCULARITY'
    'CONVEX_HULL_AREA'
    'SOLIDITY'
    'EULER_NUMBER'
    'FRACT_DIM_BOXCOUNT'
    'FRACT_DIM_PERIMETER'
    'DIAMETER_MIN_ENCLOSING_CIRCLE'
    'DIAMETER_CIRCUMSCRIBING_CIRCLE'
    'DIAMETER_INSCRIBING_CIRCLE'
    'ROI_RADIUS_MEAN'
    'ROI_RADIUS_MAX'
    'ROI_RADIUS_MEDIAN'
    'EXTREMA_P1_X'
    'EXTREMA_P1_Y'
    'EXTREMA_P2_X'
    'EXTREMA_P2_Y'
    'EXTREMA_P3_X'
    'EXTREMA_P3_Y'
    'EXTREMA_P4_X'
    'EXTREMA_P4_Y'
    'EXTREMA_P5_X'
    'EXTREMA_P5_Y'
    'EXTREMA_P6_X'
    'EXTREMA_P6_Y'
    'EXTREMA_P7_X'
    'EXTREMA_P7_Y'
    'EXTREMA_P8_X'
    'EXTREMA_P8_Y'
    };
end

function [checkedCount, passCount] = report_shape_features(features, truth, results)
checkedCount = 0;
passCount = 0;
fracTolerance = 1000.0;

fprintf('2D Shape/Morphology\n');
fprintf('%-32s %18s %18s %12s %8s\n', 'Feature', 'Expected', 'MATLAB', 'AbsDiff', 'Pass');
fprintf('%s\n', repmat('-', 1, 94));

for i = 1:numel(features)
    featureName = features{i};
    expected = truth.(featureName);
    actual = results.(featureName);
    pass = agrees_gt_matlab(actual, expected, fracTolerance);

    checkedCount = checkedCount + 1;
    if pass
        passCount = passCount + 1;
    end

    fprintf('%-32s %18.12g %18.12g %12.6g %8s\n', ...
        featureName, expected, actual, abs(actual - expected), pass_fail(pass));
end
end

function results = compute_shape_features(roiPixels)
x = roiPixels(:, 1);
y = roiPixels(:, 2);
inten = roiPixels(:, 3);
n = numel(x);

results = struct();
results.AREA_PIXELS_COUNT = n;
results.AREA_UM2 = n * 4.0;

cenx_sum = sum(x);
ceny_sum = sum(y);
results.COMPACTNESS = sample_std(sqrt((x - cenx_sum).^2 + (y - ceny_sum).^2)) / n;

minX = min(x);
maxX = max(x);
minY = min(y);
maxY = max(y);
width = maxX - minX + 1;
height = maxY - minY + 1;

results.BBOX_XMIN = minX;
results.BBOX_YMIN = minY;
results.BBOX_WIDTH = width;
results.BBOX_HEIGHT = height;

results.CENTROID_X = mean(x);
results.CENTROID_Y = mean(y);

mass = sum(inten);
if mass > 0
    results.WEIGHTED_CENTROID_X = sum((x + 1) .* inten) / mass;
    results.WEIGHTED_CENTROID_Y = sum((y + 1) .* inten) / mass;
else
    results.WEIGHTED_CENTROID_X = 0.0;
    results.WEIGHTED_CENTROID_Y = 0.0;
end

results.MASS_DISPLACEMENT = hypot(results.WEIGHTED_CENTROID_X - results.CENTROID_X, ...
    results.WEIGHTED_CENTROID_Y - results.CENTROID_Y);
results.DIAMETER_EQUAL_AREA = n / pi * 4.0;
results.EXTENT = n / (width * height);
results.ASPECT_RATIO = double(width) / double(height);

[contourPts, contourInt] = nyxus_ordered_contour_2d([x, y], inten);
[results.PERIMETER, results.DIAMETER_EQUAL_PERIMETER] = contour_perimeter(contourPts);
results.EDGE_MEAN_INTENSITY = mean(contourInt);
results.EDGE_STDDEV_INTENSITY = sample_std(contourInt);
results.EDGE_MAX_INTENSITY = max(contourInt);
results.EDGE_MIN_INTENSITY = min(contourInt);
results.EDGE_INTEGRATED_INTENSITY = sum(contourInt);

[results.CONVEX_HULL_AREA, results.SOLIDITY, results.CIRCULARITY] = convex_hull_features(x, y, n, results.PERIMETER);

[results.MAJOR_AXIS_LENGTH, results.MINOR_AXIS_LENGTH, results.ECCENTRICITY, ...
    results.ELONGATION, results.ORIENTATION, results.ROUNDNESS] = ellipse_features(x, y, n, results.CENTROID_X, results.CENTROID_Y);

extrema = extrema_features(x, y, minX, maxX, minY, maxY);
for i = 1:8
    results.(sprintf('EXTREMA_P%d_X', i)) = extrema(i, 1);
    results.(sprintf('EXTREMA_P%d_Y', i)) = extrema(i, 2);
end

results.ROI_RADIUS_MEAN = mean(min_sqdist_to_cloud([x y], contourPts));
roiMinSqDists = min_sqdist_to_cloud([x y], contourPts);
results.ROI_RADIUS_MAX = max(roiMinSqDists);
results.ROI_RADIUS_MEDIAN = median(roiMinSqDists);

results.EULER_NUMBER = calculate_euler_from_pixels(x, y, minX, maxX, minY, maxY);
results.FRACT_DIM_BOXCOUNT = fract_dim_boxcount(x, y, inten, minX, maxX, minY, maxY);
results.FRACT_DIM_PERIMETER = fract_dim_perimeter(contourPts);

[results.DIAMETER_MIN_ENCLOSING_CIRCLE, results.DIAMETER_INSCRIBING_CIRCLE, results.DIAMETER_CIRCUMSCRIBING_CIRCLE] = ...
    circle_features(contourPts, results.CENTROID_X, results.CENTROID_Y);
end

function [perimeter, diameterEqPerim] = contour_perimeter(contourPts)
prevPts = contourPts([end, 1:end-1], :);
perimeter = sum(sqrt(sum((contourPts - prevPts).^2, 2)));
diameterEqPerim = perimeter / pi;
end

function [hullArea, solidity, circularity] = convex_hull_features(x, y, roiArea, perimeter)
hull = nyxus_convex_hull([x, y]);
hullArea = polygon_area(hull);
solidity = roiArea / hullArea;
circularity = sqrt(4.0 * pi * roiArea / (perimeter * perimeter));
end

function hull = nyxus_convex_hull(points)
points = sortrows(points, [1, 2]);
n = size(points, 1);
if n < 2
    hull = zeros(0, 2);
    return;
end

upper = points(1:2, :);
for i = 3:n
    while size(upper, 1) > 1 && ~right_turn(upper(end - 1, :), upper(end, :), points(i, :))
        upper(end, :) = [];
    end
    upper(end + 1, :) = points(i, :); %#ok<AGROW>
end

lower = points(n:-1:n-1, :);
for i = 3:n
    p = points(n - i + 1, :);
    while size(lower, 1) > 1 && ~right_turn(lower(end - 1, :), lower(end, :), p)
        lower(end, :) = [];
    end
    lower(end + 1, :) = p; %#ok<AGROW>
end

hull = upper;
for i = 1:size(lower, 1)
    if ~any(hull(:, 1) == lower(i, 1) & hull(:, 2) == lower(i, 2))
        hull(end + 1, :) = lower(i, :); %#ok<AGROW>
    end
end
end

function tf = right_turn(p1, p2, p3)
tf = ((p3(1) - p1(1)) * (p2(2) - p1(2)) - ...
      (p3(2) - p1(2)) * (p2(1) - p1(1))) > 0;
end

function area = polygon_area(vertices)
if isempty(vertices)
    area = 0.0;
    return;
end

area = 0.0;
n = size(vertices, 1);
for i = 1:(n - 1)
    p1 = vertices(i, :);
    p2 = vertices(i + 1, :);
    area = area + p1(1) * p2(2) - p1(2) * p2(1);
end
p1 = vertices(1, :);
p2 = vertices(n, :);
area = area + p1(1) * p2(2) - p1(2) * p2(1);
area = abs(area) / 2.0;
end

function [majorAxisLength, minorAxisLength, eccentricity, elongation, orientation, roundness] = ellipse_features(x, y, area, centroidX, centroidY)
diffX = centroidX - x;
diffY = centroidY - y;

uxx = mean(diffX .^ 2) + 1.0 / 12.0;
uyy = mean(diffY .^ 2) + 1.0 / 12.0;
uxy = mean(diffX .* diffY);

common = sqrt((uxx - uyy) ^ 2 + 4.0 * uxy * uxy);
majorAxisLength = 2.0 * sqrt(2.0) * sqrt(uxx + uyy + common);
minorAxisLength = 2.0 * sqrt(2.0) * sqrt(uxx + uyy - common);
eccentricity = sqrt(1.0 - minorAxisLength * minorAxisLength / (majorAxisLength * majorAxisLength));
elongation = minorAxisLength / majorAxisLength;
roundness = (4.0 * area) / (pi * majorAxisLength * majorAxisLength);

if uyy > uxx
    num = uyy - uxx + sqrt((uyy - uxx) ^ 2 + 4 * uxy * uxy);
    den = 2 * uxy;
else
    num = 2 * uxy;
    den = uxx - uyy + sqrt((uxx - uyy) ^ 2 + 4 * uxy * uxy);
end

if uxy == 0.0
    if uxx >= uyy
        orientation = 0.0;
    else
        orientation = 90.0;
    end
else
    orientation = 180.0 / pi * atan(num / den);
end
end

function extrema = extrema_features(x, y, leftMost, rightMost, topMost, lowest)
topX = x(y == topMost);
bottomX = x(y == lowest);
leftY = y(x == leftMost);
rightY = y(x == rightMost);

extrema = [
    min(topX), topMost
    max(topX), topMost
    rightMost, min(rightY)
    rightMost, max(rightY)
    max(bottomX), lowest
    min(bottomX), lowest
    leftMost, max(leftY)
    leftMost, min(leftY)
    ];
end

function dists = min_sqdist_to_cloud(points, cloud)
dists = zeros(size(points, 1), 1);
for i = 1:size(points, 1)
    dists(i) = nyxus_min_sqdist_ordered_2d(points(i, :), cloud);
end
end

function eulerNumber = calculate_euler_from_pixels(x, y, minX, maxX, minY, maxY)
height = maxY - minY + 1;
width = maxX - minX + 1;
I = zeros(height, width);
for i = 1:numel(x)
    I(y(i) - minY + 1, x(i) - minX + 1) = 1;
end

Px = uint8([ ...
    8, 4, 2, 1, ...
    7, 11, 13, 14, ...
    9, 6]);

C1 = 0;
C3 = 0;
Cd = 0;
for row = 2:height
    for col = 2:width
        Imq = uint8(0);
        if I(row - 1, col - 1) > 0, Imq = bitor(Imq, uint8(8)); end
        if I(row - 1, col) > 0, Imq = bitor(Imq, uint8(4)); end
        if I(row, col - 1) > 0, Imq = bitor(Imq, uint8(2)); end
        if I(row, col) > 0, Imq = bitor(Imq, uint8(1)); end

        idx = find(Px == Imq, 1);
        if isempty(idx)
            continue;
        end
        if idx <= 4
            C1 = C1 + 1;
        elseif idx <= 8
            C3 = C3 + 1;
        elseif idx == 9 && idx == 10
            % Mirrors current Nyxus behavior: diagonal Cd patterns are not
            % counted. A separate diagonal-pattern regression should cover
            % any future Euler connectivity bug fix.
            Cd = Cd + 1;
        end
    end
end

eulerNumber = (C1 - C3 - (2 * Cd)) / 4;
end

function fd = fract_dim_boxcount(x, y, inten, minX, maxX, minY, maxY)
width = maxX - minX + 1;
height = maxY - minY + 1;
bigSide = max(width, height);
paddedSide = closest_pow2(bigSide);

padOffsetX = floor((paddedSide - width) / 2);
padOffsetY = floor((paddedSide - height) / 2);
I = zeros(paddedSide, paddedSide);
for i = 1:numel(x)
    px = x(i) - minX + padOffsetX + 1;
    py = y(i) - minY + padOffsetY + 1;
    I(py, px) = inten(i) + 1;
end

coverage = zeros(0, 2);
s = paddedSide;
while s > 1
    cnt = 0;
    nTiles = paddedSide / s;
    for r = 1:nTiles
        for c = 1:nTiles
            rows = (r - 1) * s + 1:r * s;
            cols = (c - 1) * s + 1:c * s;
            if any(I(rows, cols) > 0, 'all')
                cnt = cnt + 1;
            end
        end
    end
    coverage(end + 1, :) = [s, cnt]; %#ok<AGROW>
    s = floor(s / 2);
end

fd = -calc_lyapunov_slope(coverage);
end

function fd = fract_dim_perimeter(contourPts)
coverage = zeros(0, 2);
conLen = size(contourPts, 1);
s = floor(conLen / 4);
while s > 0
    p = 0.0;
    for idx = (s + 1):s:conLen
        px1 = contourPts(idx - s, :);
        px2 = contourPts(idx, :);
        p = p + sqrt(sum((px1 - px2) .^ 2));
    end
    tail = mod(conLen, s);
    if tail > 0
        px1 = contourPts(conLen - tail + 1, :);
        px2 = contourPts(1, :);
        p = p + sum((px1 - px2) .^ 2);
    end
    coverage(end + 1, :) = [s, fix(p)]; %#ok<AGROW>
    s = floor(s / 2);
end

fd = calc_lyapunov_slope(coverage);
end

function slope = calc_lyapunov_slope(coverage)
if size(coverage, 1) < 2
    slope = 0.0;
    return;
end

X = log(coverage(:, 1));
Y = log(coverage(:, 2));
Dn = diff(Y);
Dr = diff(X);
lambda = Dn ./ Dr;

n = numel(lambda);
xs = (0:n-1)';
sum_x = sum(xs);
sum_y = sum(lambda);
sum_xy = sum(xs .* lambda);
sum_x2 = sum(xs .* xs);
slope = (sum_xy * n - sum_x * sum_y) / (sum_x2 * n - sum_x * sum_x);
end

function [diamMinEnclosing, diamInscribing, diamCircumscribing] = circle_features(contourPts, centroidX, centroidY)
[center, radius] = min_enclosing_circle(contourPts);
diamMinEnclosing = 2 * radius;

xCentroid2 = centroidX - 1;
yCentroid2 = centroidY - 1;
dists = sqrt((contourPts(:, 1) - xCentroid2) .^ 2 + (contourPts(:, 2) - yCentroid2) .^ 2);
diamInscribing = 2 * min(dists);
diamCircumscribing = 2 * max(dists);
end

function [bestCenter, bestRadius] = min_enclosing_circle(points)
pts = double(points);
n = size(pts, 1);
if n == 0
    bestCenter = [0, 0];
    bestRadius = 0;
    return;
elseif n == 1
    bestCenter = pts(1, :);
    bestRadius = eps;
    return;
end

bestRadius = inf;
bestCenter = [0, 0];

for i = 1:n
    for j = i+1:n
        center = (pts(i, :) + pts(j, :)) / 2.0;
        radius = norm(pts(i, :) - center);
        if all(vecnorm(pts - center, 2, 2) <= radius + 1e-9) && radius < bestRadius
            bestRadius = radius;
            bestCenter = center;
        end
    end
end

for i = 1:n
    for j = i+1:n
        for k = j+1:n
            [center, radius, valid] = circumcircle(pts(i, :), pts(j, :), pts(k, :));
            if valid && all(vecnorm(pts - center, 2, 2) <= radius + 1e-9) && radius < bestRadius
                bestRadius = radius;
                bestCenter = center;
            end
        end
    end
end
end

function [center, radius, valid] = circumcircle(p1, p2, p3)
d = 2 * (p1(1) * (p2(2) - p3(2)) + p2(1) * (p3(2) - p1(2)) + p3(1) * (p1(2) - p2(2)));
if abs(d) <= eps
    center = [0, 0];
    radius = 0;
    valid = false;
    return;
end

ux = ((p1(1)^2 + p1(2)^2) * (p2(2) - p3(2)) + ...
      (p2(1)^2 + p2(2)^2) * (p3(2) - p1(2)) + ...
      (p3(1)^2 + p3(2)^2) * (p1(2) - p2(2))) / d;
uy = ((p1(1)^2 + p1(2)^2) * (p3(1) - p2(1)) + ...
      (p2(1)^2 + p2(2)^2) * (p1(1) - p3(1)) + ...
      (p3(1)^2 + p3(2)^2) * (p2(1) - p1(1))) / d;
center = [ux, uy];
radius = norm(p1 - center);
valid = true;
end

function p = closest_pow2(n)
p = 1;
while p < n
    p = p * 2;
end
end

function s = sample_std(v)
if numel(v) > 2
    s = std(v, 0);
else
    s = 0.0;
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
