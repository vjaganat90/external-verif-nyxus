function results = verify_shape2d_builtin_independent(repoRoot, outDir)
% Compare Nyxus 2D shape/morphology expectations to MATLAB Image Processing
% Toolbox built-ins where comparable definitions exist.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = '/Users/jaganathv2/Code/Codex/nyxus';
end
if nargin < 2 || strlength(string(outDir)) == 0
    outDir = fileparts(mfilename('fullpath'));
end

dataPath = fullfile(repoRoot, 'tests', 'test_data.h');
headerPath = fullfile(repoRoot, 'tests', 'test_shape_morphology_2d.h');

truth = load_named_constants(headerPath);
roi = load_shape_fixture(dataPath);
[bw, inten, minX, minY] = fixture_to_images(roi);

stats = regionprops(bw, inten, ...
    'Area', 'Centroid', 'WeightedCentroid', 'BoundingBox', ...
    'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', ...
    'Orientation', 'Extent', 'ConvexArea', 'Solidity', 'EulerNumber');

if isempty(stats)
    error('No object found in shape fixture.');
end
s = stats(1);

rows = {};
rows = add_builtin(rows, truth, 'AREA_PIXELS_COUNT', s.Area, 'regionprops Area', 'Binary object area in pixels');
rows = add_builtin(rows, truth, 'AREA_UM2', s.Area * 4.0, 'regionprops Area * pixel area', 'Uses Nyxus test PIXELSIZEUM=2');
rows = add_builtin(rows, truth, 'CENTROID_X', s.Centroid(1) + minX - 1, 'regionprops Centroid', 'Converted from MATLAB 1-based column coordinate to Nyxus x');
rows = add_builtin(rows, truth, 'CENTROID_Y', s.Centroid(2) + minY - 1, 'regionprops Centroid', 'Converted from MATLAB 1-based row coordinate to Nyxus y');
rows = add_builtin(rows, truth, 'WEIGHTED_CENTROID_X', s.WeightedCentroid(1) + minX, 'regionprops WeightedCentroid', 'Nyxus weighted centroid is 1-based in current 2D tests');
rows = add_builtin(rows, truth, 'WEIGHTED_CENTROID_Y', s.WeightedCentroid(2) + minY, 'regionprops WeightedCentroid', 'Nyxus weighted centroid is 1-based in current 2D tests');
rows = add_builtin(rows, truth, 'BBOX_XMIN', floor(s.BoundingBox(1) + minX - 0.5), 'regionprops BoundingBox', 'Converted from MATLAB half-pixel bounding box');
rows = add_builtin(rows, truth, 'BBOX_YMIN', floor(s.BoundingBox(2) + minY - 0.5), 'regionprops BoundingBox', 'Converted from MATLAB half-pixel bounding box');
rows = add_builtin(rows, truth, 'BBOX_WIDTH', s.BoundingBox(3), 'regionprops BoundingBox', 'MATLAB bounding box width');
rows = add_builtin(rows, truth, 'BBOX_HEIGHT', s.BoundingBox(4), 'regionprops BoundingBox', 'MATLAB bounding box height');
rows = add_builtin(rows, truth, 'EXTENT', s.Extent, 'regionprops Extent', 'Area divided by bounding-box area');
rows = add_builtin(rows, truth, 'ASPECT_RATIO', s.BoundingBox(3) / s.BoundingBox(4), 'regionprops BoundingBox composition', 'Width divided by height');
rows = add_builtin(rows, truth, 'MAJOR_AXIS_LENGTH', s.MajorAxisLength, 'regionprops MajorAxisLength', 'Equivalent ellipse major axis');
rows = add_builtin(rows, truth, 'MINOR_AXIS_LENGTH', s.MinorAxisLength, 'regionprops MinorAxisLength', 'Equivalent ellipse minor axis');
rows = add_builtin(rows, truth, 'ECCENTRICITY', s.Eccentricity, 'regionprops Eccentricity', 'Equivalent ellipse eccentricity');
rows = add_builtin(rows, truth, 'ELONGATION', s.MinorAxisLength / s.MajorAxisLength, 'regionprops axis composition', 'Minor axis divided by major axis');
rows = add_builtin(rows, truth, 'ORIENTATION', s.Orientation, 'regionprops Orientation', 'MATLAB uses opposite sign convention from current Nyxus tests');
rows = add_builtin(rows, truth, 'EULER_NUMBER', bweuler(bw, 8), 'bweuler', 'MATLAB binary Euler number; may differ by object/background convention');

pts = [roi(:, 1), roi(:, 2)];
if size(unique(pts, 'rows'), 1) >= 3
    k = convhull(pts(:, 1), pts(:, 2));
    hullArea = polyarea(pts(k, 1), pts(k, 2));
    rows = add_builtin(rows, truth, 'CONVEX_HULL_AREA', hullArea, 'convhull + polyarea', 'Point-cloud convex hull over ROI pixel centers');
    rows = add_builtin(rows, truth, 'SOLIDITY', s.Area / hullArea, 'regionprops Area + convhull/polyarea', 'Nyxus point-cloud solidity, not regionprops Solidity');
end

% Keep MATLAB regionprops convex metrics visible as explicit convention checks.
rows = add_diagnostic(rows, 'REGIONPROPS_CONVEX_AREA', s.ConvexArea, 'regionprops ConvexArea', 'MATLAB filled-pixel convex area; not the same definition as Nyxus point-cloud hull area');
rows = add_diagnostic(rows, 'REGIONPROPS_SOLIDITY', s.Solidity, 'regionprops Solidity', 'MATLAB Area/ConvexArea; not the same definition as Nyxus point-cloud solidity');

unsupported = {
    'MASS_DISPLACEMENT', 'No single MATLAB morphology builtin; derived from Nyxus centroid conventions';
    'COMPACTNESS', 'Nyxus-specific compactness formula';
    'DIAMETER_EQUAL_AREA', 'Nyxus formula differs from MATLAB EquivDiameter convention';
    'ROUNDNESS', 'Nyxus formula depends on its major-axis convention';
    'PERIMETER', 'Nyxus ordered contour convention differs from regionprops Perimeter';
    'DIAMETER_EQUAL_PERIMETER', 'Depends on Nyxus perimeter convention';
    'EDGE_MEAN_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'EDGE_STDDEV_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'EDGE_MAX_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'EDGE_MIN_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'EDGE_INTEGRATED_INTENSITY', 'No direct MATLAB feature over Nyxus ordered contour';
    'CIRCULARITY', 'Depends on Nyxus perimeter convention';
    'FRACT_DIM_BOXCOUNT', 'No direct MATLAB built-in oracle';
    'FRACT_DIM_PERIMETER', 'No direct MATLAB built-in oracle';
    'DIAMETER_MIN_ENCLOSING_CIRCLE', 'No direct MATLAB built-in oracle in base/IPT';
    'DIAMETER_CIRCUMSCRIBING_CIRCLE', 'No direct MATLAB built-in oracle in base/IPT';
    'DIAMETER_INSCRIBING_CIRCLE', 'No direct MATLAB built-in oracle in base/IPT';
    'ROI_RADIUS_MEAN', 'Nyxus-specific contour distance convention';
    'ROI_RADIUS_MAX', 'Nyxus-specific contour distance convention';
    'ROI_RADIUS_MEDIAN', 'Nyxus-specific contour distance convention';
    };
for i = 1:size(unsupported, 1)
    rows = add_status(rows, truth, unsupported{i, 1}, NaN, 'NO_BUILTIN', '', unsupported{i, 2});
end
for i = 1:8
    rows = add_status(rows, truth, sprintf('EXTREMA_P%d_X', i), NaN, 'NO_BUILTIN', '', 'No direct MATLAB regionprops extrema feature with Nyxus ordering');
    rows = add_status(rows, truth, sprintf('EXTREMA_P%d_Y', i), NaN, 'NO_BUILTIN', '', 'No direct MATLAB regionprops extrema feature with Nyxus ordering');
end

results = rows_to_table(rows);
outPath = fullfile(outDir, 'shape2d_builtin_independent_results.csv');
writetable(results, outPath);
print_result_rows(results, 'Shape2D');
fprintf('Shape2D MATLAB built-in independent rows: %d pass, %d mismatch, %d no-builtin. Wrote %s\n', ...
    sum(string(results.Status) == "PASS"), ...
    sum(string(results.Status) == "MISMATCH"), ...
    sum(string(results.Status) == "NO_BUILTIN"), outPath);
end

function rows = add_builtin(rows, truth, feature, actual, builtinName, note)
rows = add_status(rows, truth, feature, actual, '', builtinName, note);
end

function rows = add_diagnostic(rows, feature, actual, builtinName, note)
rows(end + 1, :) = {'Shape and morphology', feature, NaN, actual, NaN, NaN, NaN, 'DIAGNOSTIC', builtinName, note}; %#ok<AGROW>
end

function rows = add_status(rows, truth, feature, actual, forcedStatus, builtinName, note)
if ~isfield(truth, feature)
    return;
end
expected = truth.(feature);
tolerancePct = 5.0;
if strlength(string(forcedStatus)) > 0
    status = forcedStatus;
    absDiff = NaN;
    relDiffPct = NaN;
else
    absDiff = abs(actual - expected);
    relDiffPct = relative_diff_pct(actual, expected);
    if relDiffPct <= tolerancePct
        status = 'PASS';
    else
        status = 'MISMATCH';
    end
end
rows(end + 1, :) = {'Shape and morphology', feature, expected, actual, absDiff, relDiffPct, tolerancePct, status, builtinName, note}; %#ok<AGROW>
end

function results = rows_to_table(rows)
results = cell2table(rows, 'VariableNames', { ...
    'Family', 'Feature', 'Expected', 'MatlabBuiltinValue', 'AbsDiff', 'RelDiffPct', 'TolerancePct', 'Status', 'BuiltinUsed', 'Note'});
end

function truth = load_named_constants(headerPath)
txt = fileread(headerPath);
block = regexp(txt, 'shape2d_truth\s*\{([\s\S]*?)\};', 'tokens', 'once');
if isempty(block)
    error('Could not find shape2d_truth in %s', headerPath);
end
tokens = regexp(block{1}, '\{"([^"]+)",\s*([^}]+)\}', 'tokens');
truth = struct();
for i = 1:numel(tokens)
    truth.(tokens{i}{1}) = str2double(tokens{i}{2});
end
end

function roi = load_shape_fixture(dataPath)
txt = fileread(dataPath);
inten = load_nyxus_pixel_array(txt, 'shape2d_morphology_intensity');
mask = load_nyxus_pixel_array(txt, 'shape2d_morphology_mask');
roi = [];
for i = 1:size(mask, 1)
    if mask(i, 3) ~= 0
        idx = find(inten(:, 1) == mask(i, 1) & inten(:, 2) == mask(i, 2), 1);
        roi(end + 1, :) = [mask(i, 1), mask(i, 2), inten(idx, 3)]; %#ok<AGROW>
    end
end
end

function triples = load_nyxus_pixel_array(txt, arrayName)
pattern = sprintf('const\\s+static\\s+NyxusPixel\\s+%s\\[\\]\\s*=\\s*\\{([\\s\\S]*?)\\};', arrayName);
block = regexp(txt, pattern, 'tokens', 'once');
if isempty(block)
    error('Could not find %s', arrayName);
end
tokens = regexp(block{1}, '\{\s*([0-9.+\-eE]+)\s*,\s*([0-9.+\-eE]+)\s*,\s*([0-9.+\-eE]+)\s*\}', 'tokens');
triples = zeros(numel(tokens), 3);
for i = 1:numel(tokens)
    triples(i, :) = [str2double(tokens{i}{1}), str2double(tokens{i}{2}), str2double(tokens{i}{3})];
end
end

function [bw, inten, minX, minY] = fixture_to_images(roi)
x = roi(:, 1);
y = roi(:, 2);
minX = min(x);
minY = min(y);
maxX = max(x);
maxY = max(y);
bw = false(maxY - minY + 1, maxX - minX + 1);
inten = zeros(size(bw));
for i = 1:size(roi, 1)
    row = roi(i, 2) - minY + 1;
    col = roi(i, 1) - minX + 1;
    bw(row, col) = true;
    inten(row, col) = roi(i, 3);
end
end

function pct = relative_diff_pct(fval, gt)
if isnan(fval) || isnan(gt)
    pct = Inf;
elseif gt == 0
    if abs(fval) <= 1e-12
        pct = 0.0;
    else
        pct = Inf;
    end
else
    pct = 100.0 * abs(fval - gt) / abs(gt);
end
end

function print_result_rows(results, label)
for i = 1:height(results)
    fprintf('%s %-42s expected=%.15g matlab=%.15g absdiff=%.15g reldiff=%.6g%% tol=%.2f%% status=%s\n', ...
        label, results.Feature{i}, results.Expected(i), results.MatlabBuiltinValue(i), ...
        results.AbsDiff(i), results.RelDiffPct(i), results.TolerancePct(i), results.Status{i});
end
end
