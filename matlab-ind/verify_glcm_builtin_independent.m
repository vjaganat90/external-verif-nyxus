function results = verify_glcm_builtin_independent(repoRoot, outDir)
% Compare the MATLAB-exposed subset of Nyxus 2D GLCM expectations to
% graycomatrix/graycoprops. MATLAB does not expose the full Nyxus GLCM family.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = '/Users/jaganathv2/Code/Codex/nyxus';
end
if nargin < 2 || strlength(string(outDir)) == 0
    outDir = fileparts(mfilename('fullpath'));
end

dataPath = fullfile(repoRoot, 'tests', 'test_data.h');
headerPath = fullfile(repoRoot, 'tests', 'test_glcm.h');

slices = load_ibsi_phantom_slices(dataPath);
truth = load_glcm_truth(headerPath);

% graycomatrix offsets are [row_delta col_delta]. Nyxus stores the tested
% directions as (dx,dy): (1,0), (1,1), (0,1), (-1,1), so convert them here.
offsets = [0 1; 1 1; 1 0; 1 -1];
sliceVals = cell(numel(slices), 1);
for i = 1:numel(slices)
    I = slices{i};
    D = bin_image_nyxus_matlab_style(I, 100);
    glcms = graycomatrix(D, 'Offset', offsets, 'NumLevels', 100, ...
        'GrayLimits', [1 100], 'Symmetric', false);
    props = graycoprops(glcms, {'Contrast', 'Correlation', 'Energy', 'Homogeneity'});
    sliceVals{i} = props;
end

rows = {};
rows = add_glcm(rows, truth, 'GLCM_CONTRAST', mean_over_slices(sliceVals, 'Contrast'), 'Nyxus binning + graycomatrix + graycoprops Contrast', 'MATLAB-defined GLCM contrast after Nyxus test discretization');
rows = add_glcm(rows, truth, 'GLCM_CORRELATION', mean_over_slices(sliceVals, 'Correlation'), 'Nyxus binning + graycomatrix + graycoprops Correlation', 'MATLAB-defined GLCM correlation after Nyxus test discretization');
rows = add_glcm(rows, truth, 'GLCM_ASM', mean_over_slices(sliceVals, 'Energy'), 'Nyxus binning + graycoprops Energy', 'MATLAB Energy is angular second moment after Nyxus test discretization');
rows = add_glcm(rows, truth, 'GLCM_ENERGY', mean_over_slices(sliceVals, 'Energy'), 'Nyxus binning + graycoprops Energy', 'Nyxus names this output ENERGY/ASM in 2D tests');
rows = add_glcm(rows, truth, 'GLCM_HOM1', mean_over_slices(sliceVals, 'Homogeneity'), 'Nyxus binning + graycoprops Homogeneity', 'MATLAB homogeneity after Nyxus test discretization');
rows = add_glcm(rows, truth, 'GLCM_CONTRAST_AVE', mean_over_slices(sliceVals, 'Contrast'), 'Nyxus binning + graycomatrix + graycoprops Contrast', 'MATLAB-defined averaged GLCM contrast after Nyxus test discretization');
rows = add_glcm(rows, truth, 'GLCM_CORRELATION_AVE', mean_over_slices(sliceVals, 'Correlation'), 'Nyxus binning + graycomatrix + graycoprops Correlation', 'MATLAB-defined averaged GLCM correlation after Nyxus test discretization');
rows = add_glcm(rows, truth, 'GLCM_ASM_AVE', mean_over_slices(sliceVals, 'Energy'), 'Nyxus binning + graycoprops Energy', 'MATLAB averaged Energy is angular second moment after Nyxus test discretization');
rows = add_glcm(rows, truth, 'GLCM_ENERGY_AVE', mean_over_slices(sliceVals, 'Energy'), 'Nyxus binning + graycoprops Energy', 'Nyxus names this averaged output ENERGY/ASM in 2D tests');
rows = add_glcm(rows, truth, 'GLCM_HOM1_AVE', mean_over_slices(sliceVals, 'Homogeneity'), 'Nyxus binning + graycoprops Homogeneity', 'MATLAB averaged homogeneity after Nyxus test discretization');

unsupported = {
    'GLCM_ACOR', 'Not exposed by graycoprops';
    'GLCM_CLUPROM', 'Not exposed by graycoprops';
    'GLCM_CLUSHADE', 'Not exposed by graycoprops';
    'GLCM_CLUTEND', 'Not exposed by graycoprops';
    'GLCM_DIFAVE', 'Not exposed by graycoprops';
    'GLCM_DIFENTRO', 'Not exposed by graycoprops';
    'GLCM_DIFVAR', 'Not exposed by graycoprops';
    'GLCM_DIS', 'Not exposed by graycoprops';
    'GLCM_ENTROPY', 'Not exposed by graycoprops';
    'GLCM_HOM2', 'Not exposed by graycoprops';
    'GLCM_ID', 'Not exposed by graycoprops';
    'GLCM_IDN', 'Not exposed by graycoprops';
    'GLCM_IDM', 'Not exposed by graycoprops';
    'GLCM_IDMN', 'Not exposed by graycoprops';
    'GLCM_INFOMEAS1', 'Not exposed by graycoprops';
    'GLCM_INFOMEAS2', 'Not exposed by graycoprops';
    'GLCM_IV', 'Not exposed by graycoprops';
    'GLCM_JAVE', 'Not exposed by graycoprops';
    'GLCM_JE', 'Not exposed by graycoprops';
    'GLCM_JMAX', 'Not exposed by graycoprops';
    'GLCM_JVAR', 'Not exposed by graycoprops';
    'GLCM_SUMAVERAGE', 'Not exposed by graycoprops';
    'GLCM_SUMENTROPY', 'Not exposed by graycoprops';
    'GLCM_SUMVARIANCE', 'Not exposed by graycoprops';
    'GLCM_VARIANCE', 'Not exposed by graycoprops';
    };
for i = 1:size(unsupported, 1)
    rows = add_status(rows, truth, unsupported{i, 1}, NaN, 'NO_BUILTIN', '', unsupported{i, 2});
end
for i = 1:size(unsupported, 1)
    feature = unsupported{i, 1};
    if strcmp(feature, 'GLCM_HOM2')
        continue;
    end
    rows = add_status(rows, truth, [feature '_AVE'], NaN, 'NO_BUILTIN', '', [unsupported{i, 2} ' for the averaged Nyxus GLCM feature']);
end

results = rows_to_table(rows);
outPath = fullfile(outDir, 'glcm_builtin_independent_results.csv');
writetable(results, outPath);
print_result_rows(results, 'GLCM');
fprintf('GLCM MATLAB built-in independent rows: %d pass, %d mismatch, %d no-builtin. Wrote %s\n', ...
    sum(string(results.Status) == "PASS"), ...
    sum(string(results.Status) == "MISMATCH"), ...
    sum(string(results.Status) == "NO_BUILTIN"), outPath);
end

function val = mean_over_slices(sliceVals, propName)
vals = [];
for i = 1:numel(sliceVals)
    vals = [vals, sliceVals{i}.(propName)]; %#ok<AGROW>
end
val = mean(vals(:));
end

function rows = add_glcm(rows, truth, feature, actual, builtinName, note)
rows = add_status(rows, truth, feature, actual, '', builtinName, note);
end

function rows = add_status(rows, truth, feature, actual, forcedStatus, builtinName, note)
if ~isfield(truth, feature)
    truthFeature = glcm_truth_key(feature);
    if ~isfield(truth, truthFeature)
        return;
    end
else
    truthFeature = feature;
end
tolerancePct = 5.0;
expected = truth.(truthFeature);
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
rows(end + 1, :) = {'GLCM', feature, expected, actual, absDiff, relDiffPct, tolerancePct, status, builtinName, note}; %#ok<AGROW>
end

function results = rows_to_table(rows)
results = cell2table(rows, 'VariableNames', { ...
    'Family', 'Feature', 'Expected', 'MatlabBuiltinValue', 'AbsDiff', 'RelDiffPct', 'TolerancePct', 'Status', 'BuiltinUsed', 'Note'});
end

function truthKey = glcm_truth_key(feature)
suffix = '_AVE';
if numel(feature) > numel(suffix) && strcmp(feature(end - numel(suffix) + 1:end), suffix)
    truthKey = feature(1:end - numel(suffix));
else
    truthKey = feature;
end
end

function truth = load_glcm_truth(headerPath)
txt = fileread(headerPath);
block = regexp(txt, 'glcm_values\s*\{([\s\S]*?)\};', 'tokens', 'once');
if isempty(block)
    error('Could not find glcm_values in %s', headerPath);
end
tokens = regexp(block{1}, '\{"([^"]+)",\s*([^}]+)\}', 'tokens');
truth = struct();
for i = 1:numel(tokens)
    truth.(tokens{i}{1}) = str2double(tokens{i}{2});
end
end

function slices = load_ibsi_phantom_slices(dataPath)
txt = fileread(dataPath);
slices = {};
for z = 1:4
    inten = load_nyxus_pixel_array(txt, sprintf('ibsi_phantom_z%d_intensity', z));
    mask = load_nyxus_pixel_array(txt, sprintf('ibsi_phantom_z%d_mask', z));
    slices{end + 1} = build_masked_aux_image(inten, mask); %#ok<AGROW>
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

function D = bin_image_nyxus_matlab_style(auxImage, nLevels)
% Nyxus GLCM tests intentionally use its MATLAB-style preprocessing:
% zero/background maps to level 1, positive intensities are scaled by the
% slice maximum, and feature formulas are then evaluated over levels 1..N.
maxI = max(auxImage(:));
if maxI <= 0
    error('Expected positive max intensity for GLCM verifier.');
end

D = zeros(size(auxImage));
slope = double(nLevels) / double(maxI);
for idx = 1:numel(auxImage)
    x = auxImage(idx);
    if x == 0
        D(idx) = 1;
    else
        y = floor(slope * double(x) + 1.0);
        D(idx) = min(max(y, 1), nLevels);
    end
end
end

function auxImage = build_masked_aux_image(intensityTriples, maskTriples)
x = maskTriples(:, 1);
y = maskTriples(:, 2);
minX = min(x);
minY = min(y);
maxX = max(x);
maxY = max(y);
auxImage = zeros(maxY - minY + 1, maxX - minX + 1);
for i = 1:size(maskTriples, 1)
    if maskTriples(i, 3) == 0
        continue;
    end
    idx = find(intensityTriples(:, 1) == maskTriples(i, 1) & intensityTriples(:, 2) == maskTriples(i, 2), 1);
    row = maskTriples(i, 2) - minY + 1;
    col = maskTriples(i, 1) - minX + 1;
    auxImage(row, col) = intensityTriples(idx, 3);
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
