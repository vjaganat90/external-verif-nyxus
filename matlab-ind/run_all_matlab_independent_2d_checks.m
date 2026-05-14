function summary = run_all_matlab_independent_2d_checks(repoRoot, outDir)
% Run all 2D MATLAB built-in independent checks that currently exist.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = default_repo_root();
end
if nargin < 2 || strlength(string(outDir)) == 0
    outDir = fileparts(mfilename('fullpath'));
end

repoRoot = char(repoRoot);
outDir = char(outDir);

fprintf('\nNyxus 2D MATLAB built-in independent verification\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Output dir: %s\n\n', outDir);

firstOrder = verify_first_order_builtin_independent(repoRoot, outDir);
shape2d = verify_shape2d_builtin_independent(repoRoot, outDir);
glcm = verify_glcm_builtin_independent(repoRoot, outDir);
unsupported = audit_unsupported_builtin_families(repoRoot, outDir);

summary = [
    summarize_family('First-order intensity', firstOrder)
    summarize_family('Shape and morphology subset', shape2d)
    summarize_family('GLCM MATLAB-exposed subset', glcm)
    ];

summary = [summary; { ...
    'Families without MATLAB built-in oracle', ...
    height(unsupported), ...
    0, ...
    0, ...
    height(unsupported), ...
    'See unsupported_builtin_families.csv'}];

summaryTable = cell2table(summary, 'VariableNames', { ...
    'Family', 'CheckedRows', 'PassRows', 'MismatchRows', 'NoBuiltinRows', 'Notes'});
summaryPath = fullfile(outDir, 'matlab_independent_2d_summary.csv');
writetable(summaryTable, summaryPath);

disp(summaryTable);
fprintf('\nWrote summary: %s\n', summaryPath);
end

function row = summarize_family(name, t)
statuses = string(t.Status);
row = { ...
    name, ...
    height(t), ...
    sum(statuses == "PASS"), ...
    sum(statuses == "MISMATCH"), ...
    sum(statuses == "NO_BUILTIN"), ...
    'MATLAB built-in oracle comparison'};
end

function repoRoot = default_repo_root()
repoRoot = '/Users/jaganathv2/Code/Codex/nyxus';
end
