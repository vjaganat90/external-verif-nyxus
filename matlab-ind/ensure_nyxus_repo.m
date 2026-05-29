function repoRoot = ensure_nyxus_repo(explicitPath)
% ENSURE_NYXUS_REPO  Locate a Nyxus checkout for V&V fixture parsing.
%
% The MATLAB verifier parses golden test constants and ROI fixtures out of the
% Nyxus C++ test headers. This helper deliberately prefers an existing checkout
% and does not fetch, pull, checkout, or reset it.
%
% Resolution order:
%   1. explicitPath argument
%   2. NYXUS_REPO_DIR
%   3. common local locations, including ../nyxus and existing .cache/nyxus
%   4. optional one-time clone into .cache/nyxus when NYXUS_ALLOW_CLONE=1
%
% Optional clone knobs:
%   NYXUS_REPO_URL : git URL to clone from, default HTTPS maintainer remote.
%   NYXUS_BRANCH   : branch to clone; omitted means the remote default branch.

if nargin < 1
    explicitPath = '';
end

% ---- 1. explicit path argument --------------------------------------------
if ~isempty(explicitPath)
    repoRoot = require_nyxus_checkout(explicitPath, 'explicit Nyxus path');
    return;
end

% ---- 2. NYXUS_REPO_DIR ------------------------------------------------------
envDir = getenv('NYXUS_REPO_DIR');
if ~isempty(envDir)
    repoRoot = require_nyxus_checkout(envDir, 'NYXUS_REPO_DIR');
    return;
end

% ---- 3. common local checkouts ---------------------------------------------
thisDir = fileparts(mfilename('fullpath'));   % .../matlab-ind
rootDir = fileparts(thisDir);                 % external-verif repo root
parentDir = fileparts(rootDir);
cacheDir = fullfile(rootDir, '.cache', 'nyxus');

candidates = {
    fullfile(parentDir, 'nyxus')
    pwd
    fullfile(fileparts(pwd), 'nyxus')
    cacheDir
    };

for i = 1:numel(candidates)
    candidate = candidates{i};
    if looks_like_nyxus_checkout(candidate)
        repoRoot = candidate;
        return;
    end
end

% ---- 4. optional one-time clone --------------------------------------------
if allow_clone()
    if isfolder(cacheDir)
        repoRoot = require_nyxus_checkout(cacheDir, 'existing .cache/nyxus');
        return;
    end

    repoUrl = getenv('NYXUS_REPO_URL');
    if isempty(repoUrl)
        repoUrl = 'https://github.com/vjaganat90/nyxus.git';
    end

    branch = getenv('NYXUS_BRANCH');
    branchArg = '';
    if ~isempty(branch)
        branchArg = sprintf(' --branch "%s"', branch);
    end

    if ~isfolder(fileparts(cacheDir))
        mkdir(fileparts(cacheDir));
    end
    run_git(sprintf('git clone --depth 1%s "%s" "%s"', branchArg, repoUrl, cacheDir));
    repoRoot = require_nyxus_checkout(cacheDir, 'newly cloned .cache/nyxus');
    return;
end

searched = strjoin(candidates, newline + "  - ");
error('ensure_nyxus_repo:notFound', ...
    ['Could not find a Nyxus checkout for V&V fixture parsing.' newline ...
     'Pass repoRoot to run_all_matlab_independent_2d_checks, set NYXUS_REPO_DIR, ' ...
     'place a checkout at ../nyxus, or set NYXUS_ALLOW_CLONE=1 to clone once into .cache/nyxus.' newline ...
     'Searched:' newline '  - %s'], searched);
end

function run_git(cmd)
fprintf('[ensure_nyxus_repo] $ %s\n', cmd);
[status, out] = system(cmd);
if status ~= 0
    error('ensure_nyxus_repo:gitFailed', 'git command failed (%d):\n%s', status, out);
end
end

function tf = looks_like_nyxus_checkout(path)
tf = isfolder(path) && isfile(fullfile(path, 'tests', 'test_data.h'));
end

function repoRoot = require_nyxus_checkout(path, sourceName)
if looks_like_nyxus_checkout(path)
    repoRoot = path;
else
    error('ensure_nyxus_repo:notNyxus', ...
        '%s points to %s, but tests/test_data.h was not found.', sourceName, path);
end
end

function tf = allow_clone()
value = lower(strtrim(getenv('NYXUS_ALLOW_CLONE')));
tf = any(strcmp(value, {'1', 'true', 'yes', 'on'}));
end
