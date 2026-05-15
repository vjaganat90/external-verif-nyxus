function results = audit_unsupported_builtin_families(~, outDir)
% Record 2D feature families for which MATLAB has no comparable high-level
% built-in feature oracle.

if nargin < 2 || strlength(string(outDir)) == 0
    outDir = fileparts(mfilename('fullpath'));
end

rows = {
    'GLRLM', 32, 'NO_BUILTIN', 'MATLAB Image Processing Toolbox has no full GLRLM matrix/feature family oracle';
    'GLDZM', 18, 'NO_BUILTIN', 'MATLAB has no GLDZM feature-family oracle';
    'GLSZM', 16, 'NO_BUILTIN', 'MATLAB can label connected zones, but has no GLSZM feature-family oracle';
    'GLDM', 14, 'NO_BUILTIN', 'MATLAB has no GLDM feature-family oracle';
    'NGLDM', 19, 'NO_BUILTIN', 'MATLAB has no NGLDM feature-family oracle';
    'NGTDM', 5, 'NO_BUILTIN', 'MATLAB has no NGTDM feature-family oracle';
    'Gabor score', 1, 'PARTIAL_BUILTIN', 'MATLAB has gabor/imgaborfilt but not Nyxus custom CPU score and threshold statistic';
    'Neighbor / touching', 9, 'NO_BUILTIN', 'MATLAB has geometry primitives but no Nyxus multi-ROI neighbor/touching feature oracle';
    'Advanced shape long tail', NaN, 'PARTIAL_BUILTIN', 'Some primitives exist, but Feret/Martin/Nassenstein/chords/polygonality/geodesic/thickness need separate oracles';
    'Moments / Zernike / radial', NaN, 'PARTIAL_BUILTIN', 'Some moment primitives exist in MATLAB/File Exchange ecosystems, but no single trusted built-in oracle is wired here';
    };

results = cell2table(rows, 'VariableNames', {'Family', 'FeatureCount', 'Status', 'Reason'});
outPath = fullfile(outDir, 'unsupported_builtin_families.csv');
writetable(results, outPath);
fprintf('Unsupported/partial MATLAB built-in family audit wrote %s\n', outPath);
end
