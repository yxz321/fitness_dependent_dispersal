function [dist_mat] = calcDissimilarity(B, method)
    [nSpecies, nPatch] = size(B);
    w_sp = any(B > 0, 2);  % =0 if species extinct
    
    if method == 'Go'  % Gower, J.C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 27, 857–871.
        rang_sp = max(B, [], 2);  % maximum population size of a species, used to rescale population size
        B_cell = repmat(mat2cell(B, size(B, 1), ones(1, size(B, 2))),nPatch,1);
        dist_mat = cellfun(@(B_i, B_j) ...
            nansum(w_sp .* abs(B_i - B_j) ./ rang_sp) / sum(w_sp), ...
            B_cell, B_cell');
    elseif method == 'MGo'
        B_cell = repmat(mat2cell(B, size(B, 1), ones(1, size(B, 2))),nPatch,1);
        dist_mat = cellfun(@(B_i, B_j) ...
            nansum(w_sp .* abs(B_i - B_j)) / sum(w_sp), ...
            B_cell, B_cell');
    elseif method == 'An'  % Anderson, M.J. et. al. (2006). Multivariate dispersion as a measure of beta diversity
        B(B==0) = 1e-6;  % must assign a number for fold difference between existence-nonexistence
        B_mod = log(B);
        B_mod_cell = repmat(mat2cell(B_mod, size(B, 1), ones(1, size(B, 2))),nPatch,1);
        dist_mat = cellfun(@(B_i, B_j) ...
            sum(w_sp .* abs(B_i - B_j)) / sum(w_sp), ...
            B_mod_cell, B_mod_cell');
    else
        error('Undefined dissimilarity metrics!'); 
    end
end