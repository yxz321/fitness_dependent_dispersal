function processdata(snr)
    %% define data structure: lambda * heter * m0 * rep
    nSpecies = 50;
    nPatch = 3;

    lambList = [0,1,10,1e2,1e3,1e4]; % fitness-dependency lambda
    lamb_labels = {'0','10^0','10^1','10^2','10^3','10^4'};
    heterList = [0:0.2:1]; % habitat heterogeneity H
    mList = [0,10.^(-4:1/3:0)]; % migration rate m0
    if strcmp(snr, '')
        repList = 1:500;
    else
        repList = 1:50;
    end
    nLam = numel(lambList);
    nHet = numel(heterList);
    nMig = numel(mList);
    nRep = numel(repList);

    % data available on : (:, :, 11, :) + (:, :, 1, :) and (:, 4, :, :)
    % (:, :, 1, :) + (:, :, 11, :): m0=0 or 0.1, focus on lambda (+isolate) * heter interaction
    % (:, 4, :, :): H=0.6, focus on lambda * m0 (migration rate) interaction
    metrixKey = {'alphaDiv', 'betaDiv', 'gammaDiv'};
    for fn = metrixKey
        res.(fn{1}) = nan(nLam, nHet, nMig, nRep);
    end

    % in simulation data strategy are linearly numbered
    % straIdx & supStraIdx indecates the corresponding strategy no. for each position in new structure
    % straIdx is for simulations with different H
    % raw data files: '%srep%03d_diff_H.mat' % snr, irep
    straIdx = nan(nLam, nHet, nMig);
    straIdx(:, :, 1) = repmat([1:7:42], 6, 1); % m0 = 0, lambda doesn't matter
    straIdx(:, :, 11) = [2:7]' + [0:5] * 7; % m0 = 0.1

    % supStraIdx is for simulations with different m0
    % raw data files: '%srep%03d_diff_m0.mat' % snr, irep
    supStraIdx = nan(nLam, nHet, nMig);
    supStraIdx(:, 4, 2:end) = reshape(1:78, 6, 1, 13); % m = 0.01, 0.5, 1. in datafiles sup_rep***

    % functions for calculating corresponding indexes
    cutoff = 1e-6;  % global cutoff
    cutoff_local = 1e-1;  % local cutoff
    func_rescale = @(x_B) x_B ./ (max(x_B, [], 2) + eps);
    func.gammaDiv = @(x_B) sum(max(x_B, [], 2) > cutoff);
    func_mean_dist = @(dist) mean(dist(triu(true(size(dist)), 1)));
    func.betaDiv = @(x_B) func_mean_dist(calcDissimilarity(x_B, 'Go'));
    func.alphaDiv = @(x_B) mean(sum(func_rescale(x_B) > cutoff_local));

    %% load diff_H and diff_m0 data files, calculate diversity
    for ir = repList
        % temporary summary data, for one replication only
        flag_incomplete_rep = 0;
        clear temp
        for fn = metrixKey
            temp.(fn{1}) = nan(nLam, nHet, nMig);
        end

        % load varying heterogeneity data file (dispersal strategy * heter)
        clear BfinalSum
        try
            load(sprintf('./rawdata/%srep%03d_diff_H.mat', snr, ir))
        catch ME
            % some reps might be missing due to integration error
            if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                continue
            else
                rethrow(ME)
            end
        end
        % summarize diff_H data files
        for is = 1:42
            if isempty(BfinalSum{is})  % integration failure
                flag_incomplete_rep = 1;
                break
            end
            BfinalSum{is} = BfinalSum{is} .* (BfinalSum{is} > cutoff);
            for fn = metrixKey
    %             temp.(fn{1})(straIdx == is) = func.(fn{1})(BfinalSum{is});
                temp.(fn{1})(straIdx == is) = func.(fn{1})(BfinalSum{is});
            end
        end
        if flag_incomplete_rep
            continue
        end
        B_iso = BfinalSum{22};  % H = 0.6
        range_B_iso = nanmax(B_iso, [], 2);
        % load varying migration rate data file (dispersal strategy * migration rate)
    %     clear BfinalSum
        try
            load(sprintf('./rawdata/%srep%03d_diff_m0.mat', snr, ir))
        catch ME
            % some reps might be missing due to integration error
            if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                for fn = metrixKey
                    res.(fn{1})(:,:,:,ir) = temp.(fn{1});
                end
                continue
            else
                rethrow(ME)
            end
        end
        % summarize diff_m0 data file
        for is = 1:78
            if isempty(BfinalSum{is})
                flag_incomplete_rep = 1;
                break
            end
            BfinalSum{is} = BfinalSum{is} .* (BfinalSum{is} > cutoff);
            for fn = metrixKey
                temp.(fn{1})(supStraIdx == is) = func.(fn{1})(BfinalSum{is});
            end
        end

        if flag_incomplete_rep
            continue
        end

        % incorporate temporary summary matrix to final result summary matrix
        for fn = metrixKey
            res.(fn{1})(:,:,:,ir) = temp.(fn{1});
        end
    end

    % save(sprintf('./%sres.mat',snr), 'res')
end