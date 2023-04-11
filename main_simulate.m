function main_simulate(snr, vary, repList)
% snr - simulation scenario (see main.m)
% vary - {'H', 'm0'}
%   'H' - varying habitat heterogeneity H (and fix m0=0.1)
%   'm0' - varying habitat conectivity m0 (and fix H=0.6)

    %% setup different simulation scenarios

    % main scenario
    nSpecies = 50;
    nPatch = 3;
    pltMigStr = 'stra'; % plant also disperse with fitness-dependency/strategy
    epsilon = 1e-5;

    % modified scenarios
    switch(snr)
        case ''
        case 'pltiso_'
            pltMigStr = 'iso';
        case 'pltrand_'
            pltMigStr = 'rand';
        case 'p6_'
            nPatch = 6;
        case 'noepsilon_'
            epsilon = 0;
        otherwise
            error('Undefined simulation scenario: %s', snr)
    end

    %% generate or load randPars (random basal foodwebs)

    % % if generates new random basal foodwebs
    % randPars = parasetgenerator(nRep, 'nSpecies', nSpecies, 'nPatch', nPatch);
    % save(sprintf('randPars_s%dp%dr%d_%s.mat', nSpecies, nPatch, nRep, ...
    %     datetime('now','format','MMddHHmmss')),'randPars');

    % if use existing basal foodwebs
    if nPatch == 3
        load('./randPars_s50p3r1000_0218062416.mat');
    elseif nPatch == 6
        load('./randPars_s50p6r100_0815225637.mat');
    end

    if strcmp(vary, 'H')
        %% varying habitat heterogeneity level H (and fix m0=0.1)
        % define scenarios of spatial process: as dispersal mode (7, 6*lambda+isolate) * heter (6) array
        pltMigStr = 'stra'; % stra, iso, rand
        nRep = numel(repList);
        nLam = 6;
        nHet = 6;

        % potential migration rate / habitat connectivity m0
        mArr = repmat([0;ones(nLam,1)*0.1],1,nHet);
        % fitness-dependency lambda
        lambArr = repmat([0, 0, 1, 10, 1e2, 1e3, 1e4]',1,nHet);
        % habitat heterogeneity H
        heterArr = repmat([0:0.2:1],1+nLam,1); 

        for ir = repList
            % basal foodweb parameters
            parBase = randPars{ir};
            parBase.epsilon = epsilon;  % tendency to stay

            [BfinalSum, BfinalCVSum] = sim_scenarios(parBase, lambArr, mArr, heterArr, pltMigStr);
            save(sprintf('./rawdata/%srep%03d_diff_H.mat',snr,ir),...
                'BfinalSum', 'BfinalCVSum', 'parBase')
            clear BfinalSum BfinalCVSum

            % progress bar
            fprintf('.')
            if (mod(ir,100) == 0)
                fprintf('%d\n',ir);
            end
        end
    elseif strcmp(vary, 'm0')
       %% varying habitat connectivity level m0 (and fix H=0.6)
        % define scenarios of spatial process: as lambda (6) * m0 (13) array
        mArr = repmat(10.^(-4:1/3:0),6,1); % migration rate m0
        lambArr = repmat([0, 1, 10, 1e2, 1e3, 1e4]',1,13); % fitness dependency lambda
        heterArr = repmat([0.6],6,13); % habitat heterogeneity H

        for ir = repList
            % basal foodweb parameters
            parBase = randPars{ir};
            parBase.epsilon = epsilon;

            % plant dispersal: fitness dependent
            [BfinalSum, BfinalCVSum] = sim_scenarios(parBase, lambArr, mArr, heterArr, pltMigStr);
            save(sprintf('./rawdata/%srep%03d_diff_m0.mat',snr,ir),...
                'BfinalSum', 'BfinalCVSum', 'parBase')
            clear BfinalSum BfinalCVSum

            % progress bar
            fprintf('.')
            if (mod(ir,100) == 0)
                fprintf('%d\n',ir);
            end
        end
    end

end