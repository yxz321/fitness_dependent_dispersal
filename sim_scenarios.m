function [BfinalSum, BfinalCVSum] = sim_scenarios(parBase, lambArr, mArr, heterArr, pltMigStr)
% [BfinalSum, BfinalCVSum] = sim_scenarios(parBase, lambArr, mArr, heterArr, pltMigStr)
% simulate one foodweb (defined by parBase) across different dispersal scenarios
% i.e. different lambda, m, H, and plant dispersal scenario
% ---
% Input parameters
% parBase - structure. foodweb parameters + initial disperal parameters
%       the latter would be modified during run
% lambArr, mArr, heterArr - 
%       list/matrix (would be flattened) of dispersal parameters
%       should be of the same shape (NOTE: not list of unique values!)
%       each entry stands for the dispersal parameters for one dispersal mode
%       e.g. for position "is", lambArr[is] + mArr[is] + heterArr[is]
%            together define one dispersal mode
%
%       lambArr - fitness dependency, lambda
%       mArr - migration rate, m0
%       heterArr - habitat heterogeneity, H
% ---
% Output parameters
% BfinalSum - cell(nStr,1), nStr is the flattend size of lambArr/mArr/heterArr
%       cell "is" is the equilibrium biomass (nSpecies * nPatch matrix) for dispersal mode "is"
% BfinalCVSum = cell(nStr,1), nStr is the flattend size of lambArr/mArr/heterArr
%       cell "is" is the coefficient of variance of biomass in last 1000
%           time steps (nSpecies * nPatch matrix) for dispersal mode "is"

    % change integration failure warnings into errors for handling
    s1 = warning('error','MATLAB:ode15s:IntegrationTolNotMet');
    s2 = warning('error','MATLAB:illConditionedMatrix');

    % option used for ode solver (simulation to equilibrium)
    odeopt = odeset('AbsTol',1e-11,'RelTol',1e-4,'NonNegative',...
        ones(parBase.nSpecies * parBase.nPatch,1));
    B0 = ones(parBase.nSpecies * parBase.nPatch,1) * 1; % initial biomass

    nStr = numel(lambArr);
    % equilibrium biomass
    BfinalSum = cell(nStr,1); 
    % coefficient of variance of biomass in last 1000 time steps
    BfinalCVSum = cell(nStr,1);
    
    for is = 1:nStr
       %% generate parameters of spatial foodwebs under each scenario
        par = parBase;
        if strcmp(pltMigStr,'iso')  % plant dispersal scenario: isolate
            par.m = mArr(is) * par.isAnm;
        else
            par.m = mArr(is);
        end
        if strcmp(pltMigStr,'rand') % plant dispersal scenario: random
            par.lambda = lambArr(is) * par.isAnm;
        else
            par.lambda = lambArr(is);
        end
        par.heter = heterArr(is);
        % apply heterogeneity. K * (1 + heter * K_mod)
        par.K = cellfun(@(dK) par.K .* (1 + par.heter * dK),...
            par.K_mod,'UniformOutput',false);
        par.r = cellfun(@(dr) par.r .* (1 + par.heter * par.isPlt .* dr),...
            par.r_mod,'UniformOutput',false);
       
        %% Simulate to equilibrium
        try
            % whenever there's global extinction, stop integration, set population to 0, then start again
            par.extinSpe = zeros(par.nSpecies * par.nPatch,1);
            [~,b_end] = ode15s_withevent(B0,par,[0 1e5],odeopt,'fun',@odefunc);
            % record the last few time points to check for convergence
            [~,b] = ode15s_withevent(b_end(end,:)',par,[1:1000],odeopt,'fun',@odefunc);

            Bfinal = reshape(mean(b,1),par.nSpecies,par.nPatch);
            BfinalCV = reshape(var(b,1),par.nSpecies,par.nPatch);
            BfinalSum{is} = Bfinal;
            BfinalCVSum{is} = BfinalCV;
        catch ME
            if strcmp(ME.identifier,'MyError:TimeLimitReached') || strcmp(ME.identifier,'MATLAB:ode15s:IntegrationTolNotMet')...
                    || strcmp(ME.identifier,'MATLAB:illConditionedMatrix')
                fprintf('/')
                break %  in case of integration failure, stop and abandon the current rep but continue the whole program
            else
                rethrow(ME)
            end
        end
    end
    % restoring back to warnings
    warning(s1)
    warning(s2)
end
