function paraSet = parasetgenerator(nrep, varargin)
% Rely on files: foodwebgenerator
% nrep = 100

% set default values
par.nSpecies = 50; % Number of total species
par.nPatch = 3; % Number of habitat patches
range_C = [0.01, 0.05]; % range of foodweb connectance C
range_r_p = [0.10 0.30]; % range of plant growth rate r_p
range_r_a = -[0.02 0.06]; % range of animal growth rate r_a
range_K = [4 6]; % carrying capacity K
range_alpha = [0 0.2]; % plants' interspecies competition coefficient
range_a_a = [0.10 0.45]; % attack rate
range_a_p = [0.40 0.90];
range_e_a = [0.3 0.8]; % assimilation coefficient
range_e_p = [0.1 0.4];
range_q = [1.3 2]; % Hill coefficient
range_h = [0.1 0.5]; % handling time
par.extin = 1e-10; % Extinction threshold
B0 = ones(par.nSpecies * par.nPatch,1) * 1; % Initial species biomass
% handle options
if ~isempty(varargin)
    varkeys = varargin(1:2:numel(varargin));
    varvalues = varargin(2:2:numel(varargin));
    if numel(varkeys) > numel(varvalues)
        warning(['Argument key without value, use default: ', varkeys(end), '.'])
    end
    for ik = 1:numel(varvalues)
        key = varkeys{ik};
        value = varvalues{ik};
        if ~ischar(varkeys{ik})
            error('Argument key %d format error. Must be string.', ik)
        end
        switch key
            case 'nSpecies'
                par.nSpecies = value;
            case 'nPatch'
                par.nPatch = value;
            otherwise
                warning('Invalid argument key:%s', key)
        end
    end
end
            
% for each rep generate spatial foodwebs parameters
rng('shuffle');
randPars = cell(nrep,1);
for irep = 1:nrep
    % Connectance of food web
    connectance = range_C(1) + (range_C(2) - range_C(1)) * rand(1);
    % food web structure matrix
    par.L =  foodwebgenerator(par.nSpecies, connectance);
    par.isPlt = (sum(par.L,2) == 0);
    par.isAnm = (sum(par.L,2) > 0);
    
    % Innate growth rate
    par.r_mod = cell(par.nPatch,1); % final par.r = par.r * (1 + heter * par.r_mod)
    for ii = 1:par.nPatch
        par.r_mod{ii} = rand(par.nSpecies,1) * 2 - 1; % [-1, 1]
    end
    par.r = par.isPlt .* (range_r_p(1) +...
            (range_r_p(2) - range_r_p(1)) * rand(par.nSpecies,1))...
            + par.isAnm .* (range_r_a(1) +...
            (range_r_a(2) - range_r_a(1)) * rand(par.nSpecies,1));
    % Carrying capacity for each patch
    par.K_mod = cell(par.nPatch,1); % final par.K = par.K * (1 + heter * par.K_mod)
    for ii = 1:par.nPatch
        par.K_mod{ii} = rand(par.nSpecies,1) * 2 - 1; % [-1 1]
    end
    par.K = range_K(1) + (range_K(2) - range_K(1)) * rand(par.nSpecies,1);
    % competitive coefficients between plants
    par.alpha = range_alpha(1) + rand(par.nSpecies) * (range_alpha(2) - range_alpha(1));
    par.alpha(logical(eye(par.nSpecies))) = 1; % intraspecies
    % Attack rate
    par.a = cell(par.nPatch,1); % predator * prey
    par.a(:) = {par.L .* (par.isPlt' .* (range_a_p(1) +...
        (range_a_p(2) - range_a_p(1)) * rand(par.nSpecies))...
        + par.isAnm' .* (range_a_a(1) +...
        (range_a_a(2) - range_a_a(1)) * rand(par.nSpecies)))};
    % assimilation coefficient
    par.e = cell(par.nPatch,1);% pdt * prey
    par.e(:) = {par.L .* (par.isPlt' .* (range_e_p(1) +...
        (range_e_p(2) - range_e_p(1)) * rand(par.nSpecies))...
        + par.isAnm' .* (range_e_a(1) +...
        (range_e_a(2) - range_e_a(1)) * rand(par.nSpecies)))};
    % Hill coefficient
    par.q = cell(par.nPatch,1);
    par.q(:) = {range_q(1) + (range_q(2) - range_q(1)) * rand(par.nSpecies)};
    % Handling time
    par.h = cell(par.nPatch,1);
    par.h(:) = {range_h(1) + (range_h(2) - range_h(1)) * rand(par.nSpecies)};
    % Prey preference
    par.w = par.L ./ (sum(par.L, 2) + eps);
    % Migration strategy parameter: seek growth
    par.beta = 1;
    % Migration strategy parameter: avoid death
    par.gamma = 1;
    
    randPars{irep} = par;
end
paraSet = randPars;
end