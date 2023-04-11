function varargout = odefunc(B,par)
    % B - population abundance, concatinated from all patches
    % par - structure storing all parameters
    
    B(logical(par.extinSpe)) = 0;
    B(B < 0) = 0;
    Bmat = reshape(B,par.nSpecies,par.nPatch);
    %% functional response
    % F_ij - functional response of predator i to prey j
    F = cell(par.nPatch,1);
    for ii = 1:par.nPatch
        F{ii} = par.L .* (par.w .* par.a{ii} .* (repmat(Bmat(:,ii)',par.nSpecies,1) .^ par.q{ii}) ./...
            (1 + sum(par.w .* par.a{ii} .* par.h{ii} .* (repmat(Bmat(:,ii)',par.nSpecies,1) .^ par.q{ii}),2)));
    end
    %% r_birth - per capita birth rate
    r_birth = cell(par.nPatch,1);
    for ii = 1:par.nPatch
        % with plant competition
        r_birth{ii} = par.isPlt .* (par.r{ii} .* (1 - par.alpha * (Bmat(:,ii) .* par.isPlt) ./ par.K{ii}))...
            + par.isAnm .* sum(par.e{ii} .* F{ii},2);
    end
%     r_birth_mat = cell2mat(r_birth);
    %% r_death - per capita death rate
    r_death = cell(par.nPatch,1);
    for ii = 1:par.nPatch
        r_death{ii} = sum(-F{ii} .* Bmat(:,ii),1)' ./ (Bmat(:,ii) + eps)... % eaten by predator
            + par.isAnm .* par.r{ii}; % intrinsic death, animal only
    end
%     r_death_mat = cell2mat(r_death);
    %% migration matrix. with migratory cost considered but not implemented.
    %m_j(k1->k2)
    logm0 = zeros(par.nSpecies,par.nPatch,par.nPatch);
    for k1 = 1:par.nPatch
        for k2 = 1:par.nPatch
            if k1 == k2
                logm0(:,k1,k2) = (par.lambda .* ...
                    (par.beta * r_birth{k1} + par.gamma * r_death{k1}));
            else
                logm0(:,k1,k2) = (- par.lambda * par.epsilon);
            end
        end
    end
    m0 = exp(logm0 - max(logm0,[],3)); % exp(large x) would go to Inf. mathematiclaly same as exp(logm0), but computational doable
    m = par.m .* m0 ./ repmat(sum(m0,3),1,1,par.nPatch); % normalizeb
    %% r_migrt
    dBdt_migrt = cell(par.nPatch,1);
    for ii = 1:par.nPatch
        dBdt_migrt{ii} = -par.m .* Bmat(:,ii)... % emmigration from patch ii
            + sum(squeeze(m(:,:,ii)) .* Bmat,2); % immigration into patch ii
    end
    %% overall per capita growth (r_real) & dBdt
    r_local = cellfun(@(r_b,r_d) r_b+r_d,r_birth,r_death,'UniformOutput',false);
    dBdt = cell2mat(r_local) .* B + cell2mat(dBdt_migrt);
    %% output
    % species * patch. proportion of species stay in original patch
    stayRate = (1 - par.m) + m(:,eye(par.nPatch)==1);
    stayRateMeta = sum(stayRate .* Bmat,2) ./ (sum(Bmat,2) + eps);
    % except for dBdt, other terms of varargout are for debug/inspection
    varargout = {dBdt; 1 - stayRate; 1 - stayRateMeta; ...
        1-mean(stayRateMeta(stayRateMeta > 0)); 1-median(stayRateMeta(stayRateMeta > 0)); ...
        cell2mat(r_death); cell2mat(r_birth); cell2mat(dBdt_migrt); m; F};
end