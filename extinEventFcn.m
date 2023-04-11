function [value,isterminal,direction] = extinEventFcn(t,B,par)
    Bmat = reshape(B,par.nSpecies,par.nPatch);
    
    par.m = par.m .* ones(par.nSpecies,1);
    value = Bmat - par.extin;
    for is = 1:par.nSpecies
        if par.m(is) > 0
            value(is,:) = max(value(is,:)); % if connected by migration, then only consider global extinction
        end
    end
    
    value(logical(par.extinSpe)) = 1; % species that are already set to extinction
    isterminal = ones(par.nSpecies,par.nPatch);
    direction = zeros(par.nSpecies,par.nPatch);
end
