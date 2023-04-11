function [res_t,res_b] = ode15s_withevent(B0,par,tspan,odeopt,varargin)
    % default
    noise = zeros(1,floor(tspan(end))+1); % no environmental noise
    fun = @odefunc;
    % arguments
    for j= 1:2:length(varargin)
        string= lower(varargin{j});
        switch string(1:min(3,length(string)))
            case 'noi'
                noise = varargin{j+1};
            case 'fun'
                fun = varargin{j+1};
            otherwise
                error('Unknown parameter name');
        end
    end
    
    % solve with extinction event
    te = tspan(1);
    res_t = [];
    res_b = [];
    par.extinSpe = zeros(par.nSpecies * par.nPatch,1);
    while te < tspan(end)
        odeopt = odeset(odeopt,'Events',@(t,B) extinEventFcn(t,B,par));
%         [t,b,te,be,ie] = ode15s_md(@(t,B) fun(B,par), [te(1) tmax], B0, odeopt);
        [t,b,te,be,ie] = ode15s_md(@(t,B) fun(B,par) + B .* noise(:,floor(t)+1), ...
            [te(1),tspan(find(tspan>te(1),1):end)], B0, odeopt);
        res_t = [res_t;t];
        res_b = [res_b;b];
        if isempty(te)
            break
        end
        B0 = be(1,:)';
        par.extinSpe(ie) = 1;
        B0(logical(par.extinSpe)) = 0;
        B0(B0 < 0) = 0;
    end
end