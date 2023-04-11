% only run with enough reps of raw abundance data in folder "rawdata",
% otherwise might produce results for very few rep instead of average over all

%% shared parameters
nSpecies = 50;
nPatch = 3;

% resDic{key} is a data metrix corresponding to key from metrixKey
% data metrix structure: lambda * heter * m * rep
metrixKey = {'alphaDiv', 'betaDiv', 'gammaDiv'};
metrixValue = cellfun(@(fn) res.(fn), metrixKey, 'uniformoutput', false);
metrixName = {'\alpha diversity', '\beta diversity', '\gamma diversity'};
resDic = containers.Map(metrixKey, metrixValue);
nameDic = containers.Map(metrixKey, metrixName);

lambList = [0,1,10,1e2,1e3,1e4]; % density dependency lambda
lamb_labels = {'0','10^0','10^1','10^2','10^3','10^4'};
heterList = [0:0.2:1]; % habitat heterogeneity H
mList = [0,10.^(-4:1/3:0)]; % migration rate m
if strcmp(snr, '')
    repList = 1:500;
else
    repList = 1:50;
end
nLam = numel(lambList);
nHet = numel(heterList);
nMig = numel(mList);
nRep = numel(repList);

% functions for calculating corresponding indexes
cutoff = 1e-6;  % global cutoff
cutoff_local = 1e-1;  % local cutoff
func_rescale = @(x_B) x_B ./ (max(x_B, [], 2) + eps);
func.gammaDiv = @(x_B) sum(max(x_B, [], 2) > cutoff);
func_mean_dist = @(dist) mean(dist(triu(true(size(dist)), 1)));
func.betaDiv = @(x_B) func_mean_dist(calcDissimilarity(x_B, 'Go'));
func.alphaDiv = @(x_B) mean(sum(func_rescale(x_B) > cutoff_local));

% load color configuration file
colorcollection;

%% fig S1, food web structure
figure()
load('randPars_s50p3r1000_0218062416.mat')
for ir = 1:6  %
    subplot(2, 3, ir)
    plotfoodweb(randPars{ir}.L')
    title(sprintf('rep%d', ir))
end
set(gcf, 'position', [506 178 960 543], 'PaperOrientation', 'landscape',...
    'renderer', 'Painters')
saveas(gcf, './figures/sup_fig1_foodweb_struc.pdf')

%% fig S2 (sup to fig1, emmigration rate vs. r_net for all lambda)
figure();
hold on
r_list = linspace(-0.5,0.5,1000);
nP = 3;
epsilon = 1e-5;
hline = [];
for il = 1:6
    lambda = lambList(il);
    m0 = 1 - 1 ./ (1 + (nP-1) * exp(lambda * (epsilon - r_list)));
    hline = [hline;plot(r_list, m0, 'color',[Blues(il,:) 0.9], 'linewidth',5)];
end
xlabel('r_{net} in orinigal patch', 'fontsize', 20)
ylabel('emigration rate', 'fontsize', 20)
set(gca, 'fontsize', 16, 'TickLabelInterpreter', 'latex', ...
    'ytick', [0, 2/3, 1], ...
    'ytickLabel', {'0', '$\frac{n_P-1}{n_P}m_0$', '$m_0$'})
hold off
hl = legend(hline,lamb_labels,'Location','SW');
title(hl,'\lambda')
hl.set('Fontsize',12)
saveas(gcf, sprintf('./figures/sup_emigration_rate.pdf'))

%% fig S3, net dispersal rate over time (till equilibrium)
lamb_labels = {'0','10^0','10^1','10^2','10^3','10^4'};
simres = cell(6,1);
figure()
colors = colormap(copper(5));
ir = 1;
load(sprintf('./rawdata/rep%03d_diff_H.mat',ir))
color_per_sp = colors(round(trophiclevel(parBase.L)), :);
set(gcf, 'DefaultAxesColorOrder',color_per_sp)
for ilam = 1:6
    ih = 4;
    istr = 1 + ilam;
    is = 7 * (ih - 1) + istr;

    B = BfinalSum{is};
    par = parBase;

    nLam = 6;
    nHet = 6;
    mList = repmat([0;ones(nLam,1)*0.1],1,nHet); % migration rate m0
    lambList = repmat([0, 0, 10, 1e2, 1e3, 1e4, 1e5]',1,nHet); % fitness dependency lambda
    heterList = repmat([0:0.2:1],1+nLam,1); % habitat heterogeneity H

    par.m = mList(is);
    par.lambda = lambList(is);
    par.heter = heterList(is);
    % apply heterogeneity. K * (1 + heter * K_mod)
    par.K = cellfun(@(dK) par.K .* (1 + par.heter * dK),...
        par.K_mod,'UniformOutput',false);
    par.r = cellfun(@(dr) par.r .* (1 + par.heter * par.isPlt .* dr),...
        par.r_mod,'UniformOutput',false);
    par.extinSpe = zeros(par.nSpecies * par.nPatch,1);

    odeopt = odeset('AbsTol',1e-11,'RelTol',1e-4,'NonNegative',...
        ones(parBase.nSpecies * parBase.nPatch,1));

    B0 = ones(parBase.nSpecies * parBase.nPatch,1) * 1; % initial biomass
    % B0 = B;
    [t,b] = ode15s_withevent(reshape(B0, [], 1),par,[0 2000],odeopt,'fun',@odefunc);

    m = zeros(size(b, 1), 50);
    for ii = 1:size(b, 1)
        [~, m_pop, m_new] = odefunc(b(ii, :), par);
        m_new(max(reshape(b(ii, :), parBase.nSpecies, parBase.nPatch), [], 2) < par.extin)...
            = nan;
        m(ii, :) = reshape(m_new, [], 1);
    end
    
%     tempres = simres{ilam};
%     t = tempres.t;
%     b = tempres.b;
%     m = tempres.m;
%     
    tempres.t = t;
    tempres.b = b;
    tempres.m = m;
    simres{ilam} = tempres;
    

    subplot(2,3,ilam)
    plot(t, m)
    ylim([0, 0.1])
    xlim([0, 1500])
    title(sprintf('$m_0=%s, \\lambda=%s$', num2str(par.m), lamb_labels{ilam}),...
        'Interpreter', 'latex')
    xlabel('Time (a.u.)', 'fontsize', 20)
    ylabel('Net emigration rate', 'fontsize', 20)
    set(gca, 'fontsize', 16)
end

set(gcf, 'Position', [245 116 1119 594],...
   'PaperSize', [11.8 6.3]) 
hc = colorbar('Ticks', 0.1:0.2:1, 'TickLabels', 1:5,...
    'Position', [0.93, 0.2, 0.0147, 0.6118]);
set(get(hc,'label'),'String','Trophic Level')

% saveas(gcf, './figures/sup_fig3_net_mig_rate.svg')

%% fig ST3.1 - gamma different cutoff
%% fig ST3.2 - gamma diversity vs. cutoff
cutoff_list = 10 .^ linspace(-10, 0, 51);
gamma_array = nan(length(cutoff_list), 42, 500);
for ir = 1:500
    try
        load(sprintf('rawdata/rep%03d_diff_H.mat', ir))
    catch ME
        % some reps might be missing due to integration error
        if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
            continue
        else
            rethrow(ME)
        end
    end
    for istr = 1:42
        if isempty(BfinalSum{istr})  % integration failure
            gamma_array(:, :, ir) = nan;
            break
        end
        B = BfinalSum{istr};
        for ic = 1:length(cutoff_list)
            cutoff = cutoff_list(ic);
            gamma_array(ic, istr, ir) = sum(max(B, [], 2) > cutoff);
        end
    end
end

ih = 4;
colors = Blues;

% div vs cutoff
figure()
hold on
hline = [];
for il = 1:7
    istr = (ih - 1) * 7 + il;
    if il == 1
        hline = [plot(log10(cutoff_list),nanmean(gamma_array(:, istr, :), 3),...
            '-','linewidth',10,'color',[0.8094, 0.9251, 0.7826], 'markersize',30);hline];
    else
        hline = [hline;plot(log10(cutoff_list),nanmean(gamma_array(:, istr, :), 3),...
            '-','linewidth',3,'color',[colors(il,:) 0.6], 'markersize',10)];
    end
end
hold off

hl = legend(hline,["iso", '0','10^0','10^1','10^2','10^3','10^4'],'Location','SW');
title(hl,'\lambda')
hl.set('Fontsize',12)
ylim([20, 40])
set(gca,'fontsize',16)
xlabel('Global cutoff','fontsize',20)
ylabel('\gamma diversity','fontsize',20)
xticks = -10:0;
set(gca, 'xtick', xticks,...
    'xtickLabel', cellfun(@(lgx) sprintf('10^{%d}', lgx), num2cell(xticks), 'UniformOutput', false),...
    'fontsize',16);
saveas(gcf,sprintf('./figures/sup_text3_gammaDiv_vs_cutoff.svg'))

% div vs. lambda, different cutoff
figure()
errbarOn = true;
plotcutoffs = 10 .^ (-10:2:0); % cutoff values to plot div vs lambda
cutoff_idx = []; % indexes of these cutoff values to be plot
for c = plotcutoffs
    cutoff_idx = [cutoff_idx, find(cutoff_list == c)];
end
for ic = 1:length(plotcutoffs) % loop through three diversity metrics
    subplot(2,3,ic)
    dataKey = metrixKey{3};
    
    rawdata = reshape(squeeze(gamma_array(cutoff_idx(ic), :, :)), 7, 6, []);
    rawdata = rawdata([2:7, 1], :, :);
    dataCell = mat2cell(rawdata,ones(nLam+1,1),ones(nHet,1),nRep);
    relDataCell = mat2cell(rawdata - rawdata(end,:,:),ones(nLam+1,1),ones(nHet,1),nRep);
    % mean of diversity, including isolated scenario
    data = cellfun(@(x) mean(x,'omitnan'),dataCell);
%     % mean of diversity relative to isolated scenario
    reldata = cellfun(@(x) mean(x,'omitnan'),relDataCell);
%     % bootstrap 95% confidence interval as error bar
%     dataerr_ci = cellfun(@(x) bootci(100, @nanmean, x),...
%         relDataCell, 'UniformOutput', false);
    % 1/4 quantile as error bar
    dataerr_ci = cellfun(@(x) [prctile(x, 25), prctile(x, 75)], ...
        relDataCell, 'UniformOutput', false);
    dataerr_neg = reldata - cellfun(@(x) x(1), dataerr_ci);
    dataerr_pos = cellfun(@(x) x(2), dataerr_ci) - reldata;
    colors = Oranges(2:end,:);

    hold on
    hline = [];
    for ih = 1:nHet
        % solid line linking different lambda
        hline = [hline;plot(1:nLam,data(1:nLam,ih),...
            '.-','linewidth',3,'color',[colors(ih,:) 0.6], 'markersize',10)];
        % dashed line linking biggest lambda to isolated
        plot(nLam:nLam+1,data(nLam:end,ih),...
            '.--','linewidth',1,'color',[colors(ih,:) 0.6],'markersize',10);
        if errbarOn
            errorbar(1:nLam+1,data(:,ih),dataerr_neg(:,ih),dataerr_pos(:,ih),'.',...
                'color',[colors(ih,:) 0.6],'linewidth',1.5);
        end
    end
    hold off
    set(gca,'xtick',1:nLam+1,...
        'xtickLabel',[lamb_labels,{'iso'}],...
        'fontsize',16);
    xlabel('Fitness dependency of dispersal (\lambda)', 'fontsize',20)
    ylabel(nameDic(dataKey), 'fontsize',20)
    title(sprintf('global cutoff = 10^{%d}', round(log10(plotcutoffs(ic)))))
%     if errbarOn
%         saveas(gcf,sprintf('./%s/fig2_diff_H_%s_witherr.pdf',figpath,dataKey))
%     else
%         saveas(gcf,sprintf('./%s/fig2_diff_H_%s.pdf',figpath,dataKey))
%     end
end
set(gcf, 'Units', 'inches', 'Position', [3.1053, 1.0421, 19.3053, 8.6316],...
    'PaperSize', [19.3053, 8.6316], 'PaperUnits', 'inches')
saveas(gcf, './figures/sup_text3_fig_gamma_diff_cutoff.svg')
%% fig ST3.4 alpha diversity vs. cutoff
cutoff_list = 10 .^ linspace(-5, 0, 51);
alpha_array = nan(length(cutoff_list), 42, 500);
for ir = 1:500
    try
        load(sprintf('rawdata/rep%03d_diff_H.mat', ir))
    catch ME
        % some reps might be missing due to integration error
        if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
            continue
        else
            rethrow(ME)
        end
    end
    for istr = 1:42
        if isempty(BfinalSum{istr})  % integration failure
            alpha_array(:, :, ir) = nan;
            break
        end
        B = BfinalSum{istr};
        B = B ./ (max(B, [], 2) + eps);
        for ic = 1:length(cutoff_list)
            cutoff = cutoff_list(ic);
            alpha_array(ic, istr, ir) = mean(sum(B > cutoff));
        end
%         plot(log10(cutoff_list), alpha_list)
    end
end

ih = 4;
colors = Blues;

figure()
hold on
hline = [];
for il = 1:7
    istr = (ih - 1) * 7 + il;
    if il == 1
        hline = [plot(log10(cutoff_list),nanmean(alpha_array(:, istr, :), 3),...
            '-','linewidth',10,'color',[0.8094, 0.9251, 0.7826], 'markersize',30);hline];
    else
        hline = [hline;plot(log10(cutoff_list),nanmean(alpha_array(:, istr, :), 3),...
            '-','linewidth',3,'color',[colors(il,:) 0.6], 'markersize',10)];
    end
end
hold off

hl = legend(hline,["iso", '0','10^0','10^1','10^2','10^3','10^4'],'Location','SW');
title(hl,'\lambda')
hl.set('Fontsize',12)
ylim([20, 40])
set(gca,'fontsize',16)
xlabel('Local cutoff','fontsize',20)
ylabel('\alpha diversity','fontsize',20)
xticks = -6:0;
set(gca, 'xtick', xticks,...
    'xtickLabel', cellfun(@(lgx) sprintf('10^{%d}', lgx), num2cell(xticks), 'UniformOutput', false),...
    'fontsize',16);
saveas(gcf,sprintf('./figures/sup_text3_alphaDiv_vs_cutoff.svg'))

%% div vs. alpha, different cutoff
figure()
errbarOn = true;
plotcutoffs = [10 .^ (-5:-1), cutoff_list(46)]; % cutoff values to plot div vs lambda
cutoff_idx = []; % indexes of these cutoff values to be plot
for c = plotcutoffs
    cutoff_idx = [cutoff_idx, find(cutoff_list == c)];
end
cutoff_labels = {'0.001%', '0.01%', '0.1%', '1%', '10%', '30%'};
for ic = 1:length(plotcutoffs) % loop through three diversity metrics
    subplot(2,3,ic)
    dataKey = metrixKey{1};
    
    rawdata = reshape(squeeze(alpha_array(cutoff_idx(ic), :, :)), 7, 6, []);
    rawdata = rawdata([2:7, 1], :, :);
    dataCell = mat2cell(rawdata,ones(nLam+1,1),ones(nHet,1),nRep);
    relDataCell = mat2cell(rawdata - rawdata(end,:,:),ones(nLam+1,1),ones(nHet,1),nRep);
    % mean of diversity, including isolated scenario
    data = cellfun(@(x) mean(x,'omitnan'),dataCell);
%     % mean of diversity relative to isolated scenario
    reldata = cellfun(@(x) mean(x,'omitnan'),relDataCell);
%     % bootstrap 95% confidence interval as error bar
%     dataerr_ci = cellfun(@(x) bootci(100, @nanmean, x),...
%         relDataCell, 'UniformOutput', false);
    % 1/4 quantile as error bar
    dataerr_ci = cellfun(@(x) [prctile(x, 25), prctile(x, 75)], ...
        relDataCell, 'UniformOutput', false);
    dataerr_neg = reldata - cellfun(@(x) x(1), dataerr_ci);
    dataerr_pos = cellfun(@(x) x(2), dataerr_ci) - reldata;
    colors = Oranges(2:end,:);

    hold on
    hline = [];
    for ih = 1:nHet
        % solid line linking different lambda
        hline = [hline;plot(1:nLam,data(1:nLam,ih),...
            '.-','linewidth',3,'color',[colors(ih,:) 0.6], 'markersize',10)];
        % dashed line linking biggest lambda to isolated
        plot(nLam:nLam+1,data(nLam:end,ih),...
            '.--','linewidth',1,'color',[colors(ih,:) 0.6],'markersize',10);
        if errbarOn
            errorbar(1:nLam+1,data(:,ih),dataerr_neg(:,ih),dataerr_pos(:,ih),'.',...
                'color',[colors(ih,:) 0.6],'linewidth',1.5);
        end
    end
    hold off
    set(gca,'xtick',1:nLam+1,...
        'xtickLabel',[lamb_labels,{'iso'}],...
        'fontsize',16);
    xlabel('Fitness dependency of dispersal (\lambda)', 'fontsize',20)
    ylabel(nameDic(dataKey), 'fontsize',20)
    title(sprintf('local cutoff = %s', cutoff_labels{ic}))
%     if errbarOn
%         saveas(gcf,sprintf('./%s/fig2_diff_H_%s_witherr.pdf',figpath,dataKey))
%     else
%         saveas(gcf,sprintf('./%s/fig2_diff_H_%s.pdf',figpath,dataKey))
%     end
end
set(gcf, 'Units', 'inches', 'Position', [3.1053, 1.0421, 19.3053, 8.6316],...
    'PaperSize', [19.3053, 8.6316], 'PaperUnits', 'inches')
saveas(gcf, './figures/sup_text3_fig_alpha_diff_cutoff.svg')

%% fig ST3. - species abundance distribution (local and global) different strategy
nRep = 500;
nSpecies = 50;
nPatch = 3;
abundance = nan(42, nRep, nSpecies, nPatch);
snr = '';
extin = 1e-10;
for ir = 1:nRep
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
    if size(BfinalSum{end}) == 0
        continue
    end
    for is = 1:42
        BfinalSum{is} = BfinalSum{is} .* (BfinalSum{is} > 0);
        relabundance(is, ir, :, :) = BfinalSum{is} ./ (max(BfinalSum{is}, [], 2) + eps);
        abundance(is, ir, :, :) = BfinalSum{is};
        maxabundance(is, ir, :, :) = (max(BfinalSum{is}, [], 2) + eps);
        sumabundance(is, ir, :, :) = sum(BfinalSum{is}, 2);
    end
end

lambList = [0,1,10,1e2,1e3,1e4]; % density dependency lambda
lamb_labels = {'0','10^0','10^1','10^2','10^3','10^4'};
heterList = [0:0.2:1]; % habitat heterogeneity H
mList = [0,10.^(-4:1/3:0)]; % migration rate m

% relative abundance, separate H and lambda
figure()
ih = 2;
il_use = [1, 4, 6];
for ii = 1:3
    subplot(3,1,ii)
    hold on
    il = il_use(ii);
%     %isolate
%     histogram(log10(reshape(relabundance((ih-1) * 7 + 1 + 0, :, :, :), [], 1)), ...
%         linspace(-4,1,20), 'Normalization', 'prob',...
%         'FaceColor', Oranges(5, :))
    %connected
    for is = (ih-1) * 7 + 1 + [il]
        histogram(log10(reshape(relabundance(is, :, :, :), [], 1)), ...
            linspace(-4,1,20), 'Normalization', 'prob',...
            'FaceColor', Blues(il+1, :))
    end
    legend(sprintf('lambda=%s', lamb_labels{il}), 'location', 'NW')
    if ii == 3
        xlabel('log10(relative population density)')
    end
    ylabel('Probability')
    xlim([-4, 0.5])
    ylim([0, 0.7])
end
set(gcf, 'Position', [744 614 493 420])
suptitle(sprintf('H=%s', string(heterList(ih))))
saveas(gcf, sprintf('./figures/sup_text3_abun_distr_rel_H=%.1f.svg', heterList(ih)))

% absolute abundance distribution for isolated species, across H
figure()
histogram(log10(abundance(1:7:42, :, :, :)), linspace(-10,1,20), 'Normalization', 'Probability')
xlabel('log10(population density)', 'fontsize', 20)
ylabel('Probability', 'fontsize', 20)
set(gca,'fontsize',16)
xlim([-10, 2])
saveas(gcf, './figures/sup_text3_abun_distr_abs.svg')

% absolute abundance, separate H and lambda
figure()
ih = 6;
il_use = [1, 4, 6];
for ii = 1:3
    subplot(3,1,ii)
    hold on
    il = il_use(ii);
    %isolate
    histogram(log10(reshape(maxabundance((ih-1) * 7 + 1 + 0, :, :, :), [], 1)), ...
        linspace(-4,2,20), 'Normalization', 'prob',...
        'FaceColor', Oranges(5, :))
    %connected
    for is = (ih-1) * 7 + 1 + [il]
        histogram(log10(reshape(maxabundance(is, :, :, :), [], 1)), ...
            linspace(-4,2,20), 'Normalization', 'prob',...
            'FaceColor', Blues(il+1, :))
    end
    legend('iso', sprintf('lambda=%s', lamb_labels{il}), 'location', 'NW')
    if ii == 3
        xlabel('log10(maximum population density)')
    end
    ylabel('Probability')
    xlim([-4, 2])
%     ylim([0, 0.2])
end
set(gcf, 'Position', [744 614 493 420])
suptitle(sprintf('H=%s', string(heterList(ih))))
saveas(gcf, sprintf('./figures/sup_text3_abun_distr_abs_H=%.1f.svg', heterList(ih)))