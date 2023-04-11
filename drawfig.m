function drawfig(snr)
%% load calculated diversity data
load(sprintf('%sres.mat',snr))

% path to save figures for each plant dispersal scenario
figpath = sprintf('%sfigures',snr);
mkdir(figpath)

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

%% fig.2,a-c - show diversity-strategy relationship for all replicates
ih = 4;
im = 11;
useReps = 1:min(100, nRep);
for ik = 1
    dataKey = metrixKey{ik};
    data = resDic(dataKey);
    data = squeeze(cat(1,data(:, ih, im, useReps),data(1, ih, 1, useReps))); % (nLam + 1) * nRep
    reldata = data - data(7, :);
    % interpolate color to represent different migraiton rates
    Greens_mig = interp1((nMig+1)*(0:6),Greens,6*(1:(nMig+1)),'spline');
    colors = Greens_mig(2:end,:);
    figure
    hline = [];
    hold on
    hline = [hline; plot(1:nLam+1, zeros(nLam+1,1), '.-','linewidth',10,'color',[0.8094, 0.9251, 0.7826],...
        'markersize',20,'markeredgecolor','w')];
    plot([1:nLam+1]' + rand(size(data)) * 0.4 - 0.2, reldata,...
        'color', [0.1, 0.1, 0.1, 0.05], 'linewidth', 2)
    hline = [hline; plot(1:nLam, nanmean(reldata(1:nLam, :), 2),...
        '.-','linewidth',5,'color',[colors(im,:), 0.6], 'markersize',20)];
    plot(nLam:nLam+1,nanmean(reldata(nLam:end, :), 2),...
        '.--','linewidth',2,'color',[colors(im,:)],'markersize',20);
    % dataerr_ci = cellfun(@(x) bootci(100, {@nanmean, x}, 'alpha', 0.01),...
    %     mat2cell(reldata, ones(nLam+1, 1), [nRep]), 'UniformOutput', false);
    dataerr_ci = cellfun(@(x) [prctile(x, 25), prctile(x, 75)], ...
        mat2cell(reldata, ones(nLam+1, 1), [length(useReps)]), 'UniformOutput', false);
    dataerr_neg = nanmean(reldata, 2) - cellfun(@(x) x(1), dataerr_ci);
    dataerr_pos = cellfun(@(x) x(2), dataerr_ci) - nanmean(reldata, 2);
    errorbar(1:nLam+1,nanmean(reldata, 2),dataerr_neg,dataerr_pos,'.',...
        'color',[colors(im,:) 0.6],'linewidth',1.5);

    set(gca,'xtick',1:nLam+1,...
        'xtickLabel',[lamb_labels, {'iso'}],...
        'fontsize',16);
    xlabel('Fitness dependency of dispersal (\lambda)', 'fontsize',20)
    ylabel(sprintf('Relative %s', nameDic(dataKey)), 'fontsize',20)
    if ik == 3
        ylim([-5, 1])
    end
    if ik == 3
        hl = legend(hline,["m_0=0 (iso)", "m_0=0.1"],'Position',[0.7513 0.4533 0.2449 0.1368]);
        hl.set('Fontsize',12)
    end
    saveas(gcf,sprintf('./%s/fig2_all_reps_%s_witherr.pdf',figpath,dataKey))
end

%% fig 2,d-f - diversity vs. lambda (under different heterogeneity)
errbarOn = true;
for ik = 1:3 % loop through three diversity metrics
    dataKey = metrixKey{ik};
    data = resDic(dataKey);
    
    rawdata = squeeze(cat(1,data(:, :, 11, :),data(1, :, 1, :)));
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

    figure
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
    switch dataKey % show legend in suitable locations
        case {'gammaDiv'}
            hl = legend(hline,string(heterList),'Position',[0.8584 0.3647 0.1368 0.3417]);           
            title(hl,'H')
            hl.set('Fontsize',12)
        otherwise
    end
    set(gca,'xtick',1:nLam+1,...
        'xtickLabel',[lamb_labels,{'iso'}],...
        'fontsize',16);
    xlabel('Fitness dependency of dispersal (\lambda)', 'fontsize',20)
    ylabel(nameDic(dataKey), 'fontsize',20)
    if errbarOn
        saveas(gcf,sprintf('./%s/fig2_diff_H_%s_witherr.svg',figpath,dataKey))
    else
        saveas(gcf,sprintf('./%s/fig2_diff_H_%s.pdf',figpath,dataKey))
    end
end
%% fig 2, g-i - diversity vs. lambda (under different migration rate)
errbarOn = false;
for ik = 1:3 % loop through three diversity metrics
    dataKey = metrixKey{ik};

    rawdata = resDic(dataKey);
    rawdata = squeeze(rawdata(:, 4, :, :));
    reldata = rawdata - rawdata(1,1,:);
    relDataCell = mat2cell(reldata,ones(nLam,1),ones(nMig,1),nRep);
    data = mean(rawdata, 3, 'omitnan');
%     %     % bootstrap 95% confidence interval as error bar
%     dataerr_ci = cellfun(@(x) bootci(100, @nanmean, x),...
%         relDataCell, 'UniformOutput', false);
    % 1/4 quantile as error bar
    dataerr_ci = cellfun(@(x) [prctile(x, 25), prctile(x, 75)], ...
        relDataCell, 'UniformOutput', false);
    dataerr_neg = reldata - cellfun(@(x) x(1), dataerr_ci);
    dataerr_pos = cellfun(@(x) x(2), dataerr_ci) - reldata;
    % interpolate color to represent different migraiton rates
    Greens_mig = interp1((nMig+1)*(0:6),Greens,6*(1:(nMig+1)),'spline');
    colors = Greens_mig(2:end,:);

    figure
    hold on
    hline = [];
    for im = 2:nMig
        % solid line linking different lambda
        hline = [hline;plot(1:nLam,data(1:nLam,im),...
            '.-','linewidth',3,'color',[colors(im,:) 0.6], 'markersize',10)];
        if errbarOn
%             errorbar(1:nLam,data(1:nLam,im),dataerr(1:nLam,im),'.',...
%                 'color',[colors(im,:) 0.6], 'linewidth',1.5);
            errorbar(1:nLam,data(1:nLam,im),dataerr_neg(1:nLam,im),dataerr_pos(1:nLam,im),'.',...
                'color',[colors(im,:) 0.6], 'linewidth',1.5);
        end
    end
    hline = [plot(1:nLam,data(1:nLam,1),...
        '.-','linewidth',10,'color',[colors(1,:)], 'markersize',30,...
        'marker', '.', 'markeredgecolor','w');hline];
    hold off
    switch dataKey % show legend in different locations
        case {'gammaDiv'}
            hl = legend(hline([1,2:3:nMig]),string(mList([1,2:3:nMig])),...
                'Position',[0.8070 0.3533 0.1870 0.3644]);
            title(hl,'m_0')
            hl.set('Fontsize',12)
        otherwise
    end
    set(gca,'xtick',1:nLam,'xtickLabel',lamb_labels,'fontsize',16);
    xlabel('Fitness dependency of dispersal (\lambda)','fontsize',20)
    ylabel(nameDic(dataKey),'fontsize',20)
    if errbarOn
        saveas(gcf,sprintf('./%s/fig2_diff_m0_%s_witherr.pdf',figpath,dataKey))
    else
        saveas(gcf,sprintf('./%s/fig2_diff_m0_%s.pdf',figpath,dataKey))
    end
end

%% fig 3 - diversity vs. heterogeneity (under differnet lambda)
errbarOn = true;
for ik = 1:3 % loop through three diversity metrics
    dataKey = metrixKey{ik};

    rawdata = resDic(dataKey);
    rawdata = squeeze(cat(1,rawdata(:, :, 11, :),rawdata(1, :, 1, :)));
    relDataCell = mat2cell(rawdata - rawdata(end,:,:),ones(nLam+1,1),ones(nHet,1),nRep);
    data = nanmean(rawdata, 3);
    reldata = cellfun(@(x) mean(x,'omitnan'),relDataCell);
    %     %     % bootstrap 95% confidence interval as error bar
    dataerr_ci = cellfun(@(x) bootci(100, @nanmean, x),...
        relDataCell, 'UniformOutput', false);
    % 1/4 quantile as error bar
    dataerr_ci = cellfun(@(x) [prctile(x, 25), prctile(x, 75)], ...
        relDataCell, 'UniformOutput', false);
    dataerr_neg = reldata - cellfun(@(x) x(1), dataerr_ci);
    dataerr_pos = cellfun(@(x) x(2), dataerr_ci) - reldata;
    colors = Blues;

    figure
    hold on
    hline = [];
    hline = [plot(heterList,data(end,:),...
        '.-','linewidth',10,'color',[0.8094, 0.9251, 0.7826], 'markersize',30,...
        'marker', '.', 'markeredgecolor','w');hline];
    for il = 1:nLam
        % solid line linking different lambda
        hline = [hline;plot(heterList,data(il, :),...
            '.-','linewidth',3,'color',[colors(il+1,:) 0.6], 'markersize',10)];
        if errbarOn
            errorbar(heterList,data(il,:),dataerr_neg(il,:),dataerr_pos(il,:),'.',...
                'color',[colors(il+1,:) 0.6],'linewidth',1.5);
        end
    end
    hold off
    switch dataKey % show legend in different locations
        case {'betaDiv'}
            hl = legend(hline,["iso", lamb_labels],'Position',[0.1705 0.4935 0.1414 0.4552]);
            title(hl,'\lambda')
            hl.set('Fontsize',12)
        otherwise
    end
    set(gca,'fontsize',16)
    xlabel('Habitat heterogeneity (H)','fontsize',20)
    ylabel(nameDic(dataKey),'fontsize',20)
    if errbarOn
        saveas(gcf,sprintf('./%s/fig3_div_heter_%s_witherr.pdf',figpath,dataKey))
    else
        saveas(gcf,sprintf('./%s/fig3_div_heter_%s.pdf',figpath,dataKey))
    end
end
%% fig 4, diversity vs. migration rate (under differnet lambda)
errbarOn = false;

if errbarOn % when error bar is on plotting all lambda makes the figure too crowded
    lambListSlice = [1 4 6];
else
    lambListSlice = [1:nLam];
end
for ik = 1:3 % loop through three diversity metrics
    dataKey = metrixKey{ik};

    data = resDic(dataKey);
    dataerr = std(data(:, 4, :, :) - data(1, 4, 1, :), 1, 4, 'omitnan');
    data = squeeze(mean(data(:, 4, :, :), 4, 'omitnan'));
    colors = Blues;

    figure
    hold on
    hline = [];
    for il = lambListSlice
        % solid line linking different lambda
        hline = [hline;plot([-5, log10(mList(2:end))],data(il, :),...
            '.-','linewidth',3,'color',[colors(il+1,:) 0.6], 'markersize',10)];
        if errbarOn
            errorbar([-5, log10(mList(2:end))],data(il,:),dataerr(il,:),'.',...
                'color',[colors(il+1,:) 0.6], 'linewidth',1.5);
        end
    end
    hold off
    switch dataKey % show legend in different locations
        case {'betaDiv'}
            hl = legend(hline,lamb_labels(lambListSlice),'Location','SW');
            title(hl,'\lambda')
            hl.set('Fontsize',12)
        otherwise
    end
	set(gca,'xtick',[-5:0],'xtickLabel',['0',string(mList(2:3:nMig))], 'fontsize',16);
    xlabel('Habitat connectivity (m_0)','fontsize',20)
    ylabel(nameDic(dataKey),'fontsize',20)
    if errbarOn
        saveas(gcf,sprintf('./%s/fig4_div_m0_%s_witherr.pdf',figpath,dataKey))
    else
        saveas(gcf,sprintf('./%s/fig4_div_m0_%s.svg',figpath,dataKey))
    end
end
%% below are figures independent of plant dispersal scenario
if ~strcmp(snr, '')
    return
end
%% fig1 - emmigration rate as a function of r
figure();
epsilon = 1e-5;
r_list = linspace(-0.5, 0.5, 1000);
funcs = cell(3,1);
lambdaListSlice = [1, 3, 5];
rates = [0.15, -0.15, 0.45];
colors = [green(4,:);oranRed(4,:);blue(4,:)];
x = linspace(-0.1-epsilon,0.1-epsilon,1000);
nP = 3;
for il = lambdaListSlice
    subplot(2, 3, il);
    hold on
    lambda = lambList(il);
    m0 = 1 - 1 ./ (1 + (nP-1) * exp(lambda * (epsilon - r_list)));
    hline = plot(r_list, m0, 'color',[Blues(il,:) 0.5], 'linewidth',10);
%     xlabel('$r_k^{net}-r_1^{net}-\epsilon\cdot\delta_{k=1}$',...
%         'fontsize',14,'interpreter','latex')
    for ir = 1:3
        plot([1, 1] * (rates(ir)), [0, 1], '--', ...
            'color', [colors(ir,:), 0.5], 'linewidth', 3)
%         temp_x = rates(ir) - rates(1) - epsilon * (ir ~= 1);
%         plot(temp_x, funcs{ilam}(temp_x),'.','markersize',80,...
%             'markerfacecolor', colors(ir,:), 'markeredgecolor',colors(ir,:))
    end
    xlabel('r_{net} in natal patch', 'fontsize', 20)
    ylabel('emmigration rate', 'fontsize', 20)
    set(gca, 'fontsize', 14, 'TickLabelInterpreter', 'latex', ...
        'ytick', [0, 1/3, 2/3, 1], ...
        'ytickLabel', {'0', '$\frac{1}{3}m_0$', '$\frac{2}{3}m_0$', '$m_0$'})
    hold off
    xlim([-0.5, 0.5])
    legend([hline], strcat('\lambda=', lamb_labels(il)), 'location', 'northoutside')
end
set(gcf, 'Position', [169, -36, 1030, 567], 'PaperOrientation', 'landscape');
saveas(gcf, sprintf('./figures/fig1.pdf'))

%% fig 5(a) - example community structures
ir = 2; % rep
ih = 6; % heterogeneity
snr = '';
% higher heterogeneity makes the difference more obvious but the trend should be the same
load(sprintf('./rawdata/%srep%03d_diff_H.mat',snr,ir))
% B_example_array{1} = BfinalSum{7*3+1};  % isolate
B_example_array = BfinalSum(7*(ih-1)+1:7*ih);
% load(sprintf('./rawdata/%srep%03d_diff_m0.mat',snr,ir))
% B_example_array(2:7) = BfinalSum(6*12+(1:6));

figure('Resize', 'off')
nPatch = 3;
nSpecies = 50;

for is = 1:7
%     B = BfinalSum{7 * (ih - 1) + is};
    B = B_example_array{is};
    imgAlpha = ones(nSpecies, nPatch);
    imgAlpha(B < 1e-6) = 0;
    
    if is == 1
        hList(7) = subplot(1,8,7);
    else
        hList(is-1) = subplot(1,8,is-1);
    end
    
    colormap(viridis);
    im = imagesc(log10(B+eps), 'AlphaData', imgAlpha, [-4,0]);
    set(gca, 'color', 'k')
%     im = imagesc(B, 'AlphaData', imgAlpha, [0,1]);
%     title('Local population density on log scale')
%     xlabel('Patch')
%     ylabel('Species')
    set(gca,'xtick',1:nPatch, 'fontsize',12);
%     set(gca,'ytick',1:nSpecies);
    set(gca,'YDir','normal')
%     colorbar;
    if is == 2
        ylabel('Species','Fontsize',20)
%         set(gca,'ytick',[]);
    end
    if is ~= 2
        set(gca,'yticklabels',[]);
    end
    if is - 1 == 4
        xlabel('Patches','Fontsize',20)
    end
    if is == 1
        title('iso');
    else
        title(lamb_labels{is-1});
    end
end
    
hList(8) = subplot(1,8,8);
axis off
hc = colorbar('Ticks', linspace(0,1,5),...
    'TickLabels',cellfun(@(x) sprintf('10^{%s}',x),cellstr(num2str([-4:0]')),'Uni',false),...
    'fontsize',9);
% hc = colorbar('fontsize',9);
set(hc,'Location','west')
title(hc,{'Population','density'})

pause(1)

posCell = cell(1,8);
for is = 1:8
    posCell{is} = hList(is).Position;
end

for is = 1:8
    pos = posCell{is};
    if is > 1
        pos(1) = pos(1) - 0.025 * (is - 1);
        if is == 7
            pos(1) = pos(1) + 0.02;
        elseif is == 8
            pos(1) = pos(1) + 0.05;
        end
    end
    hList(is).Position = pos;
end

set(gcf,'Color',[1 1 1]); set(gca,'Color',[.8 .8 .8]); set(gcf,'InvertHardCopy','off');

saveas(gcf, sprintf('./figures/fig5_smp_comm.pdf'))
%% fig 5(b)
% lambda * heter * m * rep
figure
set(gcf, 'Position', [437.0000   45.0000  300  300.0000]);
for ik = 1:3 % loop through three diversity metrics
    dataKey = metrixKey{ik};
    data = cellfun(@(x_B) func.(dataKey)(x_B), B_example_array);
    subplot(3,1,ik)
    plot(1:6, data(2:7), '.-', 'color', [0.3 0.3 0.3 0.6], 'linewidth',3);
    hold on
    plot(6:7, data([7,1]), '.--', 'color', [0.3 0.3 0.3 0.6], 'linewidth',1);
    hold off
    if ik == 3
        set(gca,'xtick',1:nLam+1,...
            'xtickLabel',[lamb_labels,'iso'],...
            'fontsize',9);
        xlabel('Fitness dependency of dispersal (\lambda)', 'fontsize',12)
    end
    ylabel(nameDic(dataKey), 'fontsize',12)
    if ik < 3
        set(gca,'XTick',[]);
    end
    ylim([min(data) - 0.1 * (max(data)-min(data)),
        min(data) + 1.1 * (max(data)-min(data))])
    xlim([0.8, 8])
end
saveas(gcf, sprintf('./figures/fig5_smp_comm_divplot.pdf'))

%% fig S1 - heatmap of diversity/stability under 2D conditions: lambda * heterogeneity
withiso = true;
for ik = 1:3 % loop through three diversity metrics
    dataKey = metrixKey{ik};
    data = resDic(dataKey);
    if withiso
        data = [squeeze(mean(data(:, :, 11, :), 4, 'omitnan'));
            mean(data(1, :, 1, :), 4, 'omitnan')];
    else
        data = squeeze(mean(data(:, :, 11, :), 4, 'omitnan'));
    end

    figure
    imagesc(data')
    
    xlabel('Fitness dependency of dispersal (\lambda)','fontsize',20)
    ylabel('Habitat heterogeneity (H)','fontsize',20)
    title(nameDic(dataKey),'fontsize',20)
    set(gca,'ytickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},...
        'fontsize',12);
    colorbar
    if withiso
        set(gca,'xtickLabel',[lamb_labels, 'iso'],...
            'fontsize',12);
        hold on
        plot([6.5,6.5],[0.5,6.5],'k-')
        hold off
        saveas(gcf,sprintf('./%s/figs1_%s_heatmap_lh_withiso.svg',figpath,dataKey))
    else
        set(gca,'xtickLabel',lamb_labels);
        saveas(gcf,sprintf('./%s/figs4_%s_heatmap_lh.pdf',figpath,dataKey))
    end
end
%% fig S2 - heatmap of diversity under 2D conditions: lambda * m0
withiso = false;
for ik = 1:3 % loop through three diversity metrics
    dataKey = metrixKey{ik};
    data = resDic(dataKey);
    data = squeeze(mean(data(:, 4, :, :), 4, 'omitnan'));

    figure
    imagesc(data')
    
    xlabel('Fitness dependency of dispersal (\lambda)','fontsize',20)
    ylabel('Habitat connectivity (m_0)','fontsize',20)
    title(nameDic(dataKey),'fontsize',20)
    
    yticks = cellfun(@(x) sprintf('$10^{-%d\\frac{%d}{3}}$',fix(x/3),mod(x,3)),...
        num2cell(12:-1:0),'UniformOutput',false);
    yticks(1:3:end-1) = cellfun(@(x) sprintf('$10^{-%d}$',x),...
        num2cell(4:-1:1),'UniformOutput',false);
    yticks{end} = '$10^0$';
    yticks = {'0',yticks{:}};
    set(gca,'ytick',1:14,'ytickLabel',yticks,'fontsize',12,...
        'xtick',1:nLam,'xtickLabel',cellfun(@(lb) ['$',lb,'$'], lamb_labels, 'uni', false),...
        'TickLabelInterpreter','latex');
    
    colorbar
    saveas(gcf,sprintf('./%s/figs5_%s_heatmap_lm.svg',figpath,dataKey))
end
end