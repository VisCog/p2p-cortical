%temporal_fit_Winawer_data.m
%
% fits Winawer bightness data, corresponds to figure 2c

clear all; close all

T = readtable('datasets/Winawer2016_data.xlsx');
T.brightness( T.brightness==-1000) = NaN; % weird hack because otherwise brightness stuck as a cell array

tp = p2p_c.define_temporalparameters();
tp.model = 'compression';
tp.simdur = 3;
rng(11171964) % fix the random number generator, used Geoff's birthday in paper

gvals = find(T.amp>0); % only include trials where there was a simulus and a reported brightness
gvals = gvals(randperm(length(gvals)));
rsamp(1).val = gvals(1:ceil(length(gvals)/2));
rsamp(2).val  = gvals(ceil(length(gvals)/2)+1:end);

titleStr = {'Training Data', 'Test Data'};
FITFLAG = 1;
if FITFLAG
    Texp = T(rsamp(1).val,:);
    freeParams = {'power'};
    tp = fit('p2p_c.fit_brightness',tp,freeParams,Texp);
end
colorList  = [1 0 0; 0 1 0; 1 .7 0; 1 .3  1; 0 0 1 ]; % roughly match colors to Winawer paper
alpha = 0.5; sz = 150; sd = .1; % jitter for plotting

for site =1:5
    % plot model predictions vs. data
    for r = 1:2 % test and train data
        eid = intersect(rsamp(r).val, find(T.electrode==site));
        if ~isempty(eid)
            Texp = T(eid, :);
            trl= p2p_c.loop_convolve_model(tp, Texp);
            subplot(1, 3, r)
            x = [trl.maxresp]; y = [Texp.brightness] ;
            x = reshape(x, length(x), 1); y = reshape(y, length(x), 1);
            x = x + sd*(rand(size(x))-0.5); y = y + sd*(rand(size(x))-0.5);
            h = scatter(x, y, sz, 'MarkerFaceColor', colorList(site, :), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha);
            hold on
            t = text(15, site/2, ['corr = ', num2str(round(corr(x, y), 3))]);
            set(t, 'Color',colorList(site, :));
            xlabel('Model Estimate');
        end
        title(titleStr{r});
    end
    % look at correlations with total charge
    subplot(1, 3, 3)
    eid =  find(T.electrode==site & ~isnan(T.brightness));
    if ~isempty(eid)
    Texp = T(eid, :);
    x = [Texp.totalcharge];  y = [Texp.brightness] ;
    x = x + sd*(rand(size(x))-0.5); y = y + sd*(rand(size(x))-0.5);
    h = scatter(x, y, sz, 'MarkerFaceColor', colorList(site, :), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha);
    hold on
    t = text(15, site/2, ['corr = ', num2str(round(corr(x, y), 3))]);
    set(t, 'Color',colorList(site, :));
    xlabel('Total Charge');
    end
end
for i = 1:3
    subplot(1, 3, i)
    ylabel('Reported Brightness');
%     set(gca, 'XLim', [0 25]);
%     set(gca, 'YLim', [0 11]);
end

