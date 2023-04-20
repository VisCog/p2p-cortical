% Winawer_Brightness_Neuron2016
%
% Simulates phosphene brightness as a function of electrode location and pulse
% train
%
% Winawer J, Parvizi J. Linking Electrical Stimulation of Human Primary Visual Cortex,
% %Size of Affected Cortical Area, Neuronal Responses, and Subjective Experience.
% Neuron. 2016 Dec 21;92(6):1213-1219. doi: 10.1016/j.neuron.2016.11.008. Epub
%  2016 Dec 8. PMID: 27939584; PMCID: PMC5182175.
%
% Written by GMB & IF
% 25/02/2023 moved into clean folder (IF)

clear
figure(1); clf

T = readtable('datasets/Winawer2016_data.xlsx');
T.brightness( T.brightness==-1000) = NaN; % weird hack because otherwise brightness stuck as a cell array
eid = find(T.amp~=0 & ~isnan(T.brightness));
T = T(eid, :);
tp = p2p_c.define_temporalparameters();
tp.model = 'compression';
tp.slope =  0.4300;
tp.power = 10;
rng; %(11171964) % fix the random number generator, used Geoff's birthday in paper

FITFLAG = 0;
if FITFLAG
    Texp = T;
    freeParams = { 'slope'};
    tp = fit('p2p_c.fit_brightness',tp,freeParams,Texp);
end
colorList  = [1 0 0; 0 1 0; 1 .7 0; .3 .3  1; 0 0 1 ]; % roughly match colors to Winawer paper
sd = .2;
for site =1:5
    eid =  find(T.electrode==site);
    if ~isempty(eid)
        Texp = T(eid, :);
        trl= p2p_c.loop_convolve_model(tp, Texp);
        x = [trl.maxresp]; y = [Texp.brightness] ;
        x= reshape(x, length(x), 1);  y = reshape(y, length(x), 1);
        ind = ~isnan(x) & ~isnan(y);
        if sum(ind)>0
            h = scatter(x+sd*randn(size(x)), y+sd*randn(size(x)), 'o', 'MarkerFaceColor', colorList(site, :),  'MarkerEdgeColor','none','MarkerFaceAlpha',.5); hold on
            corrval = corr(x(ind), y(ind));
            xlabel('Model Estimate');
            ylabel('Reported Brightness');
            set(gca, 'XLim', [0 11]);
            set(gca, 'YLim', [0 11]);
     %       plot([0 11], [0 11], 'k')
            t = text(8, site/2, ['corr = ', num2str(round(corrval, 3))]);
            set(t, 'Color',colorList(site, :) )
        end
    end
end
axis square

