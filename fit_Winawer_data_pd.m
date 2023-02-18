%Fit_Winawer_Data_pd.m
% fits Girvin data using a probability summation model
clear all; close all

sites = 1:5;

T = readtable('datasets/Winawer2016_data.xlsx');
tp.tSamp = 1000;
tp.dt  = 1e-6;
tp.nreps = 12; % the duration of the binary search to find threshold
tp.tau1 = [ 0.000475]; % fixed based on Nowak and Bullier, 1998
tp.tau1Fac = [ 1]; %.relative sensitivity of the two tau1 if using more than one
tp.refractoryPeriod = [85 10];
tp.ncascades  = 3;
tp.win = .5/tp.tSamp;
tp.probsumflag = 0;
tp.gammaflag = 1;
tp.refractoryflag = 0;
tp.tau2 =0.0868; % used if second stage is gamma
colorList  = [1 0 0; 0 1 0; .6 .4 0; .3 .3  1; 0 0 1 ]; % roughly match Winawer paper
tp.model = 'compression'; % weirdly, this non-linearity used in different places in the model for gamma and probsum
tp.scFac = 1/26; tp.sigma = .15; tp.mean =0.6; tp.power =26; %only some used, depends on tp.model
tp.thresh_resp = 0.75; % 2 alt FC detection threshold for Girvin, fixed for probsum second stage
FITFLAG = 1;

for site = 3

 eid = T.electrode ==site;
 Texp = T(eid,:);    
 Texp.brightness( Texp.brightness==-1000) = NaN; % weird hack because otherwise brightness stuck as a cell array

    if FITFLAG
        freeParams = {'power'};
        tp = fit('p2p_c.fit_brightness',tp,freeParams,Texp);
    end

    trl= p2p_c.loop_convolve_model(tp, Texp);

    subplot(2,1,1)
    x = [trl.maxresp]; y = [Texp.brightness] ; 
    x= reshape(x, length(x), 1);  y = reshape(y, length(x), 1);
    plot(x, y, 'o', 'MarkerFaceColor', colorList(site, :), 'MarkerEdgeColor', 'none'); hold on
    ind = ~isnan(x) & ~isnan(y);
     disp(corr(x(ind), y(ind)));   
     xlabel('model estimate');
    ylabel('reported Brightness')

end
