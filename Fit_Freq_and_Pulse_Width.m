%% Fit_Freq_and_Pulse_Width.m
%
% Fits a variety of frequency and pulse width data, examining how threshold varies with
% pulse width
%
% Written by GMB & IF
% 5/08/2023 combined frequency and pulse width fitting into this script,
% which uses the same scale factor for electrodes shared across
% experiments.
% 
% 25/02/2023 moved into clean folder (IF)
tp = p2p_c.define_temporalparameters();

% Calculate model response to a standard trial
clear standard_trl;
standard_trl.pw =  0.00025;  %0.00025;  why did we use .00025 for frequency data and .001 for pw?  
standard_trl.amp = 3;
standard_trl.dur = .5 ;  % why did we use 0.5 for frequency and 1 for pw?
standard_trl.freq = 50;
standard_trl.simdur = 6; %sec
standard_trl = p2p_c.define_trial(tp,standard_trl);
standard_trl= p2p_c.convolve_model(tp, standard_trl);
tp.thresh_resp = standard_trl.maxresp;% use this as the threshold response

%% Define the papers and experiments
paper(1).name = 'Dobelle74';
paper(1).electrode = {'A1','A2','A3','A4','B2','B3','B4','C1','C2','C3','C4'};
paper(1).exp = {'pd','freq'};
paper(1).legend = {'Dobelle 1974, 50Hz 1 sec duration',...
                   'Dobelle 1974, 0.5 msec pulse width, 1 sec duration'};
paper(1).x = {[],[]}; % will hold x-axis values (frequency or pw)
paper(1).y = {[],[]}; % will hold y-axis values (scaled thresholds)
paper(1).color = [0 0 1; 0 .3 .7]

paper(2).name = 'Girvin79';
paper(2).electrode = {'A'};
paper(2).exp = {'pd1','pd50','freq'};
paper(2).legend = {'Girvin 1979, 1 pulse 0.5 sec duration',...
                   'Girvin 1979, 50Hz, 0.5 sec duration',...
                   'Girvin 1979, 0.25 msec pulse width, 0.5 sec duration'};
paper(2).x = {[],[],[]};
paper(2).y = {[],[],[]};
paper(2).color = [ 0 1 0; 0 .7 0; 0 .5 0];

paper(3).name = 'Dobelle79';
paper(3).electrode = {'A'};
paper(3).exp = {'pd'};
paper(3).legend = {'Dobelle 1979, 50Hz 0.5 sec duration'};
paper(3).x = {[]};
paper(3).y = {[]};
paper(3).color = [0 .5 .5];

paper(4).name = 'Fernandez21';
paper(4).electrode = {'A'};
paper(4).exp = {'pd', 'freq'};
paper(4).legend = {'Fernandez 2021, 300Hz 0.167 sec duration', ...
                                  'Fernandez 2021, 0.17msec pw 0.167 sec duration' };
paper(4).x = {[], []};
paper(4).y = {[], []};
paper(4).color = [1 0 0 ; .5 0 0];

%%
% Compile the thresholds across papers and experiments, scaling by a common
% scale factor for each electrode.

for paperNum = 1:length(paper)
    disp(sprintf('Paper: %s',paper(paperNum).name));
    T = readtable(['datasets/', paper(paperNum).name, '_data.xlsx']);
    
    % loop through electrodes before experiments
    for electrodeNum=1:length(paper(paperNum).electrode)
        standard_thresh = [];
        amp = [];
        disp(sprintf('  Electrode: %s ',paper(paperNum).electrode{electrodeNum}));
        for expNum=1:length(paper(paperNum).exp)
            fprintf('    Experiment: %5s',paper(paperNum).exp{expNum});
            eid = strcmp(T.experiment,paper(paperNum).exp{expNum}) & ...
                strcmp(T.electrode,paper(paperNum).electrode{electrodeNum});
            disp(sprintf('        %d points',sum(eid)));
                  
            Texp = T(eid,:);
            [err, thresh] = p2p_c.loop_find_threshold(tp,Texp);
            
            % concatenate the amplitudes at threshold for this electrode
            amp = [amp;Texp.amp];
            Tsim = Texp;
            if strcmp(paper(paperNum).exp{expNum}(1:2),'pd')
                Tsim.freq = standard_trl.freq*ones(size(Tsim.freq));                
            else
                Tsim.pw = standard_trl.pw*ones(size(Tsim.freq));
            end
            Tsim.dur = standard_trl.dur*ones(size(Tsim.dur));
            
            % get model response to standard stimulus at these pd's
            [err, tmp_thresh] = p2p_c.loop_find_threshold(tp,Tsim);
            
            % concatenate the standard responses for this electrode
            standard_thresh = [standard_thresh;tmp_thresh];
        end % expNum
                
        % regression to scale the actual data.
        scFac =  standard_thresh'/amp';
        disp(sprintf('      scFac = %5.4f',scFac));
        
        % Apply this scale factor and save the results to plot later.
        
        for expNum=1:length(paper(paperNum).exp)
            eid = strcmp(T.experiment,paper(paperNum).exp{expNum}) & ...
                strcmp(T.electrode,paper(paperNum).electrode{electrodeNum});
            Texp = T(eid,:);
            
            if strcmp(paper(paperNum).exp{expNum}(1:2),'pd') % pd experiment
                paper(paperNum).x{expNum} = [paper(paperNum).x{expNum};log(Texp.pw*1000)];
                paper(paperNum).y{expNum} = [paper(paperNum).y{expNum};Texp.amp*scFac];    % note the scFac
            else  % frequency experiment
                paper(paperNum).x{expNum} = [paper(paperNum).x{expNum};log(Texp.freq)];
                paper(paperNum).y{expNum} = [paper(paperNum).y{expNum};Texp.amp*scFac];    % note the scFac
            end
        end
    end
end

%% predict smooth pw and freq curves using standard trial parameters.
%  Should go through our standard trial's point: amp = 3 at pw = 1msec and freq = 50Hz

clear standard_thresh

% make fake data table (pw)
pdList = exp(linspace(log(0.05), log(6.4), 40))/1000;
%pdList =exp([   -2.7806   -2.3026   -2.0794   -1.7720    -1.3863    -0.9163    -0.6931    -0.2231       0     0.6931   1.6094]);
%pdList = unique(pdList); pdList = sort(pdList);
freq = standard_trl.freq*ones(size(pdList));
dur = standard_trl.dur*ones(size(pdList));

Tstandard = table(pdList', freq', dur');
Tstandard.Properties.VariableNames = {'pw','freq','dur'};
% get thresholds
[err, standard_thresh.pw] = p2p_c.loop_find_threshold(tp,Tstandard);

% make fake data table (frequency)
freqList = exp(linspace(log(8), log(2000), 40));
%freqList = exp([2.4849    2.5257    3.2189   3.9120    4.6052    5.2983    5.9915  6.6846    7.3778]);
%freqList = unique(freqList); freqList = sort(freqList);
pw = standard_trl.pw*ones(size(freqList));
dur = standard_trl.dur*ones(size(freqList));

Tstandard = table(pw', freqList', dur');
Tstandard.Properties.VariableNames = {'pw','freq','dur'};
% get thresholds
[err, standard_thresh.freq] = p2p_c.loop_find_threshold(tp,Tstandard);

%% plotting
% (fast after running previous chunks)

% figure 1: pw
% figure 1: freq

sz = 150; % symbol size
xtick.pd =[.05,.1,.2,.4,.8,1.6,3.2,6.4];

xtick.freq= [8,16,32,64,128,256,512,1024];

alpha = .5;  % transparency
sd = 0.1; % jitter
fontSize = 12;
clear h

% plot the standard thresholds

figure(1); clf; hold on
plot(log(pdList*1000),standard_thresh.pw,'k-','LineWidth',2);

figure(2); clf; hold on
plot(log(freqList),standard_thresh.freq,'k-','LineWidth',2);

h.pd = [];
h.freq = [];
legStr.pd = {};
legStr.freq = {};

allpd_x= [];
allpd_y = [];
allfreq_x = [];
allfreq_y = [];
for paperNum = 1:length(paper)
    for expNum = 1:length(paper(paperNum).exp)
        if strcmp(paper(paperNum).exp{expNum}(1:2),'pd')
            figure(1)
            h.pd(end+1)= scatter(paper(paperNum).x{expNum} + ...
                2*sd*(rand(size(paper(paperNum).x{expNum}))-0.5),paper(paperNum).y{expNum} + ...
                sd*(rand(size(paper(paperNum).x{expNum}))-0.5) ,sz , 'MarkerFaceColor', paper(paperNum).color(expNum, :),...
                'MarkerEdgeColor','none','MarkerFaceAlpha',alpha);
            allpd_x = cat(1, allpd_x, paper(paperNum).x{expNum});
            allpd_y =cat(1, allpd_y, paper(paperNum).y{expNum});
            legStr.pd{end+1} = paper(paperNum).legend{expNum};            
        else % frequency experiment
            figure(2)
            h.freq(end+1)= scatter(paper(paperNum).x{expNum} + ...
                2*sd*(rand(size(paper(paperNum).x{expNum}))-0.5),paper(paperNum).y{expNum} + ...
                sd*(rand(size(paper(paperNum).x{expNum}))-0.5) ,sz , 'MarkerFaceColor', paper(paperNum).color(expNum, :),...
                'MarkerEdgeColor','none','MarkerFaceAlpha',alpha);
            legStr.freq{end+1} = paper(paperNum).legend{expNum};     
            allfreq_x = cat(1, allfreq_x, paper(paperNum).x{expNum});
            allfreq_y =cat(1, allfreq_y, paper(paperNum).y{expNum});
        end
    end    
end

% axes, legends and stuff...

figure(1)
set(gca,'XTick',log(xtick.pd)); logx2raw; ylabel('Threshold'); widen; grid;
xlabel('Pulse Width (msec)'); ylabel('Amplitude @ Threshold');
legend(h.pd,legStr.pd,'Location','NorthEast');
set(gca,'YLim',[0,12]);
set(gca,'FontSize',fontSize)

allpd_predy = interp1(log(pdList*1000), standard_thresh.pw, allpd_x);
[r, p] = corr(allpd_y, allpd_predy);
disp(['r =  ', num2str(round(r, 3)), ' p = ', num2str(round(p, 4)), ' df = ', num2str(length(allpd_y)-2)]);

figure(2)
set(gca,'XTick',log(xtick.freq)); logx2raw(exp(1),0); ylabel('Threshold'); widen; grid;
xlabel('Frequency (Hz)'); ylabel('Amplitude @ Threshold');
legend(h.freq,legStr.freq,'Location','NorthEast');
set(gca,'XLim',log([9,2000]));
set(gca,'YLim',[0,6]);
set(gca,'FontSize',fontSize);
title('Scaled model predictions');
allfreq_predy = interp1(log(freqList), standard_thresh.freq, allfreq_x);
[r, p] = corr(allfreq_y, allfreq_predy);
disp(['r =  ', num2str(round(r, 3)), ' p = ', num2str(round(p, 4)), ' df = ', num2str(length(allfreq_predy)-2)]);




