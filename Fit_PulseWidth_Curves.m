%% Fit_PulseWidth_Curves.m
%
% Fits a variety of pulse width data, examining how threshold varies with
% pulse width
%
% Written by GMB & IF 
% 25/02/2023 moved into clean folder (IF)

tp = p2p_c.define_temporalparameters();

% Calculate model response to a standard trial
clear standard_trl;
standard_trl.pw = 0.001;
standard_trl.amp = 3;
standard_trl.dur = 1 ;
standard_trl.freq = 50;
standard_trl.simdur = 3; %sec
standard_trl = p2p_c.define_trial(tp,standard_trl);
standard_trl= p2p_c.convolve_model(tp, standard_trl);
tp.thresh_resp = standard_trl.maxresp;% use that as the threshold response

% Now on to the real experiments...  Note, Henderson79 has no pd
% experiment, and Girvin79 has two.
% Data from Brindley et al. 1968 are excluded because measured thresholds
% were  unexpectedly non-monotonic as a function of frequency, suggesting
% an electronics issue
paperList = {'Dobelle74'; 'Dobelle79'; 'Henderson79' ;'Girvin79'; 'Fernandez21';};
exList = {{'pd'},{'pd'},{''},{'pd1','pd50'},{'pd'}};
legStr = {'Dobelle 1974, 50Hz 1 sec duration',...
    'Dobelle 1979, 50Hz 0.5 sec duration',...
    'Girvin 1979, 1 pulse 0.5 sec duration','Girvin 1979 50Hz, 0.5 sec duration',...
    'Fernandez 2021, 300Hz 0.167 sec duration'};

clear x y1 y2
count = 0;
for i = 1:length(paperList)
    for j=1:length(exList{i})
        disp(sprintf('Paper: %s, experiment: %s \n',paperList{i},exList{i}{j}))
        T = readtable(['datasets/', paperList{i}, '_data.xlsx']);
        eid =strcmp(T.experiment,exList{i}{j});
        if sum(eid)
            count = count+1;
            Texp = T(eid,:);
            [err, thresh] = p2p_c.loop_find_threshold(tp,Texp);

            Tsim = Texp;
            Tsim.freq = standard_trl.freq*ones(size(Tsim.freq));
            Tsim.dur = standard_trl.dur*ones(size(Tsim.dur));

            % get model response to standard stimulus at these pd's
            [err, standard_thresh] = p2p_c.loop_find_threshold(tp,Tsim);

            % regression to scale the actual data.
            scFac1 =  standard_thresh'/Texp.amp';

            % save results to plot later.
            x{count} = log(Texp.pw*1000);
            y1{count} = Texp.amp*scFac1;

            scFac2 = standard_thresh'/thresh';
            y2{count} = thresh*scFac2;
        else
            disp('No data points found.');
        end
    end
end

%% predict smooth pd curve using standard trial parameters. 
% Should go through our standard trial's point: amp = 3 at pw = 1msec

% make fake data table
pdList = exp(linspace(log(0.05), log(6.4), 40))/1000;
freq = standard_trl.freq*ones(size(pdList));
dur = standard_trl.dur*ones(size(pdList));

Tstandard = table(pdList', freq', dur');
Tstandard.Properties.VariableNames = {'pw','freq','dur'};
% get thresholds
[err, standard_thresh] = p2p_c.loop_find_threshold(tp,Tstandard);

%% Plotting real data
% (fast after running previous chunks)

sz = 150; % symbol size
xtick = [.05,.1,.2,.4,.8,1.6,3.2,6.4];
colList = {'b', 'r', 'c', 'm', 'g'};
alpha = .5;  % transparency
sd = .15; % jitter

fontSize = 12;
figure(1); clf; hold on
plot(log(pdList*1000),standard_thresh,'k-','LineWidth',2);
for i=1:length(x)
    h(i)= scatter(x{i}+2*sd*(rand(size(x{i}))-.5),y1{i},sz , 'MarkerFaceColor', colList{i},...
        'MarkerEdgeColor','none','MarkerFaceAlpha',alpha);
end
set(gca,'XTick',log(xtick)); logx2raw; ylabel('Threshold'); widen; grid;
xlabel('Pulse Width (msec)'); ylabel('Amplitude @ Threshold');
legend(h,legStr,'Location','NorthEast');
set(gca,'YLim',[0,20]);
set(gca,'FontSize',fontSize);

% Figure 2 holds the predicted thresholds for each experiment scaled
% to match the standard stimulus curve.  The curves match perfectly - the
% shape
sd = 0;
figure(2); clf; hold on
plot(log(pdList*1000),standard_thresh,'k-','LineWidth',2);
for i=1:length(x)
    h(i)= scatter(x{i}+2*sd*(rand(size(x{i}))-.5),y2{i},sz , 'MarkerFaceColor', colList{i},...
        'MarkerEdgeColor','none','MarkerFaceAlpha',alpha);
end
set(gca,'XTick',log(xtick)); logx2raw; ylabel('Threshold'); widen; grid;
xlabel('Pulse Width (msec)'); ylabel('Amplitude @ Threshold');
legend(h,legStr,'Location','NorthEast');
set(gca,'YLim',[0,20]);
set(gca,'FontSize',fontSize);
title('Scaled model predictions');



