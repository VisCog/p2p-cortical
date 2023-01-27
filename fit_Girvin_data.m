%Fit_Girvin_Data.m

IPI = 1;  % plot inter-pulse interval instead of frequency

T = readtable('datasets/Girvin79_data.xlsx');

tp.tSamp = 1000;
tp.experimentList = {'pd50','pd1','freq','dur','dur50'}; %list of sets to fit.  Each set gets it's own thresh_resp


% tau1 = 0.00035, tau2 = 0.015, power =  2.63 err= 8.4436
%       pd50: 4344.26
%        pd1: 1027.93
%       freq: 2581.42
%        dur: 3753.18
%      dur50: 1429.99


% Fitting pd1 and freq, setting tau1 to .0035 and tau2 t0 .01, thresh_resp
% .013.  and setting power and scFac free. Nice fit, but power is 1.963 >
% 
%                model: 'power'
%                 tau1: 3.5000e-04
%                 tau2: 0.0100
%                tSamp: 1000
%            asymptote: 1
%       experimentList: {'pd1'  'freq'}
%          thresh_resp: [0.0130 0.0130]
%                power: 1.9630
%            ncascades: 3
%                   dt: 1.0000e-06
%             spikeDur: 2.5000e-04
%                scFac: [713.1069 1.6308e+03]
%     refractoryPeriod: 0.0046

% freeParams = {'scFac'} gives:
%                model: 'power'
%                 tau1: 3.5000e-04
%                 tau2: 0.0100
%                tSamp: 1000
%            asymptote: 1
%       experimentList: {'pd50'  'pd1'  'freq'  'dur'  'dur50'}
%          thresh_resp: [0.0130 0.0130 0.0130 0.0130 0.0130]
%                power: 1
%            ncascades: 3
%                   dt: 1.0000e-06
%             spikeDur: 2.5000e-04
%                scFac: [0.8823 0.6689 0.8679 1.0151 0.6063]
%     refractoryPeriod: 0.0046
% err= 7.5169

%                model: 'power'
%                 tau1: 3.5000e-04
%                 tau2: 0.0152
%                tSamp: 1000
%            asymptote: 1
%       experimentList: {1×5 cell}
%          thresh_resp: [0.0130 0.0130 0.0130 0.0130 0.0130]
%                power: 1
%            ncascades: 3
%                   dt: 1.0000e-06
%             spikeDur: 2.5000e-04
%                scFac: [0.9139 1.0189 1.1427 1.3593 0.7584]
%     refractoryPeriod: 0.0046
%  err= 6.0712

% freeParams = {'scFac','tau2'} gives:  Not great fit for 4 pulses/trial
% tp.experimentList = {'pd50','pd1','freq','dur','dur50'}; %list of sets to fit.  Each set gets it's own thresh_resp
tp.scFac =  [0.8823 0.6689 0.8679 1.0151 0.6063];
tp.power = 	1;
tp.tau1 = .00035;
tp.tau2 = 0.0032;
tp.tau2 = .01;
tp.dt  = 1e-6;
tp.model = 'power';
tp.ncascades  = 3;
tp.refractoryPeriod = 1/217;
tp.thresh_resp = .0130;
freeParams = {'scFac','tau2'};
% tpBest = fit('p2p_c.getErrCronaxie',tp,freeParams,T);%
 %tp = tpBest;

tic
[err,thresh] = p2p_c.getErrCronaxie(tp,T);
toc
%% Plot the fits

%% figure 4a - 'pd50'
eid =strcmp(T.experiment,'pd50');

x = T.pw(eid)*1000;
y = T.amp(eid);
s = T.sem(eid);

figure(1)
clf
hold on
errorbar(log(x),y,s,'LineStyle','none','Color','k');
plot(log(x),y,'ko'); p2p_c.fillSymbols
plot(log(x),thresh(eid),'k-');
set(gca,'XTick',log(x));
logx2raw
xlabel('Pulse Duration (msec)');
ylabel('Threshold');
widen
grid
title(sprintf('50 Hz (%g msec IPI) pulses',1000/50));

%% figure 4b - 'pd1'
eid =strcmp(T.experiment,'pd1');

x = T.pw(eid)*1000;
y = T.amp(eid);
s = T.sem(eid);

figure(2)
clf
hold on
errorbar(log(x),y,s,'LineStyle','none','Color','k');
plot(log(x),y,'ko'); p2p_c.fillSymbols
plot(log(x),thresh(eid),'k-');
set(gca,'XTick',log(x));
logx2raw
xlabel('Pulse Duration (msec)');
ylabel('Threshold');
widen
grid
title('Single pulses');


%% figure 5 - 'freq'
eid =strcmp(T.experiment,'freq');

x = T.freq(eid);
y = T.amp(eid);
s = T.sem(eid);

figure(3)
clf
hold on
errorbar(log(x),y,s,'LineStyle','none','Color','k');
plot(log(x),y,'ko'); p2p_c.fillSymbols
plot(log(x),thresh(eid),'k-');
set(gca,'XTick',log(x));
logx2raw(exp(1),0);

if IPI
    str = num2str(1000./x);
    set(gca,'XTickLabels',str)
    xlabel('Inter-pulse interval (msec)');
else
xlabel('Frequency (Hz)');
end
ylabel('Threshold');
widen
grid
title('0.5 sec duration, .25 msec pw');
%% figure 6a - 'dur'
eid =strcmp(T.experiment,'dur');

x = T.dur(eid)*1000;
y = T.amp(eid);
s = T.sem(eid);

figure(4)
clf
hold on
errorbar(log(x),y,s,'LineStyle','none','Color','k');
plot(log(x),y,'ko'); p2p_c.fillSymbols
plot(log(x),thresh(eid),'k-');
set(gca,'XTick',log(x-.5));
logx2raw
if IPI
    str = num2str((x-.5)/4);
    set(gca,'XtickLabels',str);
    xlabel('Inter-pulse interval (msec)');
else
    xlabel('Stimulus Duration (msec)');
end
ylabel('Threshold');
widen
grid
title('4 pulses/trial');

%% figure 6b - 'dur50'
eid =strcmp(T.experiment,'dur50');

x = T.dur(eid)*1000;
y = T.amp(eid);
s = T.sem(eid);

figure(5)
clf
hold on
errorbar(log(x),y,s,'LineStyle','none','Color','k');
plot(log(x),y,'ko'); p2p_c.fillSymbols
plot(log(x),thresh(eid),'k-');
set(gca,'XTick',[log(x(1));log(x(2:end)-.5)]);
logx2raw(exp(1),0)
xlabel('Stimulus Duration (msec)');
ylabel('Threshold');
widen
grid
title(sprintf('50 Hz (%g msec IPI) pulses',1000/50));
