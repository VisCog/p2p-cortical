%Fit_Girvin_Data.m

IPI = 1;  % plot inter-pulse interval instead of frequency

T = readtable('datasets/Girvin79_data.xlsx');

tp.tSamp = 1000;
tp.dt  = 1e-6;
% tp.experimentList = {'pd1','pd50','freq','dur','dur50'}; %list of sets to fit.  Each set gets it's own thresh_resp
% pd1, pulse duration varies, 1 pulse
% pd50, pulse duration varies, 50Hz
% freq, 05s pulse train, frequency varies
% dur, 4 pulses, frequency (duration to do all 4 pulses) varies
% dur 50, 50 Hz stimulation, pulse train duration varies (i.e. number of
% pulses)

tp.tau1 = [0.0003]; % fixed based on Nowak and Bullier, 1998
tp.refractoryPeriod = [70 100];
tp.ncascades  = 3; tp.tau2 = 20/1000;
tp.thresh_resp = 14; tp.power = 1;
tp.model = 'compression'; tp.sigma = .5; tp.scFac = 1;
%figure(1); clf; y = p2p_c.nonlinearity(tp,0:.1:2); plot(0:.1:2, 0:.1:2, 'k'); hold on; plot(0:.1:2,y, 'r'); 

%% single pulse
eid =strcmp(T.experiment,'pd1'); Tpd1 = T(eid,:);  
 trl_pd1 = p2p_c.loop_convolve_model(tp, Tpd1);

%%  'freq' (figure 5)
eid =strcmp(T.experiment,'freq'); Tfreq = T(eid,:);
trl_freq = p2p_c.loop_convolve_model(tp, Tfreq);

%%  'pd1' (figure 4b)
tp.thresh_resp = 14;
%tp = fit('p2p_c.loop_find_threshold',tp,{'scFac'},Tpd1);
[err,thresh] = p2p_c.loop_find_threshold(tp,Tpd1);
x = Tpd1.pw*1000; y = Tpd1.amp; s = Tpd1.sem;
figure(1); clf; hold on
errorbar(log(x),y,s,'LineStyle','none','Color','k');
plot(log(x),y,'ko'); p2p_c.fillSymbols;
plot(log(x),thresh,'k-');
set(gca,'XTick',log(x)); logx2raw; xlabel('Pulse Duration (msec)'); ylabel('Threshold'); widen; grid; title('Single pulses');
tp.scFac = 1;

%% freq
tp.thresh_resp = 14;
%tp = fit('p2p_c.loop_find_threshold',tp,{'scFac'},Tfreq);
[err,thresh] = p2p_c.loop_find_threshold(tp,Tfreq);
x = Tfreq.freq; y = Tfreq.amp; s = Tfreq.sem;
figure(2); clf; hold on
errorbar(log(x),y,s,'LineStyle','none','Color','k');
plot(log(x),y,'ko'); p2p_c.fillSymbols
plot(log(x),thresh,'k-');
set(gca,'XTick',log(x)); logx2raw(exp(1),0);
if IPI
    str = num2str(1000./x);  set(gca,'XTickLabels',str);     xlabel('Inter-pulse interval (msec)');
else
    xlabel('Frequency (Hz)');
end
ylabel('Threshold'); widen; grid; title('0.5 sec duration, .25 msec pw');
return

%% figure 4a - 'pd50'
tp.scFac = [50];
eid =strcmp(T.experiment,'pd50'); Tpd50 = T(eid,:);
[err,thresh] = p2p_c.loop_trials(tp,Tpd50);
x = Tpd50.pw*1000; y = Tpd50.amp; s = Tpd50.sem;
figure(3); clf; hold on
errorbar(log(x),y,s,'LineStyle','none','Color','k');
plot(log(x),y,'ko'); p2p_c.fillSymbols
%plot(log(x),thresh,'k-');
set(gca,'XTick',log(x)); logx2raw; xlabel('Pulse Duration (msec)'); ylabel('Threshold'); widen; grid
title(sprintf('50 Hz (%g msec IPI) pulses',1000/50));


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
%plot(log(x),thresh(eid),'k-');
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
%plot(log(x),thresh(eid),'k-');
set(gca,'XTick',[log(x(1));log(x(2:end)-.5)]);
logx2raw(exp(1),0)
xlabel('Stimulus Duration (msec)');
ylabel('Threshold');
widen
grid
title(sprintf('50 Hz (%g msec IPI) pulses',1000/50));

