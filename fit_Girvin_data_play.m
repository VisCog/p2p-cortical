%Fit_Girvin_Data.m
clear all;
FITFLAG = 0;
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

tp.tau1 = [ 0.0003 ]; % fixed based on Nowak and Bullier, 1998
tp.tau1Fac = [ 1]; %.5 0.19
tp.nreps = 12; % the duration of the binary search to find threshold

tp.refractoryPeriod = [85 10];
tp.ncascades  = 3; tp.tau2 =12.5/1000;
tp.model  = 'linear';
tp.scFac = 1;tp.sigma = .1; tp.mean =0.5; tp.power =.8;
tp.model = 'compression'; tp.secondstate= 'probsum';
input= linspace(0, tp.power*2, 100);
figure(11); clf; y = p2p_c.nonlinearity(tp,input);
plot(input, input, 'k'); hold on;
plot(input,y, 'r');

%% single pulse

eid =strcmp(T.experiment,'pd1'); Tpd1 = T(eid,:);
trl_pd1 = p2p_c.loop_convolve_model(tp, Tpd1);
tp.thresh_resp = mean([trl_pd1.maxresp]);
disp(['intial single pulse thresh = ', num2str( tp.thresh_resp )]);
if FITFLAG
    freeParams = {'thresh_resp'};
    tp = fit('p2p_c.loop_find_threshold',tp,freeParams,Tpd1);
end
figure(2); clf; subplot(2, 1, 1)
x = Tpd1.pw*1000; y = Tpd1.amp; s = Tpd1.sem;
errorbar(log(x),y,s,'LineStyle','none','Color','r');
plot(log(x),y,'ro'); p2p_c.fillSymbols; hold on
set(gca,'XTick',log(x)); logx2raw; xlabel('Pulse Duration (msec)'); ylabel('Threshold'); widen; grid; title('Single pulses');
ylabel('Amplitude Threshold')

[err, thresh] = p2p_c.loop_find_threshold(tp,Tpd1);
plot(log(x),thresh,'k-');
subplot(2,1,2)
plot(log(x), [trl_pd1.maxresp], 'ko', 'MarkerFaceColor', 'k'); hold on
plot(log(x),tp.thresh_resp.*ones(1,length([trl_pd1])), 'r');
set(gca,'XTick',log(x)); logx2raw; xlabel('Pulse Duration (msec)'); ylabel('thresh'); widen; grid;
ylabel('Model Response @ Amplitude Threshold');
mse = sqrt(mean(([trl_pd1.maxresp]-tp.thresh_resp).^2))./tp.thresh_resp;
disp(['single pulse mse = ', num2str(mse)]);

figure(10); clf
clr =[0 0 1 ; 0 1 0];
pw = [0.0001:0.0001:001];
for i = 1:length(tp.tau1)
    f.tau = tp.tau1(i); f.amp = tp.tau1Fac(i)/1000;
    out(i, :)  = p2p_c.chronaxie(f, pw);
end
for i = 1:length(tp.tau1)
    plot(log(pw*1000), out(i, :), '--', 'Color', clr(i, :)); hold on
end
plot(log(pw*1000), min(out, [], 1), '--', 'Color', [.5 .5 .5], 'LineWidth', 2); hold on
return
% if +ve then the fitted amp threshold would be higher than the data
% if -ve the fitted amplitude threshold would be lower than the data

%%  'freq' (figure 5)
eid =strcmp(T.experiment,'freq'); Tfreq = T(eid,:);
trl_freq = p2p_c.loop_convolve_model(tp, Tfreq);
tp.thresh_resp = mean([trl_freq.maxresp]);
disp(['freq thresh = ', num2str( tp.thresh_resp )])
figure(3); clf; subplot(2, 1,1)
x = Tfreq.freq; y = Tfreq.amp; s = Tfreq.sem;
errorbar(log(x),y,s,'LineStyle','none','Color','r');
plot(log(x),y,'ro'); p2p_c.fillSymbols; hold on
set(gca,'XTick',log(x)); logx2raw(exp(1),0);     xlabel('Frequency (Hz)'); widen; grid; title('Frequency, 0.5 s dur, .25 msec pw');
ylabel('Amplitude Threshold');
[err, thresh] = p2p_c.loop_find_threshold(tp, Tfreq);
plot(log(x),thresh,'k-');

subplot(2,1,2)
plot(log(x),[trl_freq.maxresp], 'ko', 'MarkerFaceColor', 'k'); hold on
plot(log(x),tp.thresh_resp.*ones(1,length(trl_freq)), 'r');
set(gca,'XTick',log(x)); logx2raw(exp(1),0);xlabel('Frequency (Hz)');widen; grid;
ylabel('Model Response @ Amplitude Threshold')
mse = sqrt(mean(([trl_freq.maxresp]-tp.thresh_resp).^2))./tp.thresh_resp;
disp(['freq mse = ', num2str(mse)])

%%  'pd50' (figure 4a)
tp.tau1 = [ 0.00022]; % fixed based on Nowak and Bullier, 1998
tp.tau1Fac = [ 1 ]; %3
eid =strcmp(T.experiment,'pd50'); Tpd50 = T(eid,:);
trl_pd50= p2p_c.loop_convolve_model(tp, Tpd50);
tp.thresh_resp = mean([trl_pd50.maxresp]);
disp(['pd50 thresh = ', num2str( tp.thresh_resp )])
figure(4); clf; subplot(2,1,1)
x = Tpd50.pw*1000; y = Tpd50.amp; s = Tpd50.sem;
errorbar(log(x),y,s,'LineStyle','none','Color','r');
plot(log(x),y,'ro'); p2p_c.fillSymbols; hold on
set(gca,'XTick',log(x)); logx2raw; xlabel('Pulse Duration (msec)'); ylabel('Threshold'); widen; grid
title(sprintf('50 Hz (%g msec IPI) pulses',1000/50));
ylabel('Amplitude Threshold');
[err, thresh] = p2p_c.loop_find_threshold(tp, Tpd50);
plot(log(x),thresh,'k-', 'LineWidth', 2); hold on

subplot(2,1,2)
plot(log(x), [trl_pd50.maxresp], 'ko', 'MarkerFaceColor', 'k'); hold on
plot(log(x),tp.thresh_resp.*ones(1,length(trl_pd50)), 'r');
set(gca,'XTick',log(x)); logx2raw; xlabel('Pulse Duration (msec)'); ylabel('Threshold'); widen; grid;
mse = sqrt(mean(([trl_pd50.maxresp]-tp.thresh_resp).^2))./tp.thresh_resp;
disp(['pd50 mse = ', num2str(mse)])

%%  'dur/4 pulse' (figure 6a)
eid =strcmp(T.experiment,'dur'); T4pulse = T(eid,:);
trl_4pulse= p2p_c.loop_convolve_model(tp, T4pulse);
tp.thresh_resp = mean([trl_4pulse.maxresp]);
disp(['4 pulse thresh = ', num2str( tp.thresh_resp )])
figure(5); clf; subplot(2,1,1)
x = T4pulse.freq; y = T4pulse.amp; s = T4pulse.sem;
errorbar(log(x),y,s,'LineStyle','none','Color','r');
plot(log(x),y,'ro'); p2p_c.fillSymbols; hold on
set(gca,'XTick',log(x)); logx2raw; xlabel('Freq for 4 pulses (Hz)'); widen; grid ;title('4 pulses/trial');
title('4 pulses');
ylabel('Amplitude Threshold');
[err, thresh] = p2p_c.loop_find_threshold(tp, T4pulse);
plot(log(x),thresh,'k-');

subplot(2,1,2)
plot(log(x), [trl_4pulse.maxresp], 'ko', 'MarkerFaceColor', 'k'); hold on
plot(log(x),tp.thresh_resp.*ones(1,length(trl_4pulse)), 'r');
set(gca,'XTick',log(x)); logx2raw;   xlabel('Duration for 4 pulses (msec)');  widen; grid;
ylabel('Model Response @ Amplitude Threshold')
mse = sqrt(mean(([trl_4pulse.maxresp]-tp.thresh_resp).^2))./tp.thresh_resp;
disp(['pd50 mse = ', num2str(mse)]);

%%  'dur50' (figure 6b)
eid =strcmp(T.experiment,'dur50'); Tptdur = T(eid,:);
trl_ptdur= p2p_c.loop_convolve_model(tp, Tptdur);
tp.thresh_resp = mean([trl_ptdur.maxresp]);
disp(['50 Hz  vary dur thresh = ', num2str( tp.thresh_resp )])
if FITFLAG
    freeParams = {'thresh_resp', 'power'};
    tp = fit('p2p_c.loop_find_threshold',tp,freeParams,Tpd1);
end
disp(['fitted 50 Hz  vary dur thresh = ', num2str( tp.thresh_resp )])
disp(['fitted 50 Hz  power = ', num2str( tp.power )])
figure(6); clf; subplot(2,1,1)
x = Tptdur.dur*1000; y = Tptdur.amp; s = Tptdur.sem;
errorbar(log(x),y,s,'LineStyle','none','Color','r');
plot(log(x),y,'ro'); p2p_c.fillSymbols; hold on
set(gca,'XTick',log(x)); logx2raw; xlabel('Duration of 50Hz pulse train'); widen; grid ;title('4 pulses/trial');
title('50 Hz vary duration');
ylabel('Amplitude Threshold');
[err, thresh] = p2p_c.loop_find_threshold(tp, Tptdur);
plot(log(x),thresh,'k-');

subplot(2,1,2)
plot(log(x), [trl_ptdur.maxresp], 'ko', 'MarkerFaceColor', 'k'); hold on
plot(log(x),tp.thresh_resp.*ones(1,length(trl_ptdur)), 'r');
set(gca,'XTick',log(x)); logx2raw;   xlabel('Duration of 50Hz pulse train');  widen; grid;
ylabel('Model Response @ Amplitude Threshold')
mse = sqrt(mean(([trl_ptdur.maxresp]-tp.thresh_resp).^2))./tp.thresh_resp;
disp(['50 Hz  vary dur mse = ', num2str(mse)])
return
%% puzzling about why chronaxie too linear

tau = [.3 0.003];
amp = [1 100 ];
x = [.1:0.1:2]/1000;
clr =[1 0 0 ; 0 1 0]
for i = 1:length(tau)
    p.tau = tau(i);
    p.amp = amp(i);
    out(i, :)  = p2p_c.chronaxie(p, x);
    p
    plot(x*1000, out(i, :), 'Color', clr(i, :)); hold on
end

%
% if 0
%     chron = repmat([4:0.1:10]/1000, 10, 1);
%     rheo = 0;
%     x = [.1:0.1:2]/1000;
%     for i = 1:length(chron)
%         p.chron = chron(i);
%         p.rheo = rheo(i);
%         out(i, :)  = p2p_c.chronaxie(p, x);
%     end
%
%     for i = 1:length(chron)
%         plot(x*1000, out(i, :), 'Color', [.7 .7 .7]); hold on
%     end
%     plot(x*1000,min(out, [], 1), 'k-', 'LineWidth', 4); hold on
% end
%
%
%
