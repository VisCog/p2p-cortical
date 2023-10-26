% Plot_Temporal_Model.m
%
% Plots each stage of the temporal model
%
% written IF & GMB
%
% 25/02/2023 moved into clean folder (IF)
% 28/02/2023 (gmb) added leaky integrator response (R1)

tp = p2p_c.define_temporalparameters();
% Calculate model response to a standard trial
clear trl;
trl.pw = 0.002;
trl.amp = 3;
trl.dur = .05 ;
trl.freq = 75;
trl.simdur = 1.5; %sec

t_crop = round((trl.dur*1.2)/tp.dt);
trl = p2p_c.define_trial(tp,trl);
trl= p2p_c.convolve_model(tp, trl);

% R1 is not needed to compute spikes, but here it is for plotting:

ptid = find(diff(trl.pt))+1;
ptid = [ptid,ptid(end)+2000];  % recovery from last pulse
R = zeros(1,length(ptid));

% Loop through the events, calculating R1 at the end of the event
% and add impulse responses when R1 peaks and is after the refractory period.
spikeId = zeros(size(R));
trl.R1 = zeros(size(trl.t));  % will hold continuous R1 (leaky integrator) response
lastSpikeTime = -1;  % Keep track of time of last spike for refractory period.  Start 'fresh'
for i=1:(length(ptid)-1)
    tNow = trl.t(ptid(i+1));
    delta = trl.t(ptid(i+1))-trl.t(ptid(i));  % time since last 'event'
    % Closed form solution to leaky integrator that predicts
    % R(i+1) from R(i), delta and tau1:
    
    R(i+1) = exp(-delta/tp.tau1)*[R(i)-trl.pt(ptid(i))*tp.tau1] + trl.pt(ptid(i))*tp.tau1;
    ti = 1:(ptid(i+1)-ptid(i));
    trl.R1(ptid(i)+ti) =  exp(-ti*tp.dt/tp.tau1)*[R(i)-trl.pt(ptid(i))*tp.tau1] + trl.pt(ptid(i))*tp.tau1;
end
tplot = 1000*(trl.t(1:t_crop)-trl.lag);
xlim = 1000*[-trl.lag,trl.dur];

figure(1)
clf
% pulse train
subplot(4,1,1);
plot(tplot, trl.pt(1:t_crop), 'k-');
set(gca, 'YLim', trl.amp*1.2*[-1 1]);
set(gca,'XLim',xlim);
ylabel('Current Amplitude');
xlabel('Time (msec)');
title('Pulse Train');


% R1
subplot(4,1,2)
plot(tplot,1000*trl.R1(1:t_crop),'k-');
ylabel('R_1');
set(gca,'XLim',xlim);
set(gca,'YLim',[-1,1]);

xlabel('Time (msec)');
title('Leaky Integrator Response');

% spikes without refractory period
subplot(4, 1, 3)
spikes = zeros(1, t_crop);
spikes(trl.spikeWhen) = trl.spikes_norefrac;
plot(tplot, spikes(1:t_crop), '--', 'Color', [.7 .7 .7], 'LineWidth', 2); hold on
set(gca,'XLim',xlim);

ylabel('Spikes');
xlabel('Time (msec)');
title('Spikes with No Refractory Period');


% spikes with refractory period
subplot(4, 1, 3)
spikes = zeros(1, t_crop);
spikes(trl.spikeWhen) = trl.spikes;
plot(tplot, spikes(1:t_crop), 'k', 'LineWidth',2);
set(gca,'XLim',xlim);
ylabel('Spikes')
xlabel('Time (msec)');
title('Spikes with Refractory Period');

% response before nonliearity

% subplot(5, 1, 5)
% plot(trl.t(1:tp.tSamp:end), trl.resp_nononlin, 'r');
% ylabel('R_2');
% xlabel('Time (sec)');
% title('Response before nonlinearity')

% response after nonliearity
subplot(5, 1, 5)
hold on
plot([0,trl.dur],[0,0],'k-','LineWidth',4);
plot(trl.t(1:tp.tSamp:end)-trl.lag, trl.resp, 'k-');
set(gca,'YLim',[0,10]);
set(gca,'XLim',[-trl.lag,.75]);
set(gca,'YTick',[0:2:8]);
ylabel('Final Response');
xlabel('Time (sec)');
title('Final Response')



figure(2); clf
subplot(1,3,1)
h = p2p_c.gamma(tp.ncascades,tp.tau2,trl.t);
idx = trl.t(trl.t<1);
plot(trl.t(1:tp.tSamp:800000), h(1:tp.tSamp:800000), 'k-');
xlabel('Time (sec)');
ylabel('Gamma filter response');
title('Gamma filter');


interspike = 1./linspace(1000, 1, 1000);
refrac =(1-exp(-tp.refrac*(interspike+tp.delta)));
subplot(1,3,2)
plot(log10(interspike), (1-refrac), 'k-');
set(gca, 'XTick', log10([0.001 0.01  0.1 1]));
set(gca, 'XTickLabel', [0.001 0.01  0.1 1]*1000);
set(gca, 'XLim', log10([ 0.001 1]));
xlabel('Time since last spike (sec)');
ylabel('Refract attenuation');
title('Refractory Period filter');
logx2raw(10); 

subplot(1, 3,3)
x = 0:.1:100;
y = p2p_c.nonlinearity(tp,x);
plot(x, y, 'k-'); hold on
plot(x, x, '--', 'Color', [.5 .5 .5]);
xlabel('Neural response')
ylabel('Neural Compression');
title('Compressive Nonlinearity');
set(gca, 'YLim', [0 22]);








