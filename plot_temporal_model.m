
%% fitPulseWidth.m
%
% Fits a variety of pulse width data, examining how threshold varies with
% pulse width
clear all; close all;

tp = p2p_c.define_temporalparameters();

% Calculate model response to a standard trial

clear trl;
trl.pw = 0.001;
trl.amp = 3;
trl.dur = .2 ;
trl.freq = 50;
trl.simdur = 1; %sec
t_crop = round((trl.dur*1.2)/tp.dt);
trl = p2p_c.define_trial(tp,trl);
trl= p2p_c.convolve_model(tp, trl);

subplot(5,1,1); 
 plot(trl.t(1:t_crop), trl.pt(1:t_crop), 'k');
set(gca, 'YLim', trl.amp*1.05*[-1 1]);
ylabel('Current Amplitude');
xlabel('Time');

subplot(5, 1, 2)
spikes = zeros(1, t_crop);
spikes(trl.spikeId) = trl.spikes_norefrac;
plot(trl.t(1:t_crop), spikes(1:t_crop), 'g');
ylabel('Spike no Refract');
xlabel('Time');

subplot(5, 1, 3)
spikes = zeros(1, t_crop);
spikes(trl.spikeId) = trl.spikes;
plot(trl.t(1:t_crop), spikes(1:t_crop), 'g');
ylabel('Spike with Refract');
xlabel('Time');

subplot(5, 1, 5)
plot(trl.t(1:tp.tSamp:end), trl.resp*1000, 'r');

figure(2)
x = 0:.1:20;
     y = p2p_c.nonlinearity(tp,x);
     plot(x, y, 'k-');
     xlabel('neural response')
     ylabel('compression');



return
tp.thresh_resp = trl.maxresp;% use that as the threshold response

subplot(5, 1,1); hold on



% calculate output first stage
Rtmp = 0;
ptid =[1  find(diff(trl.pt))+1 length(trl.pt)]; % times of event changes
wasRising = 0;
for i=1:max(ptid)
    if i >=ptid(2)
        ptid = ptid(2:end);
    end
    delta = trl.t(i)-trl.t(ptid(1));  % time since last 'event'
    %  Rtmp = Rtmp + tp.dt * (trl.pt(i)-Rtmp)/tp.tau1;
    Rtmp(i+1) = trl.pt(i)*tp.tau1*(1-exp(-delta/tp.tau1)) + ...
        Rtmp(i)*exp(-delta/tp.tau1);
    if Rtmp(i+1)<Rtmp(i) && wasRising
        spikeId(i) = 1; % check spike id identical in both loops IF CHECK
        wasRising = 0;  % no longer rising
    else
        wasRising =1;
    end
end

subplot(5,1,1); plot(trl.t(1:t_crop), 1000*Rtmp(1:t_crop), 'r'); hold on
ylabel('Tau1 response');
xlabel('Time');

% for i=1:length(x)
%     h(i)= scatter(x{i}+sd*randn(size(x{i})),y2{i},sz , 'MarkerFaceColor', colList{i},...
%         'MarkerEdgeColor','k','MarkerFaceAlpha',alpha);
% end
%
% set(gca,'XTick',log(xtick)); logx2raw; ylabel('Threshold'); widen; grid;
% xlabel('Pulse Width (msec)'); ylabel('Amplitude @ Threshold');
% legend(h,legStr,'Location','NorthEast');
% set(gca,'YLim',[0,20]);
% set(gca,'FontSize',fontSize);
% title('Scaled model predictions');
%
%



