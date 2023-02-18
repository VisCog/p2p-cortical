%Fit_Girvin_Data_pd.m
% fits Girvin data using a probability summation model
clear all; close all

T = readtable('datasets/Girvin79_data.xlsx');
experimentList = {'pd1','pd50','freq','dur','dur50'}; %list of sets to fit.  Each set gets it's own thresh_resp
xLabelList = {'Pulse Duration (ms)', 'Pulse Duration (ms)', 'Frequency (Hz)' , 'Frequency (Hz)', '# Pulses'};
titleList = {'Single Pulse: pulse width varies', '50 Hz Pulse Train: pulse width varies', ...
    '0.5s Pulse Train (0.25ms PD): frequency varies','4 Pulses: frequency varies','50 Hz pulse train: # pulses (dur) varies'};
% pd1, pulse duration varies, 1 pulse (fig 4b)
% pd50, pulse duration varies, 50Hz (fig 4a)
% freq, 05s pulse train, frequency varies (fig 5)
% dur, 4 pulses, frequency (duration to do all 4 pulses) varies (6a)
% dur 50, 50 Hz stimulation, pulse train duration varies (i.e. number of
% pulses) (6b)

tp.tSamp = 1000;
tp.dt  = 1e-6;
tp.nreps = 12; % the duration of the binary search to find threshold
tp.tau1 = [ 1]; % fixed based on Nowak and Bullier, 1998
tp.tau1Fac = [ 1]; %.relative sensitivity of the two tau1 if using more than one
tp.refractoryPeriod = [85 10];
tp.ncascades  = 3;
tp.win = .5/tp.tSamp;
tp.probsumflag = 1;
tp.gammaflag = 0;
tp.refractoryflag = 0;
tp.tau2 =12.5/1000; % used if second stage is gamma

tp.model = 'normcdf'; % weirdly, this non-linearity used in different places in the model for gamma and probsum
tp.scFac = 1; tp.sigma = .15; tp.mean =0.6; tp.power =1; %only some used, depends on tp.model
tp.thresh_resp = 0.75; % 2 alt FC detection threshold for Girvin, fixed for probsum second stage

% plot nonlinearity
input= linspace(0, tp.mean*3, 1000);
figure(11); clf; y = p2p_c.nonlinearity(tp,input);
plot(input, input, 'k'); hold on;
plot(input,y, 'r');  
1-exp(log(.75)/4)
disp(p2p_c.nonlinearity(tp,0.4214 )); % this is the asymptotic spike value for two spikes spike in dur50 single pulse experiment
FITFLAG = 0;

for ex = [4]
    eid =strcmp(T.experiment,experimentList{ex}); Texp = T(eid,:);
    if strcmp(xLabelList{ex}, 'Pulse Duration (ms)')
        x = Texp.pw*1000;
    elseif strcmp(xLabelList{ex}, 'Frequency (Hz)')
        x = Texp.freq;
    elseif strcmp(xLabelList{ex}, '# Pulses')
        x = ceil([Texp.freq].*[Texp.dur]); % number of pulses
    else errordlg(['Experimental variation not specified for experiment', num2str(ex)]);
    end
    figure(ex); clf; subplot(2, 1, 1); hold on
    y = Texp.amp; s = Texp.sem;
    errorbar(log(x),y,s,'LineStyle','none','Color','r');
    plot(log(x),y,'ro'); p2p_c.fillSymbols;
    set(gca,'XTick',log(x)); logx2raw; ylabel('Threshold'); widen; grid;
    xlabel(xLabelList{ex}); ylabel('Amplitude @ Threshold');
    if FITFLAG
        freeParams = {'tau1'};
        tp = fit('p2p_c.loop_find_threshold',tp,freeParams,Texp);
    end

    subplot(2,1,2)
    trl= p2p_c.loop_convolve_model(tp, Texp);
    plot(log(x), [trl.pd], 'ko', 'MarkerFaceColor', 'k'); hold on
    plot(log(x),tp.thresh_resp.*ones(1,length([trl])), 'r');
    set(gca,'XTick',log(x)); logx2raw; xlabel(xLabelList{ex});  ylabel('thresh'); widen; grid;
    ylabel('p(D) @ Amplitude Threshold');

    subplot(2,1,1)
    [err, thresh] = p2p_c.loop_find_threshold(tp,Texp); thresh = reshape(thresh, size(y)); % find the amplitude which results in 0.75 detection probability
    plot(log(x),thresh,'k-');
    % legend(ph, {'Data: Amp @ Threshold','Model Prediction'});
    mse = sqrt(mean((thresh-y).^2));
    disp(['mse = ', num2str(mse)]);
    title(titleList{ex});
end
return
% %% puzzling about why chronaxie too linear
%
% tau = [.3 0.003];
% amp = [1 100 ];
% x = [.1:0.1:2]/1000;
% clr =[1 0 0 ; 0 1 0]
% for i = 1:length(tau)
%     p.tau = tau(i);
%     p.amp = amp(i);
%     out(i, :)  = p2p_c.chronaxie(p, x);
%     p
%     plot(x*1000, out(i, :), 'Color', clr(i, :)); hold on
% end
%
% figure(10); clf
% clr =[0 0 1 ; 0 1 0];
% pw = [0.0001:0.0001:001];
% for i = 1:length(tp.tau1)
%     f.tau = tp.tau1(i); f.amp = tp.tau1Fac(i)/1000;
%     out(i, :)  = p2p_c.chronaxie(f, pw);
% end
% for i = 1:length(tp.tau1)
%     plot(log(pw*1000), out(i, :), '--', 'Color', clr(i, :)); hold on
% end
% plot(log(pw*1000), min(out, [], 1), '--', 'Color', [.5 .5 .5], 'LineWidth', 2); hold on
