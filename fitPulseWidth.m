
%% fitPulseWidth.m
%
% Fits a collection of temporal psychophysics data
%
% clear trl
% trl.amp = 1000; trl.pw = 1000* 10.^-6;
% trl.freq =  50; trl.dur =  0.5000;

% find the sensitivity scale factor so a 'standard pulse'
% has a resp of 1
% with a current input of 1000


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PULSE WIDTH
% load and plot pulse width data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


expList = {'Dobelle74'; 'Dobelle79'; 'Henderson1979' ;'Girvin79'; 'Fernandez2021';};
colList = {'b', 'r', 'c', 'm', 'g'};
pwList =  [2 1.5 1 0.8 0.75 0.5 0.4 0.25 0.17 0.125 0.1 0.062];

% reorganize data for ease of normalization
ct = 1;
for x = 1:length(expList)
    el = eval(['p2p_psy.', expList{x}, '_getDataPW();']);
    for e = 1:length(el)
        for p = 1:length(pwList)
            ind = find([el(e).trl(:).pw]*1000==pwList(p));
            if ~isempty(ind)
                data.amp(ct, p)=el(e).trl(ind).amp_n;
                data.freq(ct) = el(e).trl(ind).freq;
            else
                data.amp(ct, p) = NaN;
            end
            data.col(ct) = colList(x);
        end
        ct = ct + 1;
    end
end

%% plot the data
for x =1:size(data.amp, 1)
    plot(pwList, data.amp(x,:), 'o', 'MarkerFaceColor', data.col{x}, 'MarkerEdgeColor', 'none'); hold on
end
xlabel('pulse width');
ylabel('threshold')

%% model pw
tp = p2p_c.define_temporalparameters();
tp.tau1 = 8.899999999999999e-05;
tp.model = 'simpleleakyintegrator';
trl.freq = 50;
trl.dur = 1000*10^-3;
trl.pw = 1000* 10.^-6;
trl.amp = 1;
tp.scaleAmp = 1;
trl = p2p_c.define_trial(tp,trl);
trl = p2p_c.finite_element(tp, trl);
tp.scaleAmp = 1;
% p2p_c.find_scaleAmp(trl, tp, fitParams);
% finds the scale factor such that a current amplitude of fitParams.thr results in a
% resp value of fitParams.thr.

fitParams.tol = 0.05; fitParams.lo = 0; fitParams.hi = 75;
fitParams.thr = .5; fitParams.nreps = 10;
fList = 50; %[1 30 50 300];
for f = 1:length(fList)
    trl.freq = fList(f);
    for p = 1:length(pwList)
        disp(['finding pw threshold ', num2str(p), ' out of ', num2str(length(pwList))]);
        disp(['pulse width is ', num2str(pwList(p))]);
        trl.pw = pwList(p)/1000;
        trl = p2p_c.define_trial(tp,trl);
        data.model(f,p) = p2p_c.find_threshold(trl, tp, fitParams);
    end
end

%% plot the pw model
for f = 1:length(fList)
    plot(pwList, data.model(f,:), '-'); hold on
end
