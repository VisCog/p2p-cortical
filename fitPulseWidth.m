
%% fitPulseWidth.m
%
% Fits a variety of pulse width Data
%

tp = p2p_c.define_temporalparameters();
expList = {'Dobelle74'; 'Dobelle79'; 'Henderson79' ;'Girvin79'; 'Fernandez21';};
%Data from Brindley et al. 1968 are excluded because measured thresholds were
%  unexpectedly non-monotonic as a function of frequency, suggesting an electronics issue
colList = {'b', 'r', 'c', 'm', 'g'};
pdList = exp(linspace(log(0.01), log(2), 8)); %

% begin by creating model response to our standard pulse
clear trl;  trl.pw = 0.001;   trl.amp = 3;    trl.dur = 1 ;  trl.freq = 50;   trl.simdur = 1; %sec
trl = p2p_c.define_trial(tp,trl);
trl= p2p_c.convolve_model(tp, trl);
tp.thresh_resp = trl.maxresp;% use that as the threshold response

for ex = 1:length(expList)
    % reorganize data for ease of normalization
        T = readtable(['datasets/', expList{ex}, '_data.xlsx']);
        eid =strcmp(T.experiment,'pd'); Texp = T(eid,:);
        [err, thresh] = p2p_c.loop_find_threshold(tp,Texp); % get the model response for these pulse trains

        x = Texp.pw*1000;y = Texp.amp; 
        % now linearly scale the data, to the same sensitivity range as the model
        k = mean(y./thresh); y = y*k; % rescale y to right units
        for i =1:length(x)
            plot(log(x(i)),y(i),'o', 'MarkerFaceColor', colList{ex}, 'MarkerEdgeColor','none'); hold on
        end
        plot(log(x),[thresh],'*', 'Color', colList{ex}); hold on % plot simulation values
end
pw = exp(linspace(log(0.01), log(2), 8));freq = 50*ones(size(pw)); dur = ones(size(pw));
Tsim = table(pw, freq, dur);
[loop_trl] = p2p_c.loop_find_threshold(tp,Tsim);
plot(log(pd),[loop_trl(:).maxresp],'k-');
set(gca,'XTick',log(pdList)); logx2raw; ylabel('Threshold'); widen; grid;
xlabel('Pulse Width'); ylabel('Amplitude @ Threshold');
return
%
%     el = eval(['p2p_psy.', expList{x}, '_getDataPW();']);
%     for e = 1:length(el)
%         for p = 1:length(pwList)
%             ind = find([el(e).trl(:).pw]*1000==pwList(p));
%             if ~isempty(ind)
%                 data.amp(ct, p)=el(e).trl(ind).amp_n;
%                 data.freq(ct) = el(e).trl(ind).freq;
%             else
%                 data.amp(ct, p) = NaN;
%             end
%             data.col(ct) = colList(x);
%         end
%         ct = ct + 1;
%     end
% 
% 
% %% plot the data
% for x =1:size(data.amp, 1)
%     plot(pdList, data.amp(x,:), 'o', 'MarkerFaceColor', data.col{x}, 'MarkerEdgeColor', 'none'); hold on
% end
% xlabel('pulse width');
% ylabel('threshold')
% 
% %% model pw
% tp.model = 'simpleleakyintegrator';
% tp.tau1 = 8.899999999999999e-05;
% tp.scaleR1 = 1; tp.scaleR4 = 1;
% tp = p2p_c.define_temporalparameters(tp);
% 
% trl.freq = 50;
% trl.dur = 1000*10^-3;
% trl.pw = 1000* 10.^-6;
% trl.amp = 1;
% 
% trl = p2p_c.define_trial(tp,trl);
% trl = p2p_c.finite_element(tp, trl);
% tp.scaleAmp = 1;
% % p2p_c.find_scaleAmp(trl, tp, fitParams);
% % finds the scale factor such that a current amplitude of fitParams.thr results in a
% % resp value of fitParams.thr.
% 
% fitParams.tol = 0.05; fitParams.lo = 0; fitParams.hi = 75;
% fitParams.thr = .5; fitParams.nreps = 10;
% fList = 50; %[1 30 50 300];
% for f = 1:length(fList)
%     trl.freq = fList(f);
%     for p = 1:length(pdList)
%         disp(['finding pw threshold ', num2str(p), ' out of ', num2str(length(pdList))]);
%         disp(['pulse width is ', num2str(pdList(p))]);
%         trl.pw = pdList(p)/1000;
%         trl = p2p_c.define_trial(tp,trl);
%         data.model(f,p) = p2p_c.find_threshold(trl, tp, fitParams);
%     end
% end
% 
% %% plot the pw model
% for f = 1:length(fList)
%     plot(pdList, data.model(f,:), '-'); hold on
% end
