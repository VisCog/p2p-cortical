
%% fitPulseWidth.m
%
% Fits a variety of pulse width Data
%

figure(1)
clf
hold on

figure(2)
clf
hold on

tp = p2p_c.define_temporalparameters();

% begin by creating model response to our standard pulse
clear standard_trl;
standard_trl.pw = 0.001;
standard_trl.amp = 3;
standard_trl.dur = 1 ;
standard_trl.freq = 50;
standard_trl.simdur = 1; %sec
standard_trl = p2p_c.define_trial(tp,standard_trl);
standard_trl= p2p_c.convolve_model(tp, standard_trl);
tp.thresh_resp = standard_trl.maxresp;% use that as the threshold response

% predict smooth curve at 50Hz.  Should go through our standard pulse (amp
% = 3 at pw = 1msec)

pdList = exp(linspace(log(0.05), log(6.4), 40))/1000;
freq = standard_trl.freq*ones(size(pdList));
dur = standard_trl.dur*ones(size(pdList));

Tsim_50 = table(pdList', freq', dur');
Tsim_50.Properties.VariableNames = {'pw','freq','dur'};
[err, loop_trl_50] = p2p_c.loop_find_threshold(tp,Tsim_50);  %Tsim?

% Plot the curve in both figures 1 and 2
figure(1)
plot(log(pdList*1000),[loop_trl_50],'k-','LineWidth',2);

figure(2)
plot(log(pdList*1000),[loop_trl_50],'k-','LineWidth',2);

% Figure 2 will hold the predicted thresholds for each experiment scaled
% to match the standard stimulus curve.  The curves match perfectly - the
% shape of the curve does not depend on stimulus parameters dur or freq.

% Now on to the real experiments...

% paperList = {'Dobelle74'; 'Dobelle79'; 'Henderson79' ;'Girvin79'; 'Fernandez21';};

paperList = {'Dobelle74'; 'Dobelle79';'Girvin79'; 'Fernandez21';};

exList = {{'pd'},{'pd'},{'pd1','pd50'},{'pd'}};

legStr = {'Dobelle 1974, 50Hz 1 sec duration',...
    'Dobelle 1979, 50Hz 0.5 sec duration',...
    'Girvin 1979, 1 pulse 0.5 sec duration','Girvin 1979 50Hz, 0.5 sec duration',...
    'Fernandez 2021, 300Hz .167 sec duration'};

sz = 50;
%Data from Brindley et al. 1968 are excluded because measured thresholds were
%  unexpectedly non-monotonic as a function of frequency, suggesting an electronics issue
colList = {'b', 'r', 'c', 'm', 'g'};
xjitterFac = .2; yjitterFac = .2;
count =1;
for i = 1:length(paperList)
    for j=1:length(exList{i})
        disp(sprintf('Paper: %s, experiment: %s',paperList{i},exList{i}{j}))
        % reorganize data for ease of normalization
        T = readtable(['datasets/', paperList{i}, '_data.xlsx']);
        eid =strcmp(T.experiment,exList{i}{j});
        Texp = T(eid,:);
         [err, thresh] = p2p_c.loop_find_threshold(tp,Texp);
        
        Tsim = Texp;
        Tsim.freq = standard_trl.freq*ones(size(Tsim.freq));
        Tsim.dur = standard_trl.dur*ones(size(Tsim.dur));
        
        % get model response to standard stimulus at these pd's
        [err, standard_thresh] = p2p_c.loop_find_threshold(tp,Tsim);
        
        % regression to scale the actual data.
        scFac =  standard_thresh'/Texp.amp';
        
        y = Texp.amp*scFac;
        figure(1)
        h(count)= scatter(log(Texp.pw*1000) + xjitterFac.*(rand(size(y))-.5),y+yjitterFac.*(rand(size(y))-.5),sz);
        h(count).MarkerFaceColor =  colList{count};
        h(count).MarkerEdgeColor= 'none';
        h(count).MarkerFaceAlpha = .5;    
        h(count).MarkerEdgeAlpha = .5;
        drawnow
        
        figure(2)
        hold on
        scFac2 = standard_thresh'/thresh';
        fh= scatter(log(Texp.pw*1000),thresh*scFac2,sz);
           fh.MarkerFaceColor =  colList{count};
                fh.MarkerEdgeColor= 'none';
            fh.MarkerFaceAlpha = .5;    
         fh.MarkerEdgeAlpha = .5;
        drawnow
        
        count = count+1;
    end    
end

xtick = [.05,.1,.2,.4,.8,1.6,3.2,6.4];

figure(1)
set(gca,'XTick',log(xtick)); logx2raw; ylabel('Threshold'); widen; grid;
xlabel('Pulse Width (msec)'); ylabel('Amplitude @ Threshold');
legend(h,legStr,'Location','NorthEast');

figure(2)
set(gca,'XTick',log(xtick)); logx2raw; ylabel('Threshold'); widen; grid;
xlabel('Pulse Width (msec)'); ylabel('Amplitude @ Threshold');



