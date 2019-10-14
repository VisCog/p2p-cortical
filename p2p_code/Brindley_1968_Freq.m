% Brindley_Lewin_1968
clear all
addpath(genpath('C:\Users\Ione Fine\Documents\code\UWToolbox\UWToolbox\'));
rng(1)
c.efthr = 0.05;
v.drawthr = .015;
tp.scFac = 1; % this scales the Brindley numbers so that with similar pulse parameter
% you get similar thresholds

% electrode is square but we are ignoring this since we're just looking at
% threshold data

tp.model = 'normcdf';
tp = p2p_c.define_temporalparameters(tp);

%% plot chronaxie data

% plot first trial as a phosphene
brindley = Brindley_getData('frequency');
trl.dur = 1;
for pp=1:length(brindley)
    disp(['finding threshold for trial ', num2str(pp), ' out of ', num2str(length(brindley))]);
 
    trl = brindley(pp);
    trl = p2p_c.define_trial(tp,trl);

    % trl = p2p_c.p2p_finite_element(tp, trl);
     trl = findthreshold(trl, tp, .015);
    brindley(pp).resp_max = max(trl.R2(:));
    brindley(pp).fitted_amp = trl.amp;
end

figure(1); clf
plot(cat(1, brindley(:).freq), cat(1, brindley(:).amp), 'k-o', 'MarkerFaceColor', 'k', 'MarkerSize', 10); hold on
plot(cat(1, brindley(:).freq), cat(1, brindley(:).fitted_amp)/28, 'bo'); hold on


function trl = findthreshold(trl, tp, thr)
lo = 1; hi = 400;

for i=1:10
    
    mid = (hi+lo)/2;
    trl.amp = mid;
    trl = p2p_c.define_trial(tp,trl);
    trl = p2p_c.p2p_finite_element(tp, trl);

    if max(trl.resp(:)) > thr
        hi = mid;
    else
        lo = mid;
    end
    if abs(max(trl.resp(:))-thr)<0.0001
        i=1000;
    end
    disp(['amp = ', num2str(trl.amp), ' ', '/ resp = ', num2str(max(trl.resp(:)))]);
end
trl.amp = (hi+lo)/2;

end
        
function alltrl = Brindley_getData(expname)

if strcmp(expname, 'chronaxie')
    mat(1,:)= [8 9 9 10 13 16 19 25 28 36 56]; % threshold volts
    mat(2,:) = [ 1000 600 400 300 200 100 60 40 30 20 10] * 10.^-6; % pulse width (s)
    mat(3,:) = 30; % frequency of 30 Hz
elseif strcmp(expname, 'frequency')
    mat(1,:) = [29 27 21 21 25 35 37 39 35 29]; % threshold current
    mat(2,:)= 30 * 10.^-6; % pulse duration 30 microsec, translate to s
    mat(3,:)= [  25 50 100 160 250 400 630 1000 1600 4000]; % freq
end
for i = 1:size(mat,2)
    alltrl(i).amp = 1000 * mat(1,i)/3000; % estimated resistance of 3000
    alltrl(i).pw = mat(2,i);
    alltrl(i).freq = mat(3,i);
    alltrl(i).dur = .25; % don't actually know the duration, not in the paper!
end
end