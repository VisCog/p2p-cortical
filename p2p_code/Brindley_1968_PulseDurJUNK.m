% Brindley_Lewin_1968
% Simulates chronaxie (strength-duration) data from Brindley and Lewin
% Table 1

clear all
PARALLEL = 0;

c.efthr = 0.05; % lower limit of electrical stimulation for which the phosphene is calculated
v.drawthr = .015; % visibility threshold
tp.scFac = 1; % overall scaling factor for stimulation strength
% electrode is square but we are ignoring space since we're just looking at
% threshold data

brindley = Brindley_getData('chronaxie');
tp.tau1 = 0.012* 10^-3; % this was manually fitted to the Brindley data
tp = p2p_c.define_temporalparameters(tp);
trl = [];
numCores = 0;

if PARALLEL
    numCores = feature('numcores');
    p = parpool(numCores);
end

trl.pw = 0.5;
trl.dur = 1;
 trl = p2p_c.define_trial(tp, trl);
 trl.pt = 100 *ones(size(trl.pt));
trl = p2p_c.p2p_finite_element(tp, trl);
   return
%parfor (tt=1:length(brindley), numCores * PARALLEL)
   for tt=1
    disp(['finding threshold for trial ', num2str(tt), ' out of ', num2str(length(brindley))]);
    trl = p2p_c.define_trial(tp,brindley(tt));
    trl = findthreshold(trl, tp, v.drawthr);
    tmp(tt).fitted_amp = trl.amp;
end

figure(1); clf
plot(cat(1, brindley(:).pw), cat(1, brindley(:).amp), 'k-o', 'MarkerFaceColor', 'k'); hold on
plot(cat(1, brindley(:).pw), cat(1, tmp(:).fitted_amp)/23, 'r-s','MarkerFaceColor', 'r');

function trl = findthreshold(trl, tp, thr)
tol = 0.001;
lo = 1; hi = 20000;
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
    if abs(max(trl.resp(:))-thr) < tol % close enough
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