% Brindley_Lewin_1968
clear all
addpath(genpath('C:\Users\Ione Fine\Documents\code\UWToolbox\UWToolbox\'));
rng(1)
c.efthr = 0.05;
v.drawthr = .015;
v.e.ecc = 5; % electrode 19
v.e.ang = 20;
tp.scFac = 45; % this scales the Brindley numbers so that with similar pulse parameter
% you get similar thresholds
c.e.radius = 0.4; % 8mm length of side
% electrode is square but we are ignoring this since we're just looking at
% threshold data

c.cortexSize = [80,100]; c.pixpermm = 6; c.efthr = .1;
v.retinaSize = [60,60];v.pixperdeg = 5;  %10

c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
c = p2p_c.define_electrodes(c, v);

tp = p2p_c.define_temporalparameters(tp);
c = p2p_c.generate_ef(c);
v = p2p_c.generate_rfmap(c, v);

%% plot chronaxie data

% plot first trial as a phosphene
brindley = Brindley_getData('chronaxie');

trl = p2p_c.define_trial(tp,brindley(1));
trl = p2p_c.generate_phosphene(v, tp, trl);
figure(1); clf; figure(2); clf

brindley = brindley(6);% , brindley([6, 11]);
p2p_c.plotretgrid(trl.maxphos(:, :, 1)*500, v, gray(256),1,[ 'title(''phosphene'')';]); drawnow; 
  ampVals = [3];%  40]; %6.3333 8.3333]
  for aa = 1:length(ampVals)
     for pp=1:length(brindley)
         disp(['finding threshold for trial ', num2str(pp), ' out of ', num2str(length(brindley))]);
         
         trl = brindley(pp);
         %trl.amp = ampVals(aa);
         trl = p2p_c.define_trial(tp,trl);
         tp.model = 'normcdf';
     
      %    trl = p2p_c.p2p_finite_element(tp, trl);
        % trl = p2p_c.generate_phosphene(v, tp, trl);
         
       %  figure(2);
         
%          alltrl(aa, pp).resp_max = max(trl.resp);
%          alltrl(aa, pp).R1 = max(trl.R1);
%          alltrl(aa, pp).R2 = max(trl.R2);
%          
%          alltrl(aa, pp).R3 = max(trl.R3);
%          alltrl(aa, pp).R4 = max(trl.resp);
%                subplot(4, 1, 1)
%                plot(trl.R1); hold on
%                  subplot(4, 1, 2)
%                plot(trl.R2); hold on
%                  subplot(4, 1, 3)
%                plot(trl.R3); hold on
%                subplot(4, 1, 4)
%                plot(trl.resp); hold on
%             text(12000, trl.resp(12000), num2str(trl.pw*1000));
%          
                trl = p2p_c.findthreshold(tp, v, trl);
                 alltrl(aa,pp).thresh_amp = trl.thresh_amp;
     end
 end

 
%% more plot it
figure(3); clf
clr = hsv(9);

for aa = 1:length(ampVals)
subplot(4,1,1)
plot(cat(1,brindley.pw), cat(1,(alltrl(aa, :).R1)), '-s', ...
    'MarkerEdgeColor', clr(aa, :), 'MarkerFaceColor', clr(aa, :)); hold on

subplot(4,1,2)
plot(cat(1,brindley.pw), cat(1,(alltrl(aa, :).R2)), '-s', ...
    'MarkerEdgeColor', clr(aa, :), 'MarkerFaceColor', clr(aa, :)); hold on

subplot(4,1,3)
plot(cat(1,brindley.pw), cat(1,(alltrl(aa, :).R3)), '-s', ...
    'MarkerEdgeColor', clr(aa, :), 'MarkerFaceColor', clr(aa, :)); hold on
subplot(4,1,4)
plot(cat(1,brindley.pw), cat(1,(alltrl(aa, :).R4)), '-s', ...
    'MarkerEdgeColor', clr(aa, :), 'MarkerFaceColor', clr(aa, :)); hold on
plot(cat(1,brindley.pw),ones(size(cat(1,brindley.pw))).*v.drawthr, 'k--');
xlabel('pw'); 
end
return
%% plot it
subplot(1,2,1)
plot(cat(1,brindley.pw), cat(1,brindley.amp), 's','MarkerEdgeColor', ...
    'k', 'MarkerFaceColor', 'k'); hold on
plot(cat(1,brindley.pw), cat(1,alltrl(aa, :).thresh_amp), ...
    's','MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); hold on



% ylabel('amplitude at threshold');
% subplot(1,2,2)
% plot(cat(1,alltrl(6:end).pw), cat(1, alltrl(6:end).resp_max), 's-','MarkerEdgeColor', clr, 'MarkerFaceColor', clr); hold on
% xlabel('pw'); ylabel('model maximum response at threshold');
return

for i = 1:length(alltrl)
    tmp = p2p_c.define_trial(tp,alltrl(i));
    tmp = p2p_c.generate_phosphene(v, tp, tmp);
end
xlabel('Duration (microsec)')

%% plot frequency data
alltrl = Brindley_getData('frequency');
figure(2); subplot(1,2,2)
plot(cat(1,alltrl.freq), cat(1,alltrl.amp),'s','MarkerEdgeColor', clr, 'MarkerFaceColor', clr); hold on
xlabel('Frequency (Hz)')

function alltrl = Brindley_getData(expname)

if strcmp(expname, 'chronaxie')
    mat(1,:)= [8 9 9 19 13 16 19 25 28 36 56]; % threshold volts
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