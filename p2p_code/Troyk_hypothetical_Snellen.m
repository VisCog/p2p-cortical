% Troy_hypothetical.m

clear all; close all; clear mex

rng(1)  % fix the random number generator. This affects the ocular dominance/orientation maps
c.efthr = 0.05; % what magnitude of electric field goes through the model, just a speed thing so one doesn't bother processing non-active regions of cortex
tp.scFac = 1;  % scaling of the strength of the pulse train
v.drawthr = .015; % what phosphene strength is visible, threshold of visibility/oval drawing
filename = 'Troyk_Snellen_';

% set up the locations in terms of their positions in retinotopic space
r = 2;
x_off = 3;
y_off = 2;
sampFac = .2;
x = -r+x_off:sampFac:r+x_off
for i=1:length(x)
    y(i, 1)=sqrt(r^2-((x(i)-x_off).^2))+y_off;
    y(i,2)=-sqrt(r^2-((x(i)-x_off).^2))+y_off;
end

x = repmat(x',2,1); y = cat(1, y(:, 1), y(:, 2));
y = round(y*5, 0)/5;
xy = unique([x y], 'rows');
plot(xy(:, 1),xy(:, 2), 'o'); hold on
for i = 1:length(xy)
    [theta, v.e(i).ecc] = cart2pol(xy(i, 1), xy(i, 2));
    v.e(i).ang = theta*180/pi;
    v.e(i).x = xy(i, 1);
    v.e(i).y = xy(i, 2);
    c.e(i).radius = 1/1000;
end

% define cortex & retina
c.cortexSize = [40,60]; % degrees top to bottom, degrees LR, divide by 2 to get the actual mm that are useful
c.pixpermm = 6; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased
v.retinaSize = [30,30]; v.pixperdeg = 5;  %visual field map size and samping
c = p2p_c.define_cortex(c); % define the properties of the cortical map
c = p2p_c.define_electrodes(c, v); % defines properties for each electrode in retinal space
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
c = p2p_c.generate_ef(c); % generate map of the electric field on cortical surface

tp = p2p_c.define_temporalparameters(); % define the temporal model
v = p2p_c.generate_rfmap(c, v); 

%% generate percepts for a pulse train
clist = parula(length(v.e));
    trl.lag = 100*10^-3;
trl.expname = 'Troyk_hypothetical';

trl.e = 1;
trl = p2p_c.define_trial(tp,trl);
[trl,v] = p2p_c.generate_phosphene(v, tp, trl);

%% make movie
close all
figure(1); clf
vid = VideoWriter([[filename], '.avi']);
vid.FrameRate = 30;
open(vid);
itpl= round(linspace(1,length(trl.pt), 30*trl.trialdur));
img = zeros(size(squeeze(trl.maxphos(:, :, 1))));
for ii=1:length(v.e)
    if isempty(find(isnan(v.e(ii).rfmap)))
        img = img + (v.e(ii).rfmap(:, :, 1).*max(trl.resp));
    end
end
p2p_c.plotretgrid(img*5, v, gray(256), 1,['';]);




