% Troy_hypothetical.m

clear all; close all; clear mex

rng(1)  % fix the random number generator. This affects the ocular dominance/orientation maps
c.efthr = 0.85; % what magnitude of electric field goes through the model, just a speed thing so one doesn't bother processing non-active regions of cortex
tp.scFac = 1;  % scaling of the strength of the pulse train
v.drawthr = .025; % what phosphene strength is visible, threshold of visibility/oval drawing
filename = 'Troyk_Snellen_';

% set up the locations in terms of their positions in retinotopic space
r = .05;
x_off = 1.5;
y_off = 0;


v.e(1).x =1;
v.e(1).y = 1;
[ang,ecc] = cart2pol( v.e(1).x, v.e(1).y);
v.e(1).ang = ang*180/pi;
v.e(1).ecc = ecc;
c.e(1).radius = 1/5000;


% define cortex & retina
c.cortexSize = [40,50]; % degrees top to bottom, degrees LR, divide by 2 to get the actual mm that are useful
c.pixpermm = 30; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased
v.retinaSize = [8,8]; v.pixperdeg = 15;  %visual field map size and samping
c = p2p_c.define_cortex(c); % define the properties of the cortical map
c = p2p_c.define_electrodes(c, v); % defines properties for each electrode in retinal space
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
c = p2p_c.generate_ef(c); % generate map of the electric field on cortical surface

figure(1); clf
ef = c.e(1).ef;
for e=2:length(c.e)
    ef = ef + c.e(e).ef;
end
p2p_c.plotcortgrid(256.*(ef./max(ef(:))), c, gray(256), 1,['title(''electric field'')'])

tp = p2p_c.define_temporalparameters(); % define the temporal model
 c.rftype = 'rf';
v = p2p_c.generate_rfmap(c, v); 

%% generate percepts for a pulse train
clist = parula(length(v.e));
trl.lag = 100*10^-3;
trl.expname = 'Troyk_hypothetical';

trl.e = 1;
trl = p2p_c.define_trial(tp,trl);
[trl,v] = p2p_c.generate_phosphene(v, tp, trl);

%% make visual image

figure(3); clf
img = zeros(size(squeeze(trl.maxphos(:, :, 1))));

for ii=1:length(v.e)
    if isempty(find(isnan(v.e(ii).rfmap)))
        img = img + v.e(ii).rfmap(:, :, 1);
    end
end
subplot(1,2,1)
p2p_c.plotretgrid((img./max(img(:)))*256, v, gray(256), 3,['';]);

img = zeros(size(squeeze(trl.maxphos(:, :, 1))));

for ii=1:length(v.e)
    if isempty(find(isnan(v.e(ii).rfmap)))
        img = img + v.e(ii).rfmap(:, :, 2);
    end
end
subplot(1,2,2)
p2p_c.plotretgrid((img./max(img(:)))*256, v, gray(256), 3,['';]);





