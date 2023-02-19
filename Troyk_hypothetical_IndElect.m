% Troy_hypothetical.m

clear all; close all; clear mex

rng(1171964)  % fix the random number generator. This affects the ocular dominance/orientation maps
c.efthr = 0.05; % what magnitude of electric field goes through the model, just a speed thing so one doesn't bother processing non-active regions of cortex
v.drawthr = 1; % what phosphene strength is visible, threshold of visibility/oval drawing
filename = 'Troyk_Snellen_';

% define cortex & retina
c.cortexSize = [40,50]; % degrees top to bottom, degrees LR, divide by 2 to get the actual mm that are useful
c.cortexCenter = [30,0];
c.pixpermm = 30; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased

c = p2p_c.define_cortex(c); % define the properties of the cortical map
% set up the electrode locations in terms of their positions in the teeny
% array
e = create_WFMA();
c.e = e;

% transform to visual space
v.retinaSize = [8,8]; v.pixperdeg = 15;  %visual field map size and samping
v = p2p_c.define_visualmap(v); % defines the visual map
v = p2p_c.c2v_define_electrodes(c,v); % convert electrode locations from cortex to visual space


tp = p2p_c.define_temporalparameters();

return
r = .05;
x_off = 1.5;
y_off = 0;


v.e(1).x =1;
v.e(1).y = 1;
[ang,ecc] = cart2pol( v.e(1).x, v.e(1).y);
v.e(1).ang = ang*180/pi;
v.e(1).ecc = ecc;
c.e(1).radius = 1/5000;





c = p2p_c.define_electrodes(c, v); % defines properties for each electrode in retinal space

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

function e = create_WFMA()
%Troyk, P., Bredeson, S., Cogan, S., Romero-Ortega, M., Suh, S., Hu, Z., ... & Bak, M. 
% % (2015, April). In-vivo tests of a 16-channel implantable wireless neural stimulator. 
% % 2015 7th International IEEE/EMBS Conference on Neural Engineering (NER) (pp. 474-477). IEEE.
e(1).x = .2; e(1).y = 0; e(1).radius =700/1e+7;
e(2).x = .6; e(2).y = 0; e(2).radius =500/1e+7;
e(3).x = 1; e(3).y = 0; e(2).radius =700/1e+7;

e(4).x = 0; e(4).y = .4; e(4).radius =500/1e+7;
e(5).x = .4; e(5).y = .4; e(5).radius =700/1e+7;
e(6).x = .8; e(6).y = .4; e(6).radius =500/1e+7;
e(7).x = 1.2; e(7).y = .4; e(7).radius =700/1e+7;
e(8).x = 1.6; e(8).y = .4; e(8).radius =500/1e+7;

e(9).x = .2; e(9).y = .8; e(9).radius =500/1e+7;
e(10).x = .6; e(10).y = .8; e(10).radius =700/1e+7;
e(11).x = 1; e(11).y = .8; e(11).radius =500/1e+7;
e(12).x = 1.4; e(12).y = .8; e(12).radius =700/1e+7;

e(13).x = .2; e(13).y = .8; e(13).radius =500/1e+7;
e(14).x = .6; e(14).y = .8; e(14).radius =700/1e+7;
e(15).x = 1; e(15).y = .8; e(15).radius =500/1e+7;
e(16).x = 1.4; e(16).y = .8; e(16).radius =700/1e+7;

end


