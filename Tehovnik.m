% image2cortex.m
%
% started with p2p_Tehovnik() that if wrote in p2p_main


v.pixperdeg = 12;  %visual field map size and samping
c.pixpermm =12; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased

scFac = 8;
c.cortexHeight = [-20, 20]; % degrees top to bottom, degrees LR,
c.cortexLength = [-60, 0];
v.visfieldHeight = [-15,15];
v.visfieldWidth= [-15,15];

c = p2p_c.define_cortex(c); % define the properties of the cortical map
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface

c.e.radius = 50/1000;
v.e.ang = 0;
v.e.ecc = 3.7;


tp = p2p_c.define_temporalparameters(); % define the temporal model

trl.amp = 50; trl.freq = 1;
trl.pw = 2*10^(-4);   trl.dur= .8;
trl = p2p_c.define_trial(tp,trl);


c = p2p_c.define_electrodes(c, v);

v.target.offset = 0.5; v.target.rad = 0.2;
v = p2p_c.generate_visualtarget(v);
imagesc(v.target)

c = p2p_c.generate_ef(c);
v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode
trl = p2p_c.generate_phosphene(v, tp, trl);

% cortical projection into visual field
p2p_c.plotretgrid(v.e.rfmap(:, :, 1)*1000 + v.target.img*64, v,gray(64), 1, '');
title('rfmap');
p2p_c.plotretgrid(trl.maxphos(:, :, 1)*1000 + v.target.img*64, v, gray(64), 3, '');
title('phosphene');

p2p_c.plotcortgrid(c.e.ef * 64, c, gray(64), 4,'');
title('electric field');
[v, c] = p2p_c.generate_corticalvisualresponse(c, v)
p2p_c.plotcortgrid(4 * c.target.R, c,  gray(64), 5,'');
title('Cortical response to visual target');
p2p_c.plotcortgrid(64 * c.e.ef, c,  gray(64), 6,'');
title('Electrical response');

