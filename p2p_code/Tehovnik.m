% image2cortex.m
%
% started with p2p_Tehovnik() that if wrote in p2p_main


c.a = 0.3; 
c.cortexSize = [20,45]; 
c.pixpermm = 5; 
c.efthr = .15;
c.slope = .051; 
c.min = 0; 
c.intercept = .214;

c = p2p_c.define_cortex(c);
v.e.ang = 0; 
v.e.ecc = 3.7;
v.retinaSize = [5,15];
v.pixperdeg = 10; 

v = p2p_c.define_visualmap(v);

[c, v] = p2p_c.generate_corticalmap(c, v);


c.e.radius = 50/1000; 
c = p2p_c.define_electrodes(c, v);

v.target.offset = 0.5; v.target.rad = 0.2; 
v = p2p_c.generate_visualtarget(v);
c = p2p_c.generate_corticalresponse(c, v);
c = p2p_c.generate_ef(c);
v = p2p_c.generate_rfmap(c, v);
tp = p2p_c.define_temporalparameters();
trl.expName = 'Tehovnik';
trl = p2p_c.define_trial(tp,trl);

trl = p2p_c.generate_phosphene(v, tp, trl);


% cortical projection into visual field
p2p_c.plotretgrid(v.e.rfmap(:, :, 1)*1000 + v.target.img*64, v,gray(64), 1, '');
title('rfmap');
p2p_c.plotretgrid(v.e.rfmap_noRF*64 + v.target.img*64, v, gray(64), 2, '');
title('rfmap no RF');
p2p_c.plotretgrid(trl.maxphos(:, :, 1)*1000 + v.target.img*64, v, gray(64), 3, '');
title('phosphene');

p2p_c.plotcortgrid(c.e.ef * 64, c, gray(64), 4,'');
title('electric field');
p2p_c.plotcortgrid(4 * c.target.R, c,  gray(64), 5,'');
title('Cortical response to visual target');
p2p_c.plotcortgrid(64 * c.e.ef, c,  gray(64), 6,'');
title('Electrical response');

