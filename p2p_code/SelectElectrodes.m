% Select electrodes
clear all
c.cortexSize = [40,40];c.pixpermm = 6;
v.retinaSize = [30,30]; v.pixperdeg = 5;
c.e_spacing = 0.5;
c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v);


plot(25,0,'r.')
ii=1;
for y = -5:c.e_spacing:10
for x = 15:c.e_spacing:20
    c.e(ii).x = x;
    c.e(ii).y = y;
    ii = ii+1;
end
end

figure(1);
p2p_c.plotcortgrid(zeros(size(c.X)), c, gray(64), 1,['']);
for ii= 1:length(c.e)
    plot(c.e(ii).x, c.e(ii).y, 'r.');
end

v = p2p_c.c2v_define_electrodes(c,v);
figure(2)
p2p_c.plotretgrid(zeros(size(v.X)), v, gray(256), 2,['']);

for ii= 1:length(c.e)
    plot(v.e(ii).x, v.e(ii).y, 'r.');
end