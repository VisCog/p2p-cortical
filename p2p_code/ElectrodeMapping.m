% Select electrodes
% shows the mapping between electrodes on the cortical surface and
% electrodes in visual space

c.cortexSize = [40,40];c.pixpermm = 6;
v.retinaSize = [30,30]; v.pixperdeg = 5;

c.e_spacing = 0.5;
c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v);

plot(25,0,'r.')
ee=1;
for y = -5:c.e_spacing:10
for x = 15:c.e_spacing:20
    c.e(ee).x = x;
    c.e(ee).y = y;
    ee = ee+1;
end
end

figure(1);
p2p_c.plotcortgrid(zeros(size(c.X)), c, [], 1, 'subplot(1, 2, 1)');
for ee= 1:length(c.e)
    plot(c.e(ee).x, c.e(ee).y, 'r.');
end

v = p2p_c.c2v_define_electrodes(c,v);
figure(2)
p2p_c.plotretgrid(zeros(size(v.X)), v, [], 1, 'subplot(1, 2, 2)');

for ee= 1:length(c.e)
    plot(v.e(ee).x, v.e(ee).y, 'r.');
end