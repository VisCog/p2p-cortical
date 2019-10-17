% Select electrodes
% shows the mapping between electrodes on the cortical surface and
% electrodes in visual space

%% panel A
figNum = 1;
c.animal = 'macaque'; c.cortexSize = [60,90];c.pixpermm = 6;
v.retinaSize = [60,60]; v.pixperdeg = 5;

c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v);
p2p_c.plotcortgrid(zeros(size(c.X)), c, [], figNum, 'subplot(1, 2, 1)');
p2p_c.plotretgrid(zeros(size(v.X)), v, [], figNum, 'subplot(1, 2, 2)');

figNum = 2;
c.animal = 'human'; c.cortexSize = [60,90];c.pixpermm = 6;
v.retinaSize = [60,60]; v.pixperdeg = 5;
c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v);
p2p_c.plotcortgrid(zeros(size(c.X)), c, [], figNum, 'subplot(1, 2, 1)');
p2p_c.plotretgrid(zeros(size(v.X)), v, [], figNum, 'subplot(1, 2, 2)');

%% panel B
figNum = 3;
c.animal = 'human'; c.cortexSize = [40,60];c.pixpermm = 30;
v.retinaSize = [60,60]; v.pixperdeg = 5;
c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v);
p2p_c.plotcortgrid((c.ORmap+pi)*256/(2*pi), c, hsv(256), 4);
p2p_c.plotcortgrid(c.ODmap* 256, c, gray(256), 5);

figNum = 3;
c.animal = 'human'; c.cortexSize = [4,4];c.pixpermm = 30;
v.retinaSize = [60,60]; v.pixperdeg = 5;
c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v);
p2p_c.plotcortgrid((c.ORmap+pi)*256/(2*pi), c, hsv(256),7);
p2p_c.plotcortgrid(c.ODmap* 256, c, gray(256), 8);

%% panel C
c.efthr = 0.05;

c.animal = 'human'; c.cortexSize = [60,90];c.pixpermm = 6;
v.retinaSize = [60,60]; v.pixperdeg = 5;

c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v);

%xLoc = 200; yLoc =100;
xLoc = 120; yLoc = 350;
x0 = c.v2c.X(xLoc, yLoc); % x center
y0 = c.v2c.Y(xLoc, yLoc); % y center
                        
theta = pi-c.ORmap(xLoc, yLoc);  %orientation
sigma_x = c.RFmap(xLoc, yLoc) * c.ar; % major axis sd
sigma_y = c.RFmap(xLoc, yLoc); % minor axis sd

% things you need to go from sigmas and theta to G
aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);

% oriented 2D Gaussian
G =  exp( - (aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));

subplot(1,2,1)
ctmp = zeros(size(c.ORmap));
ctmp(xLoc-3:xLoc+3, yLoc-3:yLoc+3) = 1;
p2p_c.plotcortgrid(ctmp*256, c, gray(256),7, 'subplot(1,2,1)');
p2p_c.plotretgrid(G*256, v, gray,7, 'subplot(1,2,2)');


