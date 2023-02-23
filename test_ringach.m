% check_ringach
clear all
close all
tp = p2p_c.define_temporalparameters(); % define the temporal model
% define cortex & retina
c.cortexHeight = [-10,10]; % degrees top to bottom, degrees LR,
c.cortexLength = [-3, 20];
c.pixpermm = 20; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased

c = p2p_c.define_cortex(c); % define the properties of the cortical map
% transform to visual space
v.visfieldHeight = [-5,5];
v.visfieldWidth= [-5,5];
v.pixperdeg = 20;  %visual field map size and samping
v = p2p_c.define_visualmap(v); % defines the visual map
% define pulse train
trl.amp = 50; trl.freq = 50;
trl.pw = 2*10^(-4);   trl.dur= 1;
trl = p2p_c.define_trial(tp,trl);
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
v.e.ecc = 4; v.e.ang = 0; c.e.radius = 0.01;
c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space

% set up the electrode locations in terms of their positions in the teeny array
c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface

% this is the code that does it, but i've got the guts below 
 v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode

x0 = 1.6; % x and y position of the cell in space
y0 = -0.1;
od =  0.8528; % ocular dominance
theta = 0;  %orientation
sigma_x = .2; % major axis sd
sigma_y = .8; % minor axis sd
d = .9; % select a random d value
ef = 10000; % strength of the electric field, 100 gets us into useful units
wplus = 2.3;
wminus = 1.7;

 aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
                bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
                cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);

                % calculate the area of of the on and off fields, which is
                % needed to normalize d. on and off fields use the same
                % oriented gaussian, only their central location and their
                % amplitudes differ
                tmp  = exp( -(aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
                A = sqrt(sum(tmp(:)>0.02)/v.pixperdeg.^2);
                % now create the real on and off fields, that are centered
                % on different locations in space and have variable
                % amplitudes
                xe = x0 + (d/2).*cos(theta); ye = y0 - (d/2).*sin(theta);
                xi = x0 - (d/2).*cos(theta); yi = y0 + (d/2).*sin(theta);

                % on subunit for bright dots
                hplus_on =  exp( - (aa*(v.X-xe).^2 + 2*bb*(v.X-xe).*(v.Y-ye) + cc*(v.Y-ye).^2));
                hplus_on = hplus_on./sum(hplus_on(:));

                % off subunit with dark dots
                hplus_off = exp( - (aa*(v.X-xi).^2 + 2*bb*(v.X-xi).*(v.Y-yi) + cc*(v.Y-yi).^2));
                hplus_off = hplus_off./sum(hplus_off(:));

                % add the off component for bright and dark dots, and scale relative
                % amplitudes
                bright = wplus*hplus_on - 0.4.*wminus*hplus_off;
                dark = wminus*hplus_off - 0.4*wplus*hplus_on;
                % left eye
                RF(:, :, 1, 1)  =   od*ef*bright;
                 RF(:, :, 1, 2)  =   od*ef*dark;
                  min(RF(:))
                   max(RF(:))
                 % right eye
                RF(:, :, 2, 1)  =  (1-od)*ef*bright;
                RF(:, :, 2, 2)  =  (1-od)*ef*dark;

subplot(1, 2,1); image(RF(:,:,1, 1)+50); title('on unit');
subplot(1, 2,2); image(RF(:,:,1, 2)+50); title('off unit'); colormap(gray)
