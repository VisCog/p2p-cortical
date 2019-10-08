% Beauchamp Bioarchiv 2017.m
clear all; close all
ampList = [800]; % microAmps

rng(1)
c.efthr = 0.05; % what mag of electric field goes through the model, just a  speed thing
tp.scFac = 1;  % scaling of the strength of the pulse train
% v.drawthr = 0.15;
v.drawthr = .015;
[v,xy] = Beauchamp_getData();
figure(1); clf
for ii = 1:length(v.e)
    polarplot(v.e(ii).ang.*pi/180,  v.e(ii).ecc, 'ko'); hold on
    text(v.e(ii).ang.*pi/180, v.e(ii).ecc, num2str(ii));
    ax  = gca;
    ax.RLim = [0 13];
    c.e(ii).radius = 0.25;
end
% order upper panel
stim_order_C = [5 1 3 8 12 20 18 ];

stim_order_N = [3 8 12 20 11 5  9 17];
stim_order_S = [1 3 8  12 11 10 9 13 17 18 19 20];
stim_order_U = [8 12 20 19 14 9 5 ];
filename = 'beauchamp_U_';
stim_order = stim_order_N;
stim_method = 'sequential';
filename = [filename, stim_method];

% Set all electrode radii to .25 mm

c.cortexSize = [80,100]; c.pixpermm = 6; c.efthr = .1;
v.retinaSize = [40,40];v.pixperdeg = 5;  %10

c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
c = p2p_c.define_electrodes(c, v);
tp.model= 'normcdf';
tp = p2p_c.define_temporalparameters(tp);
c = p2p_c.generate_ef(c);

figure(2); clf; figure(1); clf
clist = parula(length(stim_order));
for ii=1:length(stim_order)
    v = p2p_c.generate_rfmap(c, v, stim_order(ii));
    disp(sprintf('Electrode %d of %d',ii,length(stim_order)));
    trl = [];
    if strcmp(stim_method, 'sequential')
        trl.lag = (ii-1)*100*10^-3;% 50ms blank period between each stimulation
    else
        trl.lag = 100*10^-3;
    end
    trl.expname = 'Beauchamp_BioRxiv';
    trl.e = stim_order(ii);
    trl = p2p_c.define_trial(tp,trl);
    figure(1);
    plot(trl.pt, 'Color', clist(ii,:)); hold on
    trl = p2p_c.generate_phosphene(v, tp, trl);
    trl.sim_radius= mean([trl.ellipse(1).sigma_x trl.ellipse(1).sigma_y]);
    p2p_c.plotretgrid(trl.maxphos(:, :, 1)*40, v, gray(256), ii+10,['title(''phosphene'')';]);
    savedtrial(ii) = trl;
    drawnow
end

%% assuming simultaneous stimulation
figure(3); clf
img = savedtrial(1).maxphos(:, :, 1)*40;
for ii=2:length(stim_order)
    img = img+(savedtrial(ii).maxphos(:, :, 1)*40); hold on
    img(99:101, 99:101) = 255;
end
image(v.x, v.y, img); colormap(gray);
%% creating movie
close all
figure(1); clf
vid = VideoWriter([filename, '.avi']);
vid.FrameRate = 30;
open(vid);
itpl= round(linspace(1,length(trl.pt), 30*trl.trialdur));
img = zeros(size(squeeze(savedtrial(ii).maxphos(:, :, 1))));
for k = 1:length(itpl)
    img = zeros(size(squeeze(savedtrial(ii).maxphos(:, :, 1))));
    for ii=1:length(stim_order)
        img = img + squeeze(savedtrial(ii).maxphos(:, :, 1))./max(savedtrial(ii).resp).*savedtrial(ii).resp(itpl(k));
    end
    img(99:101, 99:101) = 255;
    image(v.x, v.y, img*30); axis tight manual; axis off; colormap(gray);
    frame = getframe(gca);
    writeVideo(vid,frame);
end
close(vid)
save(filename, 'savedtrial', 'stim_order');
%%

function [v, xy] = Beauchamp_getData()
xy = [5	8.678; 4.0805	7.842; 3.046	7.519; 2.3565	7.495; 6.4945	7.925; 5.2875	7.366; ...
    4.4255	7.049; 3.2185	6.6625; 7.931	6.653; 6.954	6.1595; 5.9195	5.951; 4.8275	5.741; ...
    8.9655	5.539; 7.701	5.3805; 6.8965	4.9505; 5.747	4.796; 10.23	3.916; 9.138	3.9935; ...
    7.701	3.7715; 6.5515	3.7315; 10.23	3.2265; 9.483	3.028; 8.3335	2.586; 6.954	2.309];

[theta,rho] = cart2pol(xy(:, 1),xy(:, 2));
theta = theta * 180/pi;

for e = 1:size(xy, 1)
    v.e(e).ang = theta(e);
    v.e(e).ecc = rho(e);
end
end


