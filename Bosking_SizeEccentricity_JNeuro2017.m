% Bosking_Brightness_JNeuro2017
%
% Simulates phosphene brightness as a function of electrode location and pulse
% train
%
% Saturation in Phosphene Size with Increasing Current Levels Delivered to Human Visual Cortex
% William H. Bosking, Ping Sun, Muge Ozker, Xiaomei Pei, Brett L. Foster, Michael S. Beauchamp, Daniel Yoshor
%Journal of Neuroscience 26 July 2017, 37 (30) 7188-7197; DOI: 10.1523/JNEUROSCI.2896-16.2017
%
% Written by GMB & IF
% 25/02/2023 moved into clean folder (IF)

clear

Tloc = readtable('datasets/Bosking2017_data.xlsx'); % note some foveal electrodes are faked since couldn't be identified on the plot
n = size(Tloc, 1);

%% define cortical and visual space
c.cortexHeight = [-35,35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 80];
c.pixpermm = 12;
c = p2p_c.define_cortex(c);

v.visfieldHeight = [-10,10]; v.visfieldWidth= [0,60]; v.pixperdeg = 12;
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
% temporal parameters
tp = p2p_c.define_temporalparameters();

% set up plots
figure(1); hold on
p2p_c.plotretgrid(0, v, gray(64), 1, ['']);

trl.amp = [1000];
trl.dur =  200*10^-3;
trl.pw =  .1 * 10^-3;
trl.freq  = 200 * 10^-3;
trl = p2p_c.define_trial(tp, trl);
trl = p2p_c.convolve_model(tp, trl);
eccList = linspace(.1, 25, 5);
for ii=1:n
    disp(sprintf('Electrode %d of %d',ii,n));

   v.e.ecc = Tloc.ecc(ii);  v.e.ang = 0; v.drawthr = 3;
   % v.e.ecc = eccList(ii); v.e.ang = 0;
    c.e.radius = 0.25;
    c = p2p_c.define_electrodes(c, v);
    c = p2p_c.generate_ef(c);
    v = p2p_c.generate_corticalelectricalresponse(c, v);
    trl = p2p_c.generate_phosphene(v, tp, trl);
    img = mean(trl.maxphos, 3); % average over both eye-dominant cells
    p2p_c.plotretgrid((img./max(img(:)))*256, v, gray(256), 2, ['';]);    drawnow;
    trl.sim_radius= mean([trl.ellipse(1).sigma_x trl.ellipse(1).sigma_y]);
    trl.sim_diameter = 2 * trl.sim_radius; % we usually use radius/sigma, but Bosking uses diameter in this figure so we do too.
    trl.sim_brightness = max(trl.maxphos(:));
    sim_sizes(ii) =  trl.sim_diameter;
    sim_ecc(ii) = v.e.ecc;
end

%%
% Plot phosphene size as a function of phosphene eccentricity (for
% stimulation at 1000 microamps)

figure(1); clf; hold on
plot(sim_ecc, sim_sizes, 'ks', 'MarkerSize', 12,'MarkerFaceColor',[.75,.75,.75],'LineWidth',2); % Figure 4, Bosking 2017
xlabel('Eccentricity (deg)'); ylabel('Simulated Phosphene size (deg)');
set(gca,'FontSize',20);set(gca, 'YLim', [0 6]);
savefig('Bosking_SizeEccentricity1.fig');

figure(2); clf; hold on
plot(Tloc.ecc, Tloc.size, 'ks', 'MarkerSize', 12,'MarkerFaceColor',[.75,.75,.75],'LineWidth',2); % Figure 4, Bosking 2017
xlabel('Eccentricity (deg)'); ylabel('Patient Phosphene size (deg)');
set(gca,'FontSize',20);set(gca, 'YLim', [0 6]);
savefig('Bosking_SizeEccentricity2.fig');
% so the best fitting slope = 0.1787+ ecc*0.262, collected from the graph
% itself. The values here will be slightly different because we couldn't
% grab all the foveal electrodes.

figure(3); clf; hold on
plot(sim_sizes, Tloc.size, 'ks', 'MarkerSize', 12,'MarkerFaceColor',[.75,.75,.75],'LineWidth',2); hold on% Figure 4, Bosking 2017
xlabel('Simulated Phosphene size (deg)');ylabel('Patient Phosphene size (deg)');
set(gca,'FontSize',20);set(gca, 'YLim', [0 6]);set(gca, 'XLim', [0 6]);
plot([0 6], [0 6])

savefig('Bosking_SizeEccentricity3.fig');
