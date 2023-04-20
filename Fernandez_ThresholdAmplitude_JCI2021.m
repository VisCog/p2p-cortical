% Bosking_SizeAmplitude_JNeuro2017
%
% Simulates phosphene size as a function of electrode amplitude
%
% Saturation in Phosphene Size with Increasing Current Levels Delivered to Human Visual Cortex
% William H. Bosking, Ping Sun, Muge Ozker, Xiaomei Pei, Brett L. Foster, Michael S. Beauchamp, Daniel Yoshor
%Journal of Neuroscience 26 July 2017, 37 (30) 7188-7197; DOI: 10.1523/JNEUROSCI.2896-16.2017
%
% Written by GMB & IF
% 25/02/2023 moved into clean folder (IF)
% 03/03/2023 moved electrodes into table format and split into separate size vs.
% amplitude and size vs. eccentricity scripts (IF)


%% define cortical and visual space
c.cortexHeight = [-15,1]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 80];
c.pixpermm = 12;
c = p2p_c.define_cortex(c);

v.visfieldHeight = [-1,20]; v.visfieldWidth= [-1,60]; v.pixperdeg = 12;
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
% temporal parameters
tp = p2p_c.define_temporalparameters();
% set up plots
figure(1); hold on
p2p_c.plotretgrid(0, v, gray(64), 1, ['']);

amp = 0:10:100;
dur =  166*10^-3*ones(size(amp));
pw =  170 * 10^-6*ones(size(amp));
freq  = 300 *ones(size(amp));

Tsim = table(amp', pw', freq', dur');
Tsim.Properties.VariableNames = {'amp', 'pw','freq','dur'};
all_trl = p2p_c.loop_convolve_model(tp, Tsim);

v.e.ecc =15;
v.e.ang = 0; % Tloc.ang(ii);
c.e.radius = 0.25;
c = p2p_c.define_electrodes(c, v);
c = p2p_c.generate_ef(c);
all_trl = p2p_c.loop_convolve_model(tp, Tsim);
    v = p2p_c.generate_corticalelectricalresponse(c, v);
for tt=1:length(amp)
    disp(tt)
    trl = all_trl(tt);
    trl = p2p_c.generate_phosphene(v, tp, trl);
    img = mean(trl.maxphos, 3); % average over both eye-dominant cells
    trl.maxresp = max(trl.resp);
    sim_amp(tt) =  trl.maxresp;
    disp(['amp = ', num2str(amp(tt)), ' size ', num2str(sim_amp(tt))]);
    drawnow;
end

% Plot normalized brightness  as a function of Current (figure 3D)

figure(2); clf; hold on
plot([0 amp(1:end)],[0 sim_amp(:,1:end)],'k-','LineWidth',1);

xlabel('Current (mA)');
ylabel('Normalized brightness');

set(gca,'YLim',[0 1]);
set(gca,'FontSize',8);



