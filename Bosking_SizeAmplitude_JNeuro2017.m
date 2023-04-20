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


n = 5;
Tloc = readtable('datasets/Bosking2017_data.xlsx'); % note some foveal electrodes are faked since couldn't be identified on the plot
eid = randperm(size(Tloc,1)); eid = eid(1:n);
Tloc = Tloc(eid,:);

%% define cortical and visual space
c.cortexHeight = [-15,15]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 0];
c.pixpermm = 12;
c = p2p_c.define_cortex(c);

v.visfieldHeight = [-20,20]; v.visfieldWidth= [0,60]; v.pixperdeg = 12;
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
% temporal parameters
tp = p2p_c.define_temporalparameters();
tp.model = 'compression';
tp.slope = 0.4300;
tp.power = 10;

% set up plots
figure(1); hold on
p2p_c.plotretgrid(0, v, gray(64), 1, ['']);

amp = [   50 100 150 200 250 500,750,1000,1500,2000,3000,4000]/10;
dur =  200*10^-3*ones(size(amp));
pw =  .1 * 10^-3*ones(size(amp));
freq  = 200 *ones(size(amp));

Tsim = table(amp', pw', freq', dur');
Tsim.Properties.VariableNames = {'amp', 'pw','freq','dur'};
all_trl = p2p_c.loop_convolve_model(tp, Tsim);


for ii=1
    disp(sprintf('Electrode %d of %d',ii,n));
    v.e.ecc = Tloc.ecc(ii);
    v.e.ang = 0; % Tloc.ang(ii);
    c.e.radius = 0.25;
    c = p2p_c.define_electrodes(c, v);
    c = p2p_c.generate_ef(c);

    figure(1) % show electrode location in visual space
    plot([v.e.ecc].*cos([v.e.ang]),[v.e.ecc].*sin([v.e.ang]),'ko','MarkerFaceColor','w'); drawnow;
    v = p2p_c.generate_corticalelectricalresponse(c, v);
    for tt=1:length(amp)
        trl = all_trl(tt);
        trl = p2p_c.generate_phosphene(v, tp, trl);
        img = mean(trl.maxphos, 3); % average over both eye-dominant cells
        p2p_c.plotretgrid((img./max(img(:)))*256, v, gray(256), 2, ['';]);
        trl.sim_radius= mean([trl.ellipse(1).sigma_x trl.ellipse(1).sigma_y]);
        trl.sim_diameter = 2 * trl.sim_radius;
        trl.sim_brightness = max(trl.maxphos(:));
        trl.maxresp = max(trl.resp);
        sim_sizes(ii,tt) =  trl.sim_diameter;
        disp(['amp = ', num2str(amp(tt)), ' size ', num2str(sim_sizes(ii,tt))]);
        drawnow;
    end
end

% Plot normalized size as a function of Current (figure 3D)
max_size = max(sim_sizes, [], 2); % take the max over all the amplitudes
y = sim_sizes./repmat(max_size,1,length(amp));

figure(2); clf; hold on
plot([0 amp(1:end)],[0 y(:,1:end)] ,'k-','LineWidth',1);
plot(amp(1:end),y(:,1:end),'s','Color','k','MarkerFaceColor',.75*[1,1,1],'MarkerSize',18,'LineWidth',2);
xlabel('Current (mA)');
ylabel('Normalized phosphene size');
set(gca,'XLim',[0,400]);

set(gca,'YLim',[0 1]);
set(gca,'FontSize',8);



