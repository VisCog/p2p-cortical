% Bosking_JNeuro_17.m
c.e.radius = 0.25;
c.efthr = 0.05;
n_e_sample = NaN; % how many electrodes to simulate, if NaN simulates all of them
%%  begin by defining the location and size of electrodes
all_v= Bosking_getData(42);

% add some electrodes near the fovea, since we couldn't sample those from the paper
morevals = linspace(.01, 2, 8);
for i=43:50
    all_v.e(i).ecc = morevals(i-42);
    all_v.e(i).ang = (rand(1)*2*pi)-pi;
end

if ~isnan(n_e_sample)
ss = randperm(length(all_v.e), n_e_sample); % collect a random subsample of the electrodes
else
    ss = 1:length(all_v.e);
    n_e_sample = length(all_v.e);
end

%% define cortical and visual space
c.cortexHeight = [-35,35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-10, 80]; 
c.pixpermm = 6; 
c = p2p_c.define_cortex(c);

v.visfieldHeight = [-60,60]; v.visfieldWidth= [-60,60]; v.pixperdeg = 12;
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
% temporal parameters
tp = p2p_c.define_temporalparameters();

% set up plots
figure(1); hold on
p2p_c.plotretgrid(0, v, gray(64), 1, ['']);

for ii=1:n_e_sample    
    disp(fprintf('Electrode %d of %d',ii,n_e_sample));
v.e = all_v.e(ss(ii));
c.e.radius = 0.25;
c = p2p_c.define_electrodes(c, v);
c = p2p_c.generate_ef(c);

amp = [50 100 150 200 250 500,750,1000,1500,2000,3000,4000];
dur =  200*10^-3*ones(size(ampList));
pw =  .1 * 10^-3*ones(size(ampList));
freq  = 200 * 10^-3*ones(size(ampList));

Tsim = table(amp', pw', freq', dur');
Tsim.Properties.VariableNames = {'amp', 'pw','freq','dur'};
all_trl = p2p_c.loop_convolve_model(tp, Tsim);

figure(1) % show electrode location in visual space
plot([v.e.ecc].*cos([v.e.ang]),[v.e.ecc].*sin([v.e.ang]),'ko','MarkerFaceColor','w');

v = p2p_c.generate_corticalelectricalresponse(c, v);
    for tt=1:length(amp)
        trl = all_trl(tt);
        trl = p2p_c.generate_phosphene(v, tp, trl, all_trl(tt).resp);
        img = mean(trl.maxphos, 3); % average over both eye-dominant cells
        p2p_c.plotretgrid((img./max(img(:)))*256, v, gray(256), 2, ['';]);
        trl.sim_radius= mean([trl.ellipse(1).sigma_x trl.ellipse(1).sigma_y]);
        trl.sim_diameter = 2 * trl.sim_radius;
        trl.sim_brightness = max(trl.maxphos(:));
        trl.maxresp = max(trl.resp);
        sim_sizes(ii,tt) =  trl.sim_diameter;
        disp(['amp = ', num2str(ampList(tt)), ' size ', num2str(sim_sizes(ii,tt))]);
        drawnow;
    end
end

% Plot normalized size as a function of Current (figure 3D)
max_size = max(sim_sizes')';
y = sim_sizes./repmat(max_size,1,length(ampList));

figure(2); clf; hold on
plot(ampList(1:end)/1000,y(:,1:end),'k-','LineWidth',1);
plot(ampList(1:end)/1000,y(:,1:end),'s','Color','k','MarkerFaceColor',.75*[1,1,1],'MarkerSize',18,'LineWidth',2);
xlabel('Current (mA)');
ylabel('Normalized phosphene size');
set(gca,'XLim',[0,4.100]);
set(gca,'XTick',[0:4]);
set(gca,'YLIm',[0,1.1]);
set(gca,'FontSize',24);

%%
% Plot phosphene size as a function of phosphene eccentricity (for
% stimulation at 1000 microamps)

figure(3); clf; hold on
idx = find(amp == 1000);
for ii=1:n_e_sample
    plot(   all_v.e(ii).ecc, sim_sizes(ii,idx), 'ks', 'MarkerSize', 12,'MarkerFaceColor',[.75,.75,.75],'LineWidth',2); % Figure 4, Bosking 2017
if ii>42    
    plot(all_v.e(ii).ecc, sim_sizes(ii,idx), 'ks', 'MarkerSize', 12,'MarkerFaceColor',[1,.5,.5],'LineWidth',2); % Figure 4, Bosking 2017
end
end
xlabel('Eccentricity (deg)'); ylabel('Phosphene size (deg)');
set(gca,'FontSize',20);set(gca, 'YLim', [0 6])


function v = Bosking_getData(varargin)

if nargin <1 % how many electrodes to simulate
    n = 40;
else
    n = varargin{1};
end
% eccentricity and angle values from the Bosking paper
ecc = [2.186915888 2.973407843 2.373639079 2.374361692 3.485162347 3.161720782 3.903266211 2.654157433 2.842470373 ...
    3.071683206 3.535889777 4.046488101 4.091868195 4.414587147 4.833847191 4.974323153 5.392137971 3.956016957 4.467482416 4.931110897 5.58045091 4.098371712 ...
    4.608536468 5.722661143 6.047836979 2.88467097 5.307303208 5.355140187 6.143366413 5.172897196 6.056074766 6.610897004 7.446093073 7.769534637 8.557905386 ...
    3.21230369 8.461219771 9.30089604 9.441372001 10.50852683 8.608921861 9.028326428 8.193852972 6.940697562 14.13532132 13.99513441 14.22694865 14.37537335];
ecc = repmat(ecc, 1, ceil(n/length(ecc)));
ecc= ecc(randperm(length(ecc))); % select a random subset of electrodes

ang = [-2.103467982 -1.914786272 -2.00638409 -1.910671717 -1.987623682 -1.600120824 -1.629608464 -1.536051334 -1.626277634 -1.623534598 -1.424713451 ...
    -1.423537864 -1.437399993 -1.452241779 -1.265911243 -1.282320478 -1.184061 -0.365460602 -0.195931163 -0.31682069 0.664892401 0.447114914 ...
    2.172582699 2.139960161 2.455605264 2.549358325 2.489207458 2.598977892 2.19330242 -0.76158443 -0.78054077 -0.921954087 -1.030548934 ...
    -1.001649087 -1.096185873 -0.902801816 -0.980145642 -0.963344545 -1.029569278 -1.139143781 -1.578029586 -1.094422493 -0.904173334];
ang = repmat(ang, 1, ceil(n/length(ang)));
ang = ang(randperm(length(ang)))*180/pi;

for e = 1:n
    v.e(e).ang = ang(e);
    v.e(e).ecc = ecc(e);
end
end


