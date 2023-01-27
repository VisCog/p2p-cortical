% Bosking_JNeuro_17.m

n_e_sample = 10; % how many electrodes to simulate, if NaN simulates all of them
%%  begin by defining the location and size of electrodes
v = Bosking_getData(42);

% add some electrodes near the fovea, since we couldn't sample those from
% the paper
morevals = linspace(.01, 2, 8);
for i=43:50
    v.e(i).ecc = morevals(i-42);
    v.e(i).ang = (rand(1)*2*pi)-pi;
end

% Set all electrode radii to .25 mm
for i=1:length(v.e)
    c.e(i).radius = 0.25;
end

ss = randperm(length(c.e), n_e_sample); % collect a random subsample of the electrodes
c.e = c.e(ss);
v.e = v.e(ss);

%% define cortical and visual space
c.cortexSize = [80,100]; % degrees top to bottom, degrees left to right, divide by 2 to get the actual mm that are useful

c.pixpermm = 6; 
v.retinaSize = [40,40];v.pixperdeg = 10;

c.cortexSize = [80,100]; % degrees top to bottom, degrees LR, divide by 2 to get the actual mm that are useful
c.cortexCenter = [30,0];
c.pixpermm = 6; % default 8, resolution of electric field sampling, for very small electrodes may need to be decreased
fov = 30;
v.retinaSize = [fov, fov]*2.5; v.pixperdeg = 10;  %visual field map size and samping, default is 7

c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
c = p2p_c.define_electrodes(c, v);
tp = p2p_c.define_temporalparameters();
c = p2p_c.generate_ef(c);

% show phosphene locations
figure(1)
clf
p2p_c.plotretgrid(0, v, gray(64), 1, ['']);
hold on
plot([v.e.ecc].*cos([v.e.ang]),[v.e.ecc].*sin([v.e.ang]),'ko','MarkerFaceColor','w');

sim_sizes = [];

ampList = [50 100 150 200 250 500,750,1000,1500,2000,3000,4000];
v = p2p_c.generate_rfmap(c, v);
for ii=1:length(v.e)
    disp(sprintf('Electrode %d of %d',ii,length(v.e)));
    for tt=1:length(ampList)
        trl = [];
        trl.expname = 'Bosking';
        trl.e = ii;
        trl.amp = ampList(tt);
        trl.freq = NaN;
        trl = p2p_c.define_trial(tp,trl);
        trl = p2p_c.generate_phosphene(v, tp, trl);
        img = trl.maxphos(:, :, 1);
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


%%
% Plot normalized size as a function of Current (figure 3D)

max_size = max(sim_sizes')';
y = sim_sizes./repmat(max_size,1,length(ampList));

figure(2)
clf
hold on
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

figure(3)
clf
idx = find(ampList == 1000);

hold on
for ee=1:length(v.e)  
    plot(v.e(ee).ecc, sim_sizes(ee,idx), 'ks', 'MarkerSize', 18,'MarkerFaceColor',[.75,.75,.75],'LineWidth',2); % Figure 4, Bosking 2017
if ee>42    plot(v.e(ee).ecc, sim_sizes(ee,idx), 'ks', 'MarkerSize', 18,'MarkerFaceColor',[1,.5,.5],'LineWidth',2); % Figure 4, Bosking 2017
end
end
xlabel('Eccentricity (deg)'); ylabel('Phosphene size (deg)');
set(gca,'FontSize',20);set(gca, 'YLim', [0 6])


%%
% plot example pulse train and response
dt = trl.t(2)-trl.t(1);
figure(4)
clf
plot(dt:dt:(dt*length(trl.pt)),trl.pt,'k-');
set(gca,'YTick',[]);
widen
%heighten(2)
xlabel('Time (s)');
set(gca,'FontSize',16); ylabel('Amplitude');

figure(5)
clf
plot(dt:dt:(dt*length(trl.resp)),trl.resp,'k-');
set(gca,'YTick',[]); widen; 
%heighten
xlabel('Time (s)'); set(gca,'FontSize',16); ylabel('Response');

%%

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


