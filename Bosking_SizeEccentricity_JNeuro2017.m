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
debug_flag = 1;

Tloc = readtable('datasets/Bosking2017_data.xlsx'); % note some foveal electrodes are faked since couldn't be identified on the plot
n = size(Tloc, 1);

%% define cortical and visual space
c.cortexHeight = [-35,35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 80];
if debug_flag
    c.pixpermm = 4;
    v.pixperdeg = 4;
else
    c.pixpermm = 12;
    v.pixperdeg = 12;
end

for m = 1:2
    if m==1
        c.rfsizemodel = 'keliris';
    else
        c.rfsizemodel = 'winawer';
    end
    c = p2p_c.define_cortex(c);
    v.visfieldHeight = [-10,10]; v.visfieldWidth= [0,60];
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
    eccList = exp(linspace(log(.1), log(25), 5));
    for ii=1:n
        disp(sprintf('Electrode %d of %d',ii,n));
        v.e.ecc = Tloc.ecc(ii);  v.e.ang = 0; v.drawthr = 1;
        % v.e.ecc = eccList(ii); v.e.ang = 0;
        c.e.radius = 0.25;
        c = p2p_c.define_electrodes(c, v);
        c = p2p_c.generate_ef(c);
        v = p2p_c.generate_corticalelectricalresponse(c, v);
        trl = p2p_c.generate_phosphene(v, tp, trl);
        img = mean(trl.maxphos, 3); % average over both eye-dominant cells
        p2p_c.plotretgrid((img./max(img(:)))*256, v, gray(256), 2, ['';]);    drawnow;
        trl.sim_radius= mean([trl.ellipse(1).sigma_x trl.ellipse(1).sigma_y]);
        trl.sim_brightness = max(trl.maxphos(:));
        sim_sizes(m, ii) =  trl.sim_radius;
        sim_ecc(m, ii) = v.e.ecc;
    end
end
%%
% Plot phosphene size as a function of phosphene eccentricity (for
% stimulation at 1000 microamps)

figure(1); clf; hold on
plot(Tloc.ecc, Tloc.size, 'ks', 'MarkerSize', 12,'MarkerFaceColor',[.75,.75,.75],'LineWidth',1.5); % Figure 4, Bosking 2017
xlabel('Eccentricity (deg)'); ylabel('Patient Phosphene size (deg)');
set(gca,'FontSize',20);set(gca, 'YLim', [0 6]);

% now look at simulations
figure(2); clf; hold on
plot(sim_ecc, sim_sizes(1, :), 'kd', 'MarkerSize', 12,'MarkerFaceColor',[.75,.75,.75],'LineWidth',1.5); % Figure 4, Bosking 2017
plot(sim_ecc, sim_sizes(2,:), 'k^', 'MarkerSize', 12,'MarkerFaceColor',[.75,.75,.75],'LineWidth',1.5); % Figure 4, Bosking 2017
xlabel('Eccentricity (deg)'); ylabel('Simulated Phosphene size (deg)');
set(gca,'FontSize',20);set(gca, 'YLim', [0 6]);

%% best fitting line
pf = polyfit(Tloc.ecc,Tloc.size, 1);
disp('correlation for best fit line')
[r, p] = corr(Tloc.size,  polyval(pf, Tloc.ecc));
disp(['r =  ', num2str(round(r, 3)), ' p = ', num2str(round(p, 4)), ' df = ', num2str(length(Tloc.size)-2)])

for i = 1:2
    figure(i)
    plot(Tloc.ecc, polyval(pf, Tloc.ecc), 'k--');
end
disp('correlation keliris')
[r, p] = corr(Tloc.size,  sim_sizes(1, :)')
disp(['r =  ', num2str(round(r, 3)), ' p = ', num2str(round(p, 4)), ' df = ', num2str(length(Tloc.size)-2)])
disp('correlation winawer')
[r, p] = corr(Tloc.size,  sim_sizes(2, :)')
disp(['r =  ', num2str(round(r, 3)), ' p = ', num2str(round(p, 4)), ' df = ', num2str(length(Tloc.size)-2)])



% 
% % low rent function minimization
% mf= 1:0.001:4.5;
% for i = 1:length(mf)
%     y = [mf(i)*(c.intercept+ Tloc.ecc*c.slope)];
%     fmin(i)= sum((Tloc.size-y).^2);
% end
% idx = find(fmin==min(fmin));
% % scFac = mf(idx);
% kel_y = [scFac*(c.intercept+ Tloc.ecc*c.slope)];
  %  plot( Tloc.ecc, kel_y, 'k-');
% disp('correlation with real data using Keliris linear model of rf vs. ecc')
% [r, p] = corr(Tloc.size,  [scFac*(c.intercept+ Tloc.ecc*c.slope)]);
% disp(['r =  ', num2str(round(r, 3)), ' p = ', num2str(round(p, 4)), ' df = ', num2str(length(Tloc.size)-2)]);
;


% disp('correlation for best fit line with simulated data')
% [r, p] = corr( sim_sizes', polyval(pf, sim_ecc)');
% disp(['r =  ', num2str(round(r, 3)), ' p = ', num2str(round(p, 4)), ' df = ', num2str(length(sim_sizes)-2)]);


