% Winawer_ARVO_19.m
% Generates figures for 2019 ARVO talk

figure(10)
clf
subplot(1,2,1)
hold on
subplot(1,2,2)
hold on
clear trl

d = Winawer_getData(1:5);%% load the data
colorList =  [ 0.5 0 0 ;0.5 1 0.5; 1  0.8125 0;0 0.8750  1; 0  0 1.0000];

for ii=1:5
    opts.color = colorList(ii,:);
    for tt = 1:length(d(ii).duration.val)
        site(ii).trl(tt).dur = d(ii).duration.val(tt); % duration in s
        site(ii).trl(tt).pw = d(ii).pulsewidth.val(tt) *10^-6; % pulse width in ms
        site(ii).trl(tt).freq = d(ii).frequency.val(tt); % -1 for a single pulse, NaN if not using a temporal model
        site(ii).trl(tt).amp = d(ii).current.val(tt)*1000; %should be in microAmps
        site(ii).trl(tt).draw.area = d(ii).area.val(tt);
        site(ii).trl(tt).draw.radius = d(ii).polycenterr.val(tt);
        site(ii).trl(tt).draw.brightness = d(ii).brightness.val(tt);
        site(ii).trl(tt).totalCharge = d(ii).totalCharge.val(tt);
    end
    
    c.efthr = 0.01; v.drawthr = 0.05;
    % set up v and c for each electrode
    switch ii
        case 1
            v.e(1).ang = 19.8;      v.e(1).ecc = 26.6;  c.e(1).radius = 1.150 ; % in cm
            c.cortexSize = [80,100]; c.pixpermm = 5;
            v.retinaSize = [70,70]; v.pixperdeg = 5;
        case {2 , 3}
            v.e(2).ang = -166.4;    v.e(2).ecc = 9;     c.e(2).radius = 0.510 ;
            v.e(3).ang = 142.2;     v.e(3).ecc = 5.12;  c.e(3).radius = 1.150 ;
            c.cortexSize = [60,50]; c.pixpermm = 5;
            v.retinaSize = [28,28];v.pixperdeg = 12;
        case {4 , 5} % central electrodes
            v.e(4).ang = 135;       v.e(4).ecc = 1.9;   c.e(4).radius = 1.150 ;
            v.e(5).ang = 146.3;     v.e(5).ecc = 1;     c.e(5).radius = 1.150 ;
            c.cortexSize = [80,50]; c.pixpermm = 5; c.cortexCenter = [10,0];
            v.retinaSize = [8,8];v.pixperdeg = 20;
    end
    tp = p2p_c.define_temporalparameters();
    
    c = p2p_c.define_cortex(c);
    v = p2p_c.define_visualmap(v);
    [c, v] = p2p_c.generate_corticalmap(c, v);
    
    c = p2p_c.define_electrodes(c, v, ii);
    c = p2p_c.generate_ef(c, ii);
    % this is the slow part...
    v = p2p_c.generate_rfmap(c, v, ii);
    p2p_c.plotretgrid(v.e(ii).rfmap(:, :, 1)*255, v, gray(64), 4, ['subplot(2,3,', num2str(ii),  '); title(''rf map'')']);
    p2p_c.plotretgrid(v.e(ii).rfmap_noRF(:, :, 1)*64, v,gray(64), 5, ['subplot(2,3,', num2str(ii), '); title(''rf map no rf'')']);
    
    % Winawer Figure 2: Draw a phosphene for an amplitude of 50Hz, pw 1000 and dur = 1 and amp
    % = 5000 for each of the 5 elecrodes
    
    
    tmp = site(ii).trl(1);
    tmp.freq = 50;
    tmp.pw = 2*10^(-4);
    tmp.dur= 1;
    tmp.e = ii; % which electrode
    tmp = p2p_c.define_trial(tp,tmp);
    tmp = p2p_c.generate_phosphene(v, tp, tmp);
    
    figure(ii)
    clf
    hold on
    p2p_c.plotretgrid(tmp.maxphos(:, :, 1)*1000, v, gray(256), ii,[ 'subplot(1,2,1); title(''phosphene'')';]);
    p2p_c.draw_ellipse(tmp, ii,['subplot(1,2,1)'], 1)
    p2p_c.plotcortgrid(c.e(ii).ef * 100, c, gray(64), ii,['subplot(1,2,2); title(''electric field'')']);
    
    
    for tt = 1:length(site(ii).trl)
        clear tmp;
        tmp = site(ii).trl(tt);
        tmp.e = ii; % which electrode
        tmp = p2p_c.define_trial(tp,tmp);
        tmp = p2p_c.generate_phosphene(v, tp, tmp);
        tmp.draw_area = [site(ii).trl(tt).draw.area];
        
        tmp.draw_radius = [site(ii).trl(tt).draw.radius];
        tmp.draw_diameter = 2 * tmp.draw_radius;
        tmp.sim_radius = mean([tmp.ellipse(1).sigma_x tmp.ellipse(1).sigma_y]);
        tmp.sim_diameter = 2 * tmp.sim_radius;
        
        tmp.draw_brightness = [site(ii).trl(tt).draw.brightness];
        tmp.WCperTrial = site(ii).trl(tt).totalCharge;
        tmp.sim_brightness = max(tmp.maxphos(:));% max of both eyes
        trl(tt) = tmp;
    end
    
    % plot phosphene area as a function of 'Charge deposited per trial'
    %     p2p_c.figure_xyplot(trl, 'WCperTrial','draw_area',[],10, 'subplot(1,2,1)', opts);
    %     p2p_c.figure_xyplot(trl, 'WCperTrial','sim_area',[], 10, 'subplot(1,2,2)', opts);
    drawnow
    
    figure(10)
    subplot(1,2,1)
    x = d(ii).totalCharge.val;
    y = d(ii).poly_area.val;
    
    plot(log(x),log(y),'ko','MarkerFaceColor',colorList(ii,:),'MarkerSize',10);
    set(gca,'XTick',log([1,10,100,1000]));
    set(gca,'XLim',[log(1),log(1000)]);
    set(gca,'YTick',log([.01,.1,1,10,100,1000]));
    set(gca,'YLim',[log(.01),log(1000)]);
    
    
    logx2raw(exp(1),0);
    logy2raw(exp(1),0);
    xlabel('Charge Deposited per Trial (\muC)');
    ylabel('Phosphene area (deg^2)');
    
    x = d(1).totalCharge.val;
    x = [trl(:).WCperTrial];
    y = [trl(:).sim_area];
    
    subplot(1,2,2)
    plot(log(x),log(y),'ko','MarkerFaceColor',colorList(ii,:),'MarkerSize',10);
    set(gca,'XTick',log([1,10,100,1000]));
    set(gca,'XLim',[log(1),log(1000)]);
    set(gca,'YTick',log([.01,.1,1,10,100,1000]));
    set(gca,'YLim',[log(.01),log(1000)]);
    
    logx2raw(exp(1),0);
    logy2raw(exp(1),1);
    xlabel('Charge Deposited per Trial (\muC)');
    ylabel('Phosphene area (deg^2)');
    
    
    
    %    p2p_c.figure_xyplot(trl, 'CperTrial', 'draw_area','sim_area',10, ['subplot(2,2,1); title(''CperT vs. Area'')'], opts);
    %     p2p_c.figure_xyplot(trl, 'CperTrial', 'draw_radius','sim_radius',10, ['subplot(2,2,2); title(''CperT vs. Radius'')'], opts);
    %     p2p_c.figure_xyplot(trl, 'CperTrial', 'draw_brightness','sim_brightness',10, ['subplot(2,2,3); title(''CperT vs. Brightness'')', 0], opts);
    %     p2p_c.figure_xyplot(trl, 'CperPulse', 'draw_brightness', 'sim_brightness', 10, ['subplot(2,2,4); title(''CperP vs. Brightness'')'], opts);
    %
    %     p2p_c.figure_xyplot(trl, 'draw_area','sim_area',[], 11, ['subplot(2,2,1); title(''Draw vs.Sim - area'')'], opts);
    %     p2p_c.figure_xyplot(trl, 'draw_radius','sim_radius',[], 11, ['subplot(2,2,2); title(''Draw vs.Sim -Radius'')'], opts);
    %     p2p_c.figure_xyplot(trl, 'draw_brightness','sim_brightness',[], 11, ['subplot(2,2,3); title(''Draw vs.Sim - Brightness'')', 0], opts);
    %
    drawnow
end

function d = Winawer_getData(sites)

pth = fullfile(ebsRootPath, 'data', 'ebs');
for ii = sites
    fname = sprintf('trial_data_site%d', ii);
    d(ii) = load(fullfile(pth, fname));
end
end