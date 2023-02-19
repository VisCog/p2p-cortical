% Winawer Neuron 2016.m
clear all; close all

% addpath(genpath('C:\Users\Ione Fine\Documents\code\UWToolbox\UWToolbox\plotting\'));
c.efthr = 0.05;
v.drawthr = 0.25;

T = readtable('datasets/Winawer2016_data.xlsx');
colorList  = [1 0 0; 0 1 0; 1 .7 0; .3 .3  1; 0 0 1 ]; % roughly match colors to Winawer paper

v.e(1).ang = 19.8;      v.e(1).ecc = 26.6;  c.e(1).radius = 1.150; % in cm
v.e(2).ang = -166.4;    v.e(2).ecc = 9;     c.e(2).radius = 0.510;
v.e(3).ang = 142.2;     v.e(3).ecc = 5.12;  c.e(3).radius = 1.150;
v.e(4).ang = 135;       v.e(4).ecc = 1.9;   c.e(4).radius = 1.150;
v.e(5).ang = 146.3;     v.e(5).ecc = 1;     c.e(5).radius = 1.150;

% use same cortex size and center across electrodes
c.cortexSize = [80,100];
c.cortexCenter = [30,0];
c = p2p_c.define_cortex(c);
tp = p2p_c.define_temporalparameters();
%% plot individual phosphenes
for ii=1:5
    opts.color = colorList(ii,:);
    % set up v and c for each electrode
    switch ii
        case 1
            v.retinaSize = [70,70]; v.pixperdeg = 12; v.retinaCenter = [0, 0];
        case 2
            v.retinaSize = [60,60];  v.pixperdeg = 12; v.retinaCenter =[0, 0];
        case 3
            v.retinaSize = [20,20]; v.pixperdeg = 12; v.retinaCenter = [0, 0];
        case 4 % central electrodes
            v.retinaSize = [10,10]; v.pixperdeg = 12; v.retinaCenter = [0, 0];
        case 5
            v.retinaSize = [10, 10]; v.pixperdeg = 12; v.retinaCenter = [0, 0];
    end
    v.eccList = [1 2 3 5 8 13 21 34];
    v = p2p_c.define_visualmap(v);
    [c, v] = p2p_c.generate_corticalmap(c, v);
    c = p2p_c.define_electrodes(c, v);
    c = p2p_c.generate_ef(c, ii);
    % this is the slow part...
    v = p2p_c.generate_rfmap(c, v);

    % Winawer Figure 2: Draw a phosphene for an amplitude of 50Hz, pw 1000 and dur = 1 and amp
    % = 5000 for each of the 5 elecrodes
    eid = find(T.electrode==ii);
    Texp = T(eid,:);
    trl= p2p_c.loop_convolve_model(tp, Texp); % do the time
    for t = 1; %:length(trl)
        phos= p2p_c.generate_phosphene(v, tp, trl(t));
    end
    figure(ii); clf
    hold on
    p2p_c.plotretgrid(phos.maxphos(:, :, 1)*1000, v, gray(256), ii,['subplot(2,1,1); title(''phosphene'')';]);
    p2p_c.draw_ellipse(tmp, ii,['subplot(2,1,1)'], 1,colorList(ii,:))
    p2p_c.plotcortgrid(c.e(ii).ef * 256, c, gray(256), ii,['subplot(2,1,2); title(''electric field'')']);
    if strcmp(c.e(ii).hemi,'lh')
        set(gca,'YDir','reverse');
        set(gca,'XDir','reverse');
    end
end

%% calculate area as a function of pulse parameters
for ii=1:5
    disp(['site = ', num2str(ii)]);
    clear trl;
    v = p2p_c.define_visualmap(v);
    [c, v] = p2p_c.generate_corticalmap(c, v);
    c = p2p_c.define_electrodes(c, v);
    c = p2p_c.generate_ef(c, ii);
    v = p2p_c.generate_rfmap(c, v);
    tmp.e = ii; % which electrode, also called a site
    for tt = 1:length(site(ii).trl)
        disp(['trial = ', num2str(tt), ' out of ', num2str(length(site(ii).trl))])
        clear tmp;
        tmp = site(ii).trl(tt);
        tmp = p2p_c.define_trial(tp,tmp);

        tmp = p2p_c.generate_phosphene(v, tp, tmp);
        tmp.draw_area = [site(ii).trl(tt).draw.area];

        tmp.draw_radius = [site(ii).trl(tt).draw.radius];
        tmp.draw_diameter = 2 * tmp.draw_radius;
        if isfield(tmp, 'ellipse')
            tmp.sim_radius = mean([tmp.ellipse(1).sigma_x tmp.ellipse(1).sigma_y]);
            tmp.sim_diameter = 2 * tmp.sim_radius;
            tmp.draw_brightness = [site(ii).trl(tt).draw.brightness];
            tmp.WCperTrial = site(ii).trl(tt).totalCharge;
            tmp.sim_brightness = max(tmp.maxphos(:));% max of both eyes
            trl(tt) = tmp;
        end
        % save the areas for plotting
        S{ii}.area(tt) = tmp.sim_area;
    end
    drawnow
end

%% plot area as a function of charge (Figure 4A)

figure(1); clf
for ii = 1:5

    goodVals = ~isnan(W{ii}.area) & W{ii}.area>0;
    plot(log10(S{ii}.area(goodVals)), log10(W{ii}.area(goodVals)), ...
        'ko','MarkerFaceColor',colorList(ii,:),'MarkerSize',10); hold on
    plot([log10(1) log10(300)], [log10(1), log10(300)], 'k');
    set(gca,'XLim', [log10(.9) log10(350)])
    set(gca,'YLim', [log10(.9) log10(350)])
    set(gca,'XTick', log10([1 3 10 30 100 300]));
    set(gca,'YTick', log10([1 3 10 30 100 300]));
    set(gca,'XTickLabel',[1 3 10 30 100 300]);
    set(gca,'YTickLabel',[1 3 10 30 100 300]);
    logx2raw(10);
    logy2raw(10);

end

figure(10)
clf
hold on
figure(11)
clf
hold on
for ii=1:length(W)
    figure(10)
    goodVals = ~isnan(W{ii}.area) & W{ii}.area>0;
    p = polyfit(log(W{ii}.charge(goodVals)),log(W{ii}.area(goodVals)),1);
    xx = [log(1),log(1000)];
    yy = polyval(p,xx);
    plot(xx,yy,'-','Color',colorList(ii,:));
    plot(log(W{ii}.charge),log(W{ii}.area),'ko','MarkerFaceColor',colorList(ii,:),'MarkerSize',10);

    figure(11)

    goodVals = ~isnan(W{ii}.area) & S{ii}.area>0;
    p = polyfit(log(W{ii}.charge(goodVals)),log(S{ii}.area(goodVals)),1);
    xx = [log(1),log(1000)];
    yy = polyval(p,xx);
    plot(xx,yy,'-','Color',colorList(ii,:));
    plot(log(W{ii}.charge),log(S{ii}.area),'ko','MarkerFaceColor',colorList(ii,:),'MarkerSize',10);
end


%%
figure(10)
set(gca,'XTick',log([1,10,100,1000]));
set(gca,'XLim',[log(1),log(1000)]);

set(gca,'YTick',log([.01,.1,1,10,100]));
set(gca,'YLim',log([.005,800]));

logx2raw(exp(1));
logy2raw(exp(1));
xlabel('Charge Deposited per Trial (\muC)');
ylabel('Phosphene area (deg^2)');
set(gca,'FontSize',12);

figure(11)
set(gca,'XTick',log([1,10,100,1000]));
set(gca,'XLim',[log(1),log(1000)]);

set(gca,'YTick',log([.01,.1,1,10,100]));
set(gca,'YLim',log([.005,800]));

logx2raw(exp(1));
logy2raw(exp(1));
xlabel('Charge Deposited per Trial (\muC)');
ylabel('Phosphene area (deg^2)');
set(gca,'FontSize',12);
