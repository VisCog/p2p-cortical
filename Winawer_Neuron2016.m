% Winawer Neuron 2016.m
clear all; close all

% addpath(genpath('C:\Users\Ione Fine\Documents\code\UWToolbox\UWToolbox\plotting\'));
c.efthr = 0.005;% what magnitude of electric field goes through the model, just a speed thing so one doesn't bother processing non-active regions of cortex
v.drawthr =1;

T = readtable('datasets/Winawer2016_data.xlsx');
colorList  = [ 0.5    0    0;   0.5 1 0.5;   1  0.8125    0; 0   0.8750   1;     0     0    1 ]; % roughly match colors to Winawer paper

% use same cortex size and center across electrodes
c.cortexSize = [70,110];
c.cortexCenter = [30,0];
c = p2p_c.define_cortex(c);
tp = p2p_c.define_temporalparameters();
v.eccList = [1 2 3 5 8 13 21 34];

for ii=1:5
    % set up v and c for each electrode
    switch ii
        case 1
            v.e.ang = 19.8;      v.e.ecc = 26.6;  c.e.radius = 1.150; % in cm
            v.retinaSize = [70,70]; v.pixperdeg = 12; v.retinaCenter = [0, 0];
        case 2
            v.e.ang = -166.4;    v.e.ecc = 9;     c.e.radius = 0.510;
            v.retinaSize = [50,50];  v.pixperdeg = 12; v.retinaCenter =[0, 0];
        case 3
            v.e.ang = 142.2;     v.e.ecc = 5.12;  c.e.radius = 1.150;
            v.retinaSize = [20,20]; v.pixperdeg = 12; v.retinaCenter = [0, 0];
        case 4 % central electrodes
            v.e.ang = 135;       v.e.ecc = 1.9;   c.e.radius = 1.150;
            v.retinaSize = [10,10]; v.pixperdeg = 12; v.retinaCenter = [0, 0];
        case 5
            v.e.ang = 146.3;     v.e.ecc = 1;     c.e.radius = 1.150;
            v.retinaSize = [10, 10]; v.pixperdeg = 12; v.retinaCenter = [0, 0];
    end

    v = p2p_c.define_visualmap(v);
    [c, v] = p2p_c.generate_corticalmap(c, v);
    c = p2p_c.define_electrodes(c, v);
    c = p2p_c.generate_ef(c, ii);
    % this is the slow part...
    v = p2p_c.generate_rfmap(c, v);

    figure(ii); clf ; hold on
    clear trl
    trl.amp = 5000; trl.freq = 50;
    trl.pw = 2*10^(-4);   trl.dur= 1;
    trl = p2p_c.define_trial(tp,trl);
    trl = p2p_c.generate_phosphene(v, tp, trl);

    p2p_c.plotretgrid(trl.maxphos(:, :, 1)*35, v, gray(256), ii,['subplot(2,1,1); title(''phosphene'')';]);
    p2p_c.draw_ellipse(trl, ii,['subplot(2,1,1)'], 1,colorList(ii,:))
    p2p_c.plotcortgrid(c.e.ef * 30, c, gray(256), ii,['subplot(2,1,2); title(''electric field'')']);
    if strcmp(c.e.hemi,'lh')
        set(gca,'YDir','reverse');
        set(gca,'XDir','reverse');
    end
    drawnow
    clear trl
    % Winawer Figure 2: Draw a phosphene for an amplitude of 50Hz, pw 1000 and dur = 1 and amp
    % = 5000 for each of the 5 elecrodes
    eid = find(T.electrode==ii);
    Texp = T(eid,:);
    trl= p2p_c.loop_convolve_model(tp, Texp); % do the time
    clear sim*
    for t = 1:length(trl)
        tmp_trl = trl(t);
        tmp_trl= p2p_c.generate_phosphene(v, tp, tmp_trl);

        if isfield(tmp_trl, 'ellipse')
            sim_radius(t) = mean([tmp_trl.ellipse(1).sigma_x tmp_trl.ellipse(1).sigma_y]);
            sim_diameter(t) = 2 * sim_radius(t);
            sim_brightness(t) = max(tmp_trl.maxphos(:));% max of both eyes
            sim_area(t) = tmp_trl.sim_area;
        end
    end
    % save the ellipse areas for plotting
    Texp.sim_radius = sim_radius'; Texp.sim_diameter = sim_diameter';
    Texp.sim_area = sim_area';
    Texp.sim_brightness = sim_brightness';

    % plot simulated vs. actual area
    figure(6);  hold on
    goodVals = ~isnan(Texp.polyarea) & Texp.polyarea>0 & Texp.sim_area>0;
    plot(log(Texp.sim_area(goodVals)), log(Texp.polyarea(goodVals)), ...
        'ko','MarkerFaceColor',colorList(ii,:),'MarkerSize',10); hold on

    % plot areas as a function of charge
    figure(7);  hold on
    for sp = 1:2
        subplot(2,1,sp)
        if sp==1
            plot(log10(Texp.totalcharge(goodVals)), log10(Texp.polyarea(goodVals)), ...
                'ko','MarkerFaceColor',colorList(ii,:),'MarkerSize',10); hold on
            ylabel('Drawn Phosphene area (deg^2)');
        else
            plot(log10(Texp.totalcharge(goodVals)), log10(Texp.sim_area(goodVals)), ...
                'ko','MarkerFaceColor',colorList(ii,:),'MarkerSize',10); hold on
            ylabel('Simulated  Phosphene area (deg^2)');
        end   
    end
end
figure(6);
    plot([log10(1) log(10000)], [log(1), log(10000)], 'k');
    set(gca,'XLim', [log10(.0001) log10(1000000)]);  
    set(gca,'YLim', [log10(.0001) log10(1000000)])
    set(gca,'XTick', log10([.001 1 100 10000]));
    set(gca,'YTick', log10([.001 1 100 10000]));
    set(gca,'XTickLabel',[.001 1 100 10000]);
    set(gca,'YTickLabel',[.001 1 100 10000]);
xlabel('Simulated Phosphene area (deg^2)');
ylabel('Drawn Phosphene area (deg^2)');
logx2raw(10);  logy2raw(10);

figure(7)
for i = 1:2
subplot(2, 1, i)
set(gca,'XTick', log10([1 10 100 1000]));
set(gca,'YTick', log10([.001 1 100 ]));
set(gca,'XTickLabel',[ 1 10 100 1000]);
set(gca,'YTickLabel',[.001 1 100 ]);
set(gca,'XLim', [log10(1) log10(1001)]);
set(gca,'YLim', [log10(0.005) log10(2000)]);
xlabel('Charge Deposited per Trial (\muC)');
logx2raw(10);  logy2raw(10);
end



