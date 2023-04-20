% Winawer_Size_Neuron2016
%
% Simulates phosphene size as a function of electrode location and pulse
% train
%
% Winawer J, Parvizi J. Linking Electrical Stimulation of Human Primary Visual Cortex,
% %Size of Affected Cortical Area, Neuronal Responses, and Subjective Experience.
% Neuron. 2016 Dec 21;92(6):1213-1219. doi: 10.1016/j.neuron.2016.11.008. Epub
%  2016 Dec 8. PMID: 27939584; PMCID: PMC5182175.
%
% Written by GMB & IF
% 25/02/2023 moved into clean folder (IF)

clear
figure(1); clf
figure(2); clf
rng(1171964)  % fix the random number generator. This affects the ocular dominance/orientation maps

T = readtable('datasets/Winawer2016_data.xlsx');
colorList  = [ 0.5    0    0;   0.5 1 0.5;   1  0.8125    0; 0   0.8750   1;     0     0    1 ]; % roughly match colors to Winawer paper

% use same cortex size and center across electrodes
tp = p2p_c.define_temporalparameters();
ct = 1;
for ii=1:5
    v.eccList = [1 2 3 5 8 13 21 34 80];
    v. pixperdeg = 12;
    c.pixpermm = 12;
    v.drawthr = 1;

    % set up v and c for each electrode
    switch ii
        case 1
            v.e.ang = 19.8;      v.e.ecc = 26.6;  c.e.radius = 1.150; % in cm
            v.visfieldHeight = [-0,30]; v.visfieldWidth= [0,40]; % slightly wider units than Winaweer
            c.cortexLength = [-80, -40];
            c.cortexHeight = [-30,30]; % degrees top to bottom, degrees LR
        case 2
            v.e.ang = -166.4;    v.e.ecc = 9;     c.e.radius = 0.510;
            v.visfieldHeight = [-10,10]; v.visfieldWidth= [-20,0];
            c.cortexHeight = [-30,30]; % degrees top to bottom, degrees LR
            c.cortexLength = [20, 60];
        case 3
            v.e.ang = 142.2;     v.e.ecc = 5.12;  c.e.radius = 1.150;
            v.visfieldHeight = [-5,13]; v.visfieldWidth= [-12,5];
            c.cortexHeight = [-30,30]; % degrees top to bottom, degrees LR
            c.cortexLength = [15, 60];
        case 4 % central electrodes
            v.e.ang = 135;       v.e.ecc = 1.9;   c.e.radius = 1.150;
            v.visfieldHeight = [-8,10]; v.visfieldWidth= [-7,5];% same units as Winawer
            c.cortexHeight = [-30,30]; % degrees top to bottom, degrees LR
            c.cortexLength = [5, 40];
        case 5
            v.e.ang = 146.3;     v.e.ecc = 1;     c.e.radius = 1.150;
            v.visfieldHeight = [-4,4]; v.visfieldWidth= [-6,5];
             c.cortexHeight = [-30,30]; % degrees top to bottom, degrees LR
            c.cortexLength = [0, 30];
    end
    v = p2p_c.define_visualmap(v);
    c = p2p_c.define_cortex(c);
    c = p2p_c.define_electrodes(c, v);
    [c, v] = p2p_c.generate_corticalmap(c, v);
    c = p2p_c.generate_ef(c, ii);
    % this is the slow part...
    v = p2p_c.generate_corticalelectricalresponse(c, v);

    clear trl
    eid = intersect(find(T.electrode==ii), find(T.amp > 0) );
    Texp = T(eid,:);
    if ~isempty(Texp)
        trl= p2p_c.loop_convolve_model(tp, Texp); % do the time
        clear sim*
        for t = 1:length(trl)
            tmp_trl = trl(t);
            tmp_trl= p2p_c.generate_phosphene(v, tp, tmp_trl);
            p2p_c.plotretgrid(tmp_trl.maxphos(:, :, 1)*35, v, gray(256), 3,['subplot(2,1,1); title(''phosphene'')';]);
            p2p_c.draw_ellipse(tmp_trl, 3,['subplot(2,1,1)'], '',colorList(ii,:))
            p2p_c.plotcortgrid(c.e.ef * 30, c, gray(256), 3,['subplot(2,1,2); title(''electric field'')']);
            maxphos(ct) = max(tmp_trl.maxphos(:)); ct = ct+1;
            drawnow

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
        figure(1); hold on
        goodVals = ~isnan(Texp.polyarea) & Texp.polyarea>0 ;
        scatter(log10(Texp.sim_area(goodVals)), log10(Texp.polyarea(goodVals)), ...
            'MarkerFaceColor',colorList(ii,:),...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5, 'SizeData', 50); hold on
        err(ii) = sum((log10(Texp.sim_area(goodVals))-log10(Texp.polyarea(goodVals))).^2);

        % plot areas as a function of charge
        figure(2); hold on
        for sp = 1:2
            subplot(2,1,sp)
            if sp==1
                h =   scatter(log10(Texp.totalcharge(goodVals)), log10(Texp.polyarea(goodVals)), ...
                    'ko','MarkerFaceColor',colorList(ii,:),...
                    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5, 'SizeData', 50); hold on
                ylabel('Drawn Phosphene area (deg^2)');
            else
                h = scatter(log10(Texp.totalcharge(goodVals)) , log10(Texp.sim_area(goodVals)), ...
                    'ko','MarkerFaceColor',colorList(ii,:),...
                    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5, 'SizeData', 50); hold on
                ylabel('Simulated  Phosphene area (deg^2)');
            end
        end
    end
    clear v
    clear c
end
disp(['final err is ', num2str(sum(err))]);
figure(1);
plot(log10([.001 1 100 10000]),log10([.001 1 100 10000]), 'k--');
set(gca,'XLim', [log10(.005) log10(500)]);
set(gca,'YLim', [log10(.005) log10(500)])
set(gca,'XTick', log10([.01 1 100]));
set(gca,'YTick', log10([.01 1 100]));
set(gca,'XTickLabel',[.01 1 100]);
set(gca,'YTickLabel',[.01 1 100]);
xlabel('Simulated Phosphene area (deg^2)');
ylabel('Drawn Phosphene area (deg^2)');
logx2raw(10);  logy2raw(10);

figure(2)
for i = 1:2
    subplot(2, 1, i)
    set(gca,'XTick', log10([5 10 50 100 500]));
    set(gca,'YTick', log10([.01 1 100 ]));
    set(gca,'XTickLabel',[ 5 10 50 100 500 ]);
    set(gca,'YTickLabel',[.01  1 100 ]);
    set(gca,'XLim', [log10(3) log10(700)]);
    set(gca,'YLim', [log10(0.0008) log10(2000)]);
    xlabel('Charge Deposited per Trial (\muC)');
    logx2raw(10);  logy2raw(10);
end
