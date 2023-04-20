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
close all

rng(1171964)  % fix the random number generator. This affects the ocular dominance/orientation maps
T = readtable('datasets/Winawer2016_data.xlsx');
colorList  = [ 0.5    0    0;   0.5 1 0.5;   1  0.8125    0; 0   0.8750   1;     0     0    1 ]; % roughly match colors to Winawer paper

% use same cortex size and center across electrodes
tp = p2p_c.define_temporalparameters();

for ii=1:5
    v.eccList = [  1 2 3 5 8 13 21 34];
    v. pixperdeg = 18;
    c.pixpermm = 18;

    % set up v and c for each electrode
    switch ii
        case 1
            v.e.ang = 19.8;      v.e.ecc = 26.6;  c.e.radius = 1.150; % in cm
            v.visfieldHeight = [-15,25]; v.visfieldWidth= [0,35]; % same units as Winawer
            c.cortexLength = [-80, -40];
            c.cortexHeight = [-30,30]; % degrees top to bottom, degrees LR
        case 2
            v.e.ang = -166.4;    v.e.ecc = 9;     c.e.radius = 0.510;
            v.visfieldHeight = [-10,5]; v.visfieldWidth= [-15,0]; % same units as Winawer
            c.cortexHeight = [-30,30]; % degrees top to bottom, degrees LR
            c.cortexLength = [20, 60];
        case 3
            v.e.ang = 142.2;     v.e.ecc = 5.12;  c.e.radius = 1.150;
           v.visfieldHeight = [-5,10]; v.visfieldWidth= [-15,0]; % same units as Winawer
            c.cortexHeight = [-30,30]; % degrees top to bottom, degrees LR
            c.cortexLength = [15, 60];
        case 4 % central electrodes
            v.e.ang = 135;       v.e.ecc = 1.9;   c.e.radius = 1.150;
                v.visfieldHeight = [-3.5,3.5]; v.visfieldWidth= [-5,2]; % same units as Winawer
            c.cortexHeight = [-30,30]; % degrees top to bottom, degrees LR
            c.cortexLength = [5, 40];
        case 5
            v.e.ang = 146.3;     v.e.ecc = 1;     c.e.radius = 1.150;
             v.visfieldHeight = [-3.5,3.5]; v.visfieldWidth= [-5,2]; % same units as Winawer
            c.cortexLength = [0, 20];
    end
    v = p2p_c.define_visualmap(v);
    c = p2p_c.define_cortex(c);
    c = p2p_c.define_electrodes(c, v);

    [c, v] = p2p_c.generate_corticalmap(c, v);

    c = p2p_c.generate_ef(c, ii);
    % this is the slow part...
    v = p2p_c.generate_corticalelectricalresponse(c, v);

    % draw the w
   figure(ii); clf ; hold on
    clear trl
    trl.amp = 1000; trl.freq = 10;
    trl.pw = 2*10^(-4);   trl.dur= 1;
    trl = p2p_c.define_trial(tp,trl);
    trl = p2p_c.generate_phosphene(v, tp, trl);

    p2p_c.plotretgrid(trl.maxphos(:, :, 1)*25, v, gray(256), ii,['subplot(2,1,1); title(''phosphene'')';]);
  %  p2p_c.draw_ellipse(trl, ii,['subplot(2,1,1)'], 1,colorList(ii,:))
    p2p_c.plotcortgrid(c.e.ef * 30, c, gray(256), ii,['subplot(2,1,2); title(''electric field'')']);
    drawnow

    clear v
    clear c
end

