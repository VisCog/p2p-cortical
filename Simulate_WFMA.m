% Simulate_WFMA
%
% Simulation of indivual phosphenes using Wireless Floating Microelectrode
% Arrays (WFMAs)
%  Based on: Michael P Barry, Roksana Sadeghi, Vernon L Towle, Kelsey Stipp, Patricia Grant, Frank John Lane, Janet P Szlyk, Gislin Dagnelie, Philip R Troyk;
% Contributed Talk Session III: Characteristics of electrically-induced visual percepts in the first human with the Intracortical Visual Prosthesis. Journal of Vision, forthcoming 2023.
%
% Written by GMB & IF
% 25/02/2023 moved into clean folder (IF)

clear
clear all

rng(1171964)  % fix the random number generator. This affects the ocular dominance/orientation maps
c.I_k = 1000; % high electric field fall off -  only stimulating directly under the electrode

% define pulse train
tp = p2p_c.define_temporalparameters(); % define the temporal model
trl.amp = 60; trl.freq = 200;
trl.pw = 2*10^(-4);   trl.dur= .8;
trl = p2p_c.define_trial(tp,trl);

%  based on their abstract, the  sets of electrodes had the following visual field positions:
all_wfma(1).x = -3.5; all_wfma(1).y  = -3.5;
all_wfma(2).x = -4; all_wfma(2).y  = -8;
all_wfma(3).x = -30; all_wfma(3).y  = 0;

% sampling
v.pixperdeg =24;  %visual field map size and samping
c.pixpermm = 24; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased
scFac = 8;

% find out where the wfma are in cortex
for ii = 1:length(all_wfma)
    if ii ==1
        c.cortexHeight = [0, 20]; % degrees top to bottom, degrees LR,
        c.cortexLength = [20, 55];
        v.visfieldHeight = [-10,0];
        v.visfieldWidth= [-10,0];
    elseif ii==2
        c.cortexHeight = [0, 20]; % degrees top to bottom, degrees LR,
        c.cortexLength = [20, 55];
        v.visfieldHeight = [-15,0];
        v.visfieldWidth= [-15,0];
    elseif ii==3
        c.cortexHeight = [-15, 15]; % degrees top to bottom, degrees LR,
        c.cortexLength = [20, 70];
        v.visfieldHeight = [-30,30];
        v.visfieldWidth= [-45,0];
    end
    v.e.x = all_wfma(ii).x;     v.e.y = all_wfma(ii).y;
    v = p2p_c.define_visualmap(v); % defines the visual map
    c = p2p_c.define_cortex(c); % define the properties of the cortical map
    [c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
    c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
    wfma = create_WFMA(c.e.x,  c.e.y);
    wfma = wfma(1:2:end);
    ct = 1;
    for e=1:length(wfma)
        c.e.x = wfma(e).x; c.e.y = wfma(e).y;c.e.radius = .25; c.e.radius =  wfma(e).radius;
        v = p2p_c.c2v_define_electrodes(c,v); % convert electrode locations from cortex to visual space
        v = p2p_c.define_visualmap(v); % defines the visual map
        c = p2p_c.define_cortex(c); % define the properties of the cortical map
        [c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface

        c.e.radius = .75; % make electrode bigger so you can see it
        v = p2p_c.c2v_define_electrodes(c,v); % convert electrode locations from cortex to visual space
        c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
        % set up the electrode locations in terms of their positions in the teeny array
        c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface
        figure(ii); subplot(1, length(wfma), e);
        p2p_c.plotcortgrid(c.e.ef*2, c, gray(256), ii,['title(''electric field'')']); drawnow;
        savefig(['figures/Simulate_WFMA_Cortex_Fig', num2str(ii)]);

        % now do it with the real electrode size
        c.e.radius =  wfma(e).radius;
          v = p2p_c.c2v_define_electrodes(c,v); % convert electrode locations from cortex to visual space
        c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
        % set up the electrode locations in terms of their positions in the teeny array
        c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface

        % generate percepts for a pulse train
        v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode
        tmp_trl = p2p_c.generate_phosphene(v, tp, trl);
        figure(10+ii); subplot(2, length(wfma), ct)
        img = tmp_trl.maxphos(:, :, 2);
        img = (img*scFac) +125;
        p2p_c.plotretgrid(img,  v, gray(256), 10+ii,['';]);
        a = gca; set(a, 'FontSize', 6);
        t = title(['E', num2str(e), 'LE' ]);    set(t, 'FontSize', 6);
        subplot(2, length(wfma), ct+length(wfma))
        img = tmp_trl.maxphos(:, :, 2);
        img = img*scFac +125;
        p2p_c.plotretgrid(img, v, gray(256), 10+ii,['';]);
        t= title(['E', num2str(e), 'RE' ]);
        set(t, 'FontSize', 6);    a = gca; set(a, 'FontSize', 6);
         savefig(['figures/Simulate_WFMA_Visual Field_Fig', num2str(ii)]);
        ct = ct+1;
    end
end

function e = create_WFMA(xoffset, yoffset)

%Troyk, P., Bredeson, S., Cogan, S., Romero-Ortega, M., Suh, S., Hu, Z., ... & Bak, M.
% % (2015, April). In-vivo tests of a 16-channel implantable wireless neural stimulator.
% % 2015 7th International IEEE/EMBS Conference on Neural Engineering (NER) (pp. 474-477). IEEE.
r1 = sqrt(500/pi)/1e+7; r2 = sqrt(1000/pi)/1e+7;
r3 = sqrt(1500/pi)/1e+7; r4 = sqrt(2000/pi)/1e+7;
e(1).x = .2 + xoffset; e(1).y = 0 + yoffset;  e(1).radius = r1;
e(2).x = .6+ xoffset; e(2).y = 0 + yoffset;   e(2).radius =r2;
e(3).x = 1+ xoffset; e(3).y = 0 + yoffset;    e(3).radius =r3;

e(4).x = 0+ xoffset; e(4).y = .4 + yoffset; e(4).radius =r3;
e(5).x = .4+ xoffset; e(5).y = .4 + yoffset; e(5).radius =r4;
e(6).x = .8+ xoffset; e(6).y = .4 + yoffset; e(6).radius =r1;
e(7).x = 1.2+ xoffset; e(7).y = .4 + yoffset; e(7).radius =r2;
e(8).x = 1.6+ xoffset; e(8).y = .4 + yoffset; e(8).radius =r4;

e(9).x = .2+ xoffset; e(9).y = .8 + yoffset;   e(9).radius =r2;
e(10).x = .6+ xoffset; e(10).y = .8 + yoffset; e(10).radius =r3;
e(11).x = 1+ xoffset; e(11).y = .8 + yoffset; e(11).radius =r4;
e(12).x = 1.4+ xoffset; e(12).y = .8 + yoffset; e(12).radius =r1;

e(13).x = .2+ xoffset; e(13).y = .12 + yoffset; e(13).radius =r1;
e(14).x = .6+ xoffset; e(14).y = .12 + yoffset; e(14).radius =r2;
e(15).x = 1+ xoffset; e(15).y = .12 + yoffset; e(15).radius =r3;
e(16).x = 1.4+ xoffset; e(16).y = .12 + yoffset; e(16).radius =r4;

end



