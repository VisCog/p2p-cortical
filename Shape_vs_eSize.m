% Simulate_WFMA
%
% Simulation of indivual phosphenes using Wireless Floating Microelectrode
% Arrays (WFMAs)
%  Based on: Michael P Barry, Roksana Sadeghi, Vernon L Towle, Kelsey Stipp, Patricia Grant, Frank John Lane, Janet P Szlyk, Gislin Dagnelie, Philip R Troyk;
% Contributed Talk Session III: Characteristics of electrically-induced visual percepts in the first human with the Intracortical Visual Prosthesis. Journal of Vision, forthcoming 2023.
%
% Written by GMB & IF
% 25/02/2023 moved into clean folder (IF)

close all
clear all

rng(3); %(1171964)  % fix the random number generator. This affects the ocular dominance/orientation maps
% define pulse train
tp = p2p_c.define_temporalparameters(); % define the temporal model
trl.amp = 120; trl.freq = 200;
trl.pw = 2*10^(-4);   trl.dur= .8;
trl = p2p_c.define_trial(tp,trl);

% e_sizes = exp(linspace(log(0.002), log(1), 8));

e_sizes =  exp(linspace(log(.05), log(5), 10)); % esize from teeny tiny to huge
%e_sizes = exp(linspace(log(.01), log(2),2));
ecc = 2.5;

for m = 1 %:-1:1% sampling
    clear c
    clear v
    v.pixperdeg = 14;  %visual field map size and samping
    c.pixpermm =14; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased
    scFac = 8;
    c.cortexHeight = [-50, 50]; % degrees top to bottom, degrees LR,
    c.cortexLength = [-40, 0];
    v.visfieldHeight = [-2.5,2.5];
    v.visfieldWidth= [0,5];
    switch  m
        case 1
            c.rfmodel = 'ringach';
      %      c.I_k = 1000; % high electric field fall off  -  only stimulating directly under the electrode
        case 2
            c.rfmodel= 'scoreboard';
         %   c.I_k =1000; % high electric field fall off  -  only stimulating directly under the electrode
    end
    for es=1:length(e_sizes)
        rng(pi)
        c.e.radius =  e_sizes(es);
        v.e.ecc  = ecc;      v.e.ang = 0;
        c = p2p_c.define_cortex(c); % define the properties of the cortical map
        v = p2p_c.define_visualmap(v); % defines the visual map
        [c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
        [c, v] = p2p_c.define_electrodes(c,v); % convert electrode locations from cortex to visual space
        c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface
        
      
        % generate percepts for a pulse train
        v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode
        tmp_trl = p2p_c.generate_phosphene(v, tp, trl);
        img = tmp_trl.maxphos(:, :, 1);
        img = (img*scFac) +125;
        p2p_c.plotretgrid(img,  v, gray(256), es,['';]);
        a = gca; set(a, 'FontSize', 6);
        titlestr = ['e_size', num2str(es)];
        t = title([' E', num2str(es), 'LE' ]);    set(t, 'FontSize', 6);
        saveas(gcf,titlestr);

    end
end



