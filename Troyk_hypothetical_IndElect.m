% Troy_hypothetical.m

clear all; close all; clear mex

rng(1171964)  % fix the random number generator. This affects the ocular dominance/orientation maps
c.efthr = 0.85; % what magnitude of electric field goes through the model, just a speed thing so one doesn't bother processing non-active regions of cortex
v.drawthr = 1; % what phosphene strength is visible, threshold of visibility/oval drawing

% define cortex & retina
c.cortexHeight = [-10,10]; % degrees top to bottom, degrees LR,
c.cortexLength = [-3, 20];
c.pixpermm = 20; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased

c = p2p_c.define_cortex(c); % define the properties of the cortical map


tp = p2p_c.define_temporalparameters(); % define the temporal model
% transform to visual space
v.visfieldHeight = [-5,5];
v.visfieldWidth= [-5,5];
v.pixperdeg = 20;  %visual field map size and samping
v = p2p_c.define_visualmap(v); % defines the visual map
% define pulse train
trl.amp = 50; trl.freq = 50;
trl.pw = 2*10^(-4);   trl.dur= 1;
trl = p2p_c.define_trial(tp,trl);
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface


wfma= create_WFMA(10, 0);
ct = 1;
for e=1:2:length(wfma)
    c.e.radius = wfma(e).radius; c.e.x = wfma(e).x; c.e.y = wfma(e).y;
    v = p2p_c.c2v_define_electrodes(c,v); % convert electrode locations from cortex to visual space
    c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space

    % set up the electrode locations in terms of their positions in the teeny array
    c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface

    figure(1); subplot(4, 4, e);
    p2p_c.plotcortgrid(c.e.ef*2, c, gray(256), 1,['title(''electric field'')']); drawnow;
    v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode

    % generate percepts for a pulse train
    tmp_trl = p2p_c.generate_phosphene(v, tp, trl);
    figure(2)
   subplot(2, 8, ct)
   scFac = 8;
   img = tmp_trl.maxphos(:, :, 1); 
   img = img*scFac +75;
    p2p_c.plotretgrid(img,  v, gray(256), 2,['';]); 
    a = gca; set(a, 'FontSize', 6);
        t = title(['E', num2str(e), 'LE' ]);    set(t, 'FontSize', 6);
     subplot(2,8, ct+8)
    p2p_c.plotretgrid(tmp_trl.maxphos(:, :, 2)*20, v, gray(256), 2,['';]);
    t= title(['E', num2str(e), 'RE' ]);
    set(t, 'FontSize', 6);    a = gca; set(a, 'FontSize', 6);
 ct = ct+1;
end



function e = create_WFMA(xoffset, yoffset)

%Troyk, P., Bredeson, S., Cogan, S., Romero-Ortega, M., Suh, S., Hu, Z., ... & Bak, M.
% % (2015, April). In-vivo tests of a 16-channel implantable wireless neural stimulator.
% % 2015 7th International IEEE/EMBS Conference on Neural Engineering (NER) (pp. 474-477). IEEE.
r1 = sqrt(700/pi)/1e+7; r2 = sqrt(500/pi)/1e+7;
e(1).x = .2 + xoffset; e(1).y = 0 + yoffset;  e(1).radius = r1;
e(2).x = .6+ xoffset; e(2).y = 0 + yoffset;   e(2).radius =r2;
e(3).x = 1+ xoffset; e(3).y = 0 + yoffset;    e(3).radius =r1;

e(4).x = 0+ xoffset; e(4).y = .4 + yoffset; e(4).radius =r2;
e(5).x = .4+ xoffset; e(5).y = .4 + yoffset; e(5).radius =r1;
e(6).x = .8+ xoffset; e(6).y = .4 + yoffset; e(6).radius =r2;
e(7).x = 1.2+ xoffset; e(7).y = .4 + yoffset; e(7).radius =r1;
e(8).x = 1.6+ xoffset; e(8).y = .4 + yoffset; e(8).radius =r2;

e(9).x = .2+ xoffset; e(9).y = .8 + yoffset;   e(9).radius =r2;
e(10).x = .6+ xoffset; e(10).y = .8 + yoffset; e(10).radius =r1;
e(11).x = 1+ xoffset; e(11).y = .8 + yoffset; e(11).radius =r2;
e(12).x = 1.4+ xoffset; e(12).y = .8 + yoffset; e(12).radius =r1;

e(13).x = .2+ xoffset; e(13).y = .12 + yoffset; e(13).radius =r2;
e(14).x = .6+ xoffset; e(14).y = .12 + yoffset; e(14).radius =r1;
e(15).x = 1+ xoffset; e(15).y = .12 + yoffset; e(15).radius =r2;
e(16).x = 1.4+ xoffset; e(16).y = .12 + yoffset; e(16).radius =r1;

end



