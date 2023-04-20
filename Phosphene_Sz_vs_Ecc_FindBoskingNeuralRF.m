% Phosphene_Sz_vs_Ecc.
% simulates the effect of electrode size and spacing as a function of
% eccentricity
%
% creates a xlsx file of size vs. eccentriticy values as a function of electrode
% sizes called Phosphene_Sz_vs_Ecc.xlsx, may have an appendix, 'keleris' or %'bosking' depending on what model of phosphene size as a function of eccentricity % is being used by the model
%
% written IF & GMB
%
% 25/02/2023 moved into clean folder (IF)
% 06/03/2023 added bosking/keliris functionality (IF)

clear all; close all

eLoc = 3; %exp(linspace(log(.5), log(35), 15)); % the eccentricitity of the electrode
eSize =  exp(linspace(log(.05), log(5), 10)); % esize from teeny tiny to huge

%% generate each cortical surface
% define cortex & retina
c.cortexHeight = [-3,3]; % degrees top to bottom, degrees LR,
c.cortexLength = [-20, 5];
c.pixpermm = 28; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased
c.rfsizemodel = 'invented';
c = p2p_c.define_cortex(c); % define the properties of the cortical map

% transform to visual space
v.visfieldHeight = [-10,10];
v.visfieldWidth= [-5,60];
v.eccList = round(exp(linspace(log(1), log(60), 5)));
v.pixperdeg = 28;  %visual field map size and samping
v = p2p_c.define_visualmap(v); % defines the visual map

eLoc = .5;  %exp(linspace(log(.5), log(35), 15)); % the eccentricitity of the electrode
c.slope =   0.06; % Bosking data is in terms of diameter, so take the values from Figure 4 (slope = 0.2620 and intercept  = 0.1787) and divide by 2
c.intercept =0.08;
c.min = 0;
c.delta = 2; % to do with on off receptive fields

[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
tp = p2p_c.define_temporalparameters(); % define the temporal model

PLOT = 0;
% define pulse train
trl.amp = 300; trl.freq = 1;
trl.pw = 2*10^(-4);   trl.dur= 1;
trl = p2p_c.define_trial(tp,trl);

ct = 1;
% move along the cortex, calulating the size of the percept and it's shift
% in location as one moves foveal to peripheral

    for ecc = 1:length(eLoc)
        c.e.radius = .01; %25;
        v.e.ecc = eLoc(ecc);
        v.e.ang = 0;
        [c, v] = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
        c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface

        v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode
        trl = p2p_c.generate_phosphene(v, tp, trl);
        if PLOT
            figure(1); clf%subplot(length(eSize), length(eLoc), ct);
            p2p_c.plotcortgrid(c.e.ef*2, c, gray(256), 1,['title(''electric field'')']); drawnow;
            figure(2); clf%subplot(length(eSize), length(eLoc), ct);
            p2p_c.plotretgrid(trl.maxphos(:, :, 1)*20, v, gray(256), 2,['';]);
        end
        p.x = v.e.x; p.y = v.e.y;
        p.sigma = .5;
        [p, err] = fit('p2p_c.fit_phosphene', p, {'sigma'}, trl, v);
        ellipse= mean([trl.ellipse(1).sigma_x trl.ellipse(1).sigma_y]);
        disp(['diameter = ', num2str(2*ellipse)]);
        disp([ 'desired diameter = ', num2str(0.1787+ eLoc(ecc)*0.262)])
        sigma = p.sigma;
        esize = eSize*1000;
        eccentricity =eLoc(ecc);
        ct = ct + 1;
    end
    % sim_sigma_radius, the rows represent electrode sizes, the columns
% represent ecccentricity

