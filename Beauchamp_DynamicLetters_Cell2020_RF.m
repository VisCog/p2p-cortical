% Beauchamp_DynamicLetters_Cell 2020.m
%
% Replicates drawing task from Beauchamp dynamic stimulation paper
% Beauchamp, M. S., Oswalt, D., Sun, P., Foster, B. L., Magnotti, J. F., Niketeghad, S., ... & Yoshor, D. (2020).
% Dynamic stimulation of visual cortex produces form vision in sighted and blind humans. Cell, 181(4), 774-783.
% written IF & GMB
%
%  creates RF maps and saves them mat file Beauchamp_DynamicLetters_RF_Figure4.mat etc.
%
% 25/02/2023 moved into clean folder (IF)
% 06/03/2023 Split RF generation and movie generation into separate files
% (IF)  and used electrode array instead of patient report to determine the
% location of phosphenes

clear
rng;  % fix the random number generator. This affects the ocular dominance/orientation maps
expList ={'Figure 4-grid'};% {'Figure 4', 'Figure 6', 
% Figure 4 and Figure 6 replicate the corresponding figures in Beauchamp's
% paper. Figure 4 - grid uses the fact that we know the electrode array has
% 2mm electrode separation and simulated the expected location of
% phosphenes, allowing a certain amount of mashing of the cortical surface
Te = readtable("datasets/Beauchamp_2020_data.xlsx", 'Sheet','ElectrodeLocations');
To = readtable("datasets/Beauchamp_2020_data.xlsx", 'Sheet','ElectrodeOrder');

v.visfieldHeight = [-35, 35]; v.visfieldWidth= [0,35]; v.pixperdeg = 18;
v = p2p_c.define_visualmap(v);
c.cortexHeight = [-35,35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 5];
c.pixpermm = 18; 
c.squish = 0.63;
c.a = .1474;
c = p2p_c.define_cortex(c);
[c, v] = p2p_c.generate_corticalmap(c, v);
PLOT = 0;

%% define cortical and visual space
% find and plot the locations of the electrodes
for ex = 1
    disp(['simulating ', expList{ex}]);
    subplot(1, length(expList), ex);
    eid =strcmp(Te.experiment,expList{ex});
    Tloc = Te(eid, :);
    for e = 1:size(Tloc, 1)
        disp(['simulating electrode ', num2str(e)]);
        v.e.x = Tloc.x(e); v.e.y = Tloc.y(e);
        c.e.radius = Tloc.radius(e);
        c = p2p_c.define_electrodes(c, v);
        c = p2p_c.generate_ef(c);
        v = p2p_c.generate_corticalelectricalresponse(c, v); % create receptive field map for each electrode
        saved(e).rfmap = mean(v.e.rfmap, 3); % save the rf map for that electrode, averaged across both eyes
        if PLOT
            p2p_c.plotretgrid(v.e.rfmap(:, :, 1)*255, v, gray(256), 1); drawnow
        end
        rmfield(v,'e');
    end
    save(['datasets/Beauchamp_DynamicLetters_RF', expList{ex}], 'saved');
end