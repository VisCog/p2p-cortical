function Optimal_Spacing_RF(range)
%
% Creates rfmaps for optimally spaced electrodes.
% The locations of the electrodes were generated using Phosphene_Sz_2D
% and were saved in the dataset Optimal_Spacing.xlsx
%
% Once you've run this you can use these rfs to make movies using
% OptimealSpacing_Movie.m
% 25/02/2023 moved into clean folder (IF)


rng(11171964)  % fix the random number generator. This affects the ocular dominance/orientation maps

T = readtable("datasets/Optimal_Spacing.xlsx");

v.visfieldHeight = [-35, 35]; v.visfieldWidth= [-35,35];
v.pixperdeg = 15;
v = p2p_c.define_visualmap(v);
c.cortexHeight = [-40,40]; % degrees top to bottom, degrees LR,
c.rfsizemodel = 'bosking';
c.pixpermm = 15;
c.e.radius = 0.1;

%% define cortical and visual space
% find and plot the locations of the electrodes
eid = find(T.cx>0);
T = T(eid, :);

for e = range(1):range(2)
    disp(['electrode  ', num2str(e), ' out of ', num2str(range(2))]);
    c.cortexLength = [-5, 80];
    c = p2p_c.define_cortex(c);
    [c, v] = p2p_c.generate_corticalmap(c, v);

    c.e.x = T.cx(e);
    c.e.y = T.cy(e);
    v = p2p_c.c2v_define_electrodes(c, v);
    c = p2p_c.define_electrodes(c, v);
    c = p2p_c.generate_ef(c);
    v = p2p_c.generate_corticalelectricalresponse(c, v); % create receptive field map for each electrode
    % convert to uint8
    rfmap = mean(v.e.rfmap, 3);
     cx = c.e.x;cy= c.e.y; vx = v.e.x; vy =  v.e.y; v_ecc = v.e.ecc; v_ang = v.e.ang;
    save(['datasets/Optimal_Spacing/Optimal_Spacing_RFmaps_bosking_', num2str(e), '.mat'], 'rfmap', 'cx', 'cy', 'vx', 'vy', 'v_ecc', 'v_ang');
    rfmap = fliplr(rfmap); % flip for the other hemisphere to save time
    save(['datasets/Optimal_Spacing/Optimal_Spacing_RFmaps_bosking_', num2str(e+size(T,1)), '.mat'],  'rfmap', 'cx', 'cy', 'vx', 'vy', 'v_ecc', 'v_ang');
end


