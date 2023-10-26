 function Optimal_Spacing_RFmap(varargin)
%
% Creates rfmaps for optimally spaced electrodes.
% The locations of the electrodes were generated using Phosphene_Sz_2D
% and were saved in the dataset Optimal_Spacing.xlsx
%
% Once you've run this you can use these rfs to make movies using
% OptimealSpacing_Movie.m
% 25/02/2023 moved into clean folder (IF)

rng(11171964)  % fix the random number generator. This affects the ocular dominance/orientation maps

arrayStyle = 'regular';
if strcmp(arrayStyle, 'optimal')
    T = readtable("datasets/Optimal_Spacing.xlsx");
    eid = find(T.vx>0 & T.vy>0);
elseif strcmp(arrayStyle, 'regular')
    [vx, vy] = meshgrid(linspace(0, 44, 33), linspace(0, 44, 33)); % match resolution of optimal array
    vx = vx(:); vy = vy(:); T = table(vx,vy);
    eid = find(T.vx>0 & T.vy>0);
end
T = T(eid, :);
if nargin <1
    range = [1 height(T)];
else
range = varargin{1};
end
v.visfieldHeight = [0, 35]; v.visfieldWidth= [0,35];
v.pixperdeg = 15;
v = p2p_c.define_visualmap(v);
c.cortexHeight = [-40,0]; % degrees top to bottom, degrees LR
c.cortexLength = [-80, 0];
c.rfsizemodel = 'bosking';
c.pixpermm = 15;
c.e.radius = 0.1;

%% define cortical and visual space
% find and plot the locations of the electrodes

for e =range(1):range(2)
    %disp(['electrode  ', num2str(e), ' out of ', num2str(range(2))]);
    c = p2p_c.define_cortex(c);
    [c, v] = p2p_c.generate_corticalmap(c, v);

    v.e.x = T.vx(e);
    v.e.y = T.vy(e);
 %   v = p2p_c.c2v_define_electrodes(c, v);
    c = p2p_c.define_electrodes(c, v);
    c = p2p_c.generate_ef(c);
    v = p2p_c.generate_corticalelectricalresponse(c, v); % create receptive field map for each electrode
    rfmap = mean(v.e.rfmap, 3);
     cx = c.e.x;cy= c.e.y; vx = v.e.x; vy =  v.e.y;
     save(['datasets/Optimal_Spacing/RFmaps_', c.rfsizemodel, '_', arrayStyle, '_', num2str(e), '.mat'], 'rfmap', 'cx', 'cy', 'vx', 'vy');
     rfmap = fliplr(rfmap); % flip for the other hemisphere to save time
     save(['datasets/Optimal_Spacing/RFmaps_', c.rfsizemodel, '_', arrayStyle, '_', num2str(e+size(T,1)), '.mat'],  'rfmap', 'cx', 'cy', 'vx', 'vy');
end


