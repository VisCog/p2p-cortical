% SimulateMapseccList
%
% Simulation of ocular dominance, orientation, & on-off organization
% Written by GMB & IF
% 25/02/2023 moved into clean folder (IF)

clear
clear all


rng(1171964)  % fix the random number generator. This affects the ocular dominance/orientation maps

% define cortex & visual space
c.cortexHeight = [-35, 35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-5, 72];
c.cropLength  = [ 0 70];
v.eccList = [1 2 4  8 16 32];
c.pixpermm =34; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased

% transform to visual space, this should be scaled to which electrodes you
% are simulating
v.visfieldHeight = [-32,32];
v.visfieldWidth= [-32,32];
v.pixperdeg = 34;  %visual field map size and samping

v = p2p_c.define_visualmap(v); % defines the visual map
c = p2p_c.define_cortex(c); % define the properties of the cortical map
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface


figNum = 1; figure(figNum); clf;
p2p_c.plotcortgrid('', c, gray(256), figNum,[ 'title(''cortex'')']); drawnow;
savefig('figures/SimulateMaps_Fig1');
figNum = 10;figure(figNum); clf;
p2p_c.plotretgrid(zeros(size(v.X)), v, gray(256), figNum); drawnow
savefig('figures/SimulateMaps_Fig1_Inset');
% ocular dominance
figNum = 2;figure(figNum); clf;
p2p_c.plotcortgrid(c.ODmap*256, c, gray(256), figNum,[ 'title(''ocular dominance'')']); drawnow;
savefig('figures/SimulateMaps_Fig2');
figNum = 12; figure(figNum); clf;
p2p_c.plotcortgrid(c.ODmap*256, c, gray(256), figNum,[ 'title(''ocular dominance'')']); drawnow; axis square
set(gca, 'YLim', [-2.5 2.5]); set(gca, 'XLim', [5 10]); 
savefig('figures/SimulateMaps_Fig2_Inset');
% orientaton
figNum = 3;
figure(figNum); clf;
p2p_c.plotcortgrid((c.ORmap+pi)*256/(2*pi), c, phasemap(256), figNum,[ 'title(''orientation'')']); drawnow;
figNum = 13; figure(figNum); clf;
p2p_c.plotcortgrid((c.ORmap+pi)*256/(2*pi), c, phasemap(256), figNum,[ 'title(''orientation'')']); drawnow; axis square
savefig('figures/SimulateMaps_Fig3');
set(gca, 'XLim', [5 10])
set(gca, 'YLim', [-2.5 2.5]); 
savefig('figures/SimulateMaps_Fig3_Inset');
% on vs. off cells
figNum = 4; figure(figNum); clf;
p2p_c.plotcortgrid(c.ONOFFmap*256, c, gray(256), figNum,[ 'title(''on vs. off'')']); drawnow;
savefig('figures/SimulateMaps_Fig4');
figNum =14; figure(figNum); clf;
p2p_c.plotcortgrid(c.ONOFFmap*256, c, gray(256), figNum,[ 'title(''on vs. off '')']); drawnow; axis square
set(gca, 'XLim', [5 10])
set(gca, 'YLim', [-2.5 2.5]);
savefig('figures/SimulateMaps_Fig4_Inset');
% simple vs. complex
figNum = 5;figure(figNum); clf;
mnx =max(abs(c.DISTmap(:)));
p2p_c.plotcortgrid((c.DISTmap+mnx)*256/(2*mnx), c, redgreen(256), figNum,[ 'title(''simple vs. complex'')']); drawnow;
savefig('figures/SimulateMaps_Fig5');
figNum =15; figure(figNum); clf;
p2p_c.plotcortgrid((c.DISTmap+mnx)*256/(2*mnx), c, redgreen(256), figNum,[ 'title(''simple vs. complex '')']); drawnow; axis square
set(gca, 'XLim', [5 10])
set(gca, 'YLim', [-2.5 2.5]); 
savefig('figures/SimulateMaps_Fig5_Inset');
% rfsize
figNum = 6;figure(figNum); clf;
mnx =max(abs(c.DISTmap(:)));
p2p_c.plotcortgrid(c.RFsizemap*74, c, cool(256), figNum,[ 'title(''simple vs. complex'')']); drawnow;
savefig('figures/SimulateMaps_Fig6');
figNum = 16;figure(figNum); clf;  imagesc(1:256); colormap(cool(256));
colorbar('eastoutside')
savefig('figures/SimulateMaps_Fig6_Inset');

% save individual receptive fields
sz = 6*v.pixperdeg;
% cortical cells
ct = 1;
n_cells = 12; scFac = 600;
figure(7) ; clf
while ct<=n_cells
    clist = randperm(prod(size(c.X)));
    RF = p2p_c.generate_corticalcell(1000, clist(1), c, v);
    mnRF = mean(RF, 3);
    if ~isnan(sum(mnRF(:)))
        [row, col] = find(abs(mnRF)==max(abs(mnRF(:))));
        lim = sz*1.2;
        if row>lim & col>lim & row<size(c.X, 1)-lim & col<size(c.X, 2)-lim
            subplot(2, n_cells,ct); colormap(gray);
            image(127 +(scFac*RF(row-sz:row+sz, col-sz:col+sz, 1)));axis square; axis off; drawnow;
            t =  title([round(c.v.ECC(clist(1)), 2) round(c.ODmap(clist(1)), 2)]);  set(t, 'FontSize', 6)
            subplot(2,n_cells,ct+n_cells); colormap(gray);
            image(127 + (scFac*RF(row-sz:row+sz, col-sz:col+sz, 2)));       axis square; axis off; drawnow
            ct = ct+1;
        end
    end
end
set(gcf, 'Position', [0.5913    0.1100    0.0501    0.3412]);
savefig('figures/SimulateMaps_Fig7');