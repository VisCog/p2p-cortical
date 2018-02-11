function ebs_visualizeElectrodes(subj, meshData)
% Inputs 
% meshData:  'areas' | 'v1' | 'angle' | 'eccen' ;
% 
% Example: ebs_visualizeElectrodes(1, 'areas')
% Example: ebs_visualizeElectrodes(3, 'v1')

% Subject 1 is LH, subjects 2-4 are RH
if subj == 1, hemi = 'l'; else, hemi = 'r'; end

% Load electrode locations
pth = fullfile(ebsRootPath, 'data', 'anatomy', sprintf('s%d', subj));
ele = load(fullfile(pth, 'ele_xyz.mat'));

% Load T1 anatomy
ni  = MRIread(fullfile(pth, 'brain.mgz'));
surf_offsets = [ni.c_r ni.c_a ni.c_s];

% Surface mesh
surf = fullfile(pth, sprintf('%sh.white', hemi));
[vertices, faces] = freesurfer_read_surf(surf);

% Curvature (for coloring)
fname = fullfile(pth, sprintf('%sh.curv', hemi));
curv = read_curv(fname);
mshcolor.curv = zeros(size(curv));
mshcolor.curv(curv>0) = 1;

% Surface maps
%       angle, eccentrity, and area from Benson et al 
retAngle = fullfile(pth, sprintf('%sh.template_angle.mgz', hemi));
retEccen = fullfile(pth, sprintf('%sh.template_eccen.mgz', hemi));
retArea  = fullfile(pth, sprintf('%sh.template_areas.mgz', hemi));
ni       = MRIread(retAngle); mshcolor.angle = ni.vol(:);
ni       = MRIread(retEccen); mshcolor.eccen = ni.vol(:);
ni       = MRIread(retArea);  mshcolor.areas = abs(ni.vol(:));
idx      = mshcolor.areas > 0;
%       V1 probability maps from Hinds et al
retV1    = read_label([], fullfile(pth, sprintf('%sh.v1.prob', hemi)));
mshcolor.v1= zeros(size(mshcolor.areas));
mshcolor.v1(retV1(:,1)+1) = retV1(:,5);


% Color maps
switch meshData
    case 'curv',  cmap = gray(256); clim = [0 255];
    case 'angle', cmap = jet(256);  clim = [0 191];
    case 'eccen', cmap = jet(256);  clim = [0 20];
    case 'areas', cmap = jet(256);  clim = [0 4];
    case 'v1',    cmap = parula(256); clim = [0 1];   
end

lightColor = [0.5000 0.5000 0.4500];
cmap(1,:) =  2 *lightColor; % color for gyri
cmap(2,:) =  1 *lightColor; % color for sulci



%% Plot

fH = figure; pos = get(gcf, 'Position'); pos([3 4]) = [1000 800];
set(fH, 'Color', 'w', 'Position', pos); clf;

% Mesh
c = mshcolor.curv;
c(idx) = round(mshcolor.(meshData)(idx)/clim(2)*255);
mx = vertices(:,1)+surf_offsets(1);
my = vertices(:,2)+surf_offsets(2);
mz = vertices(:,3)+surf_offsets(3);
tH = trimesh(faces, mx, my, mz, c); axis equal; hold on
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',c)
colormap(cmap); set(gca, 'CLim', [0 255]);

% Electrode
plot3(ele.xyz(:,1), ele.xyz(:,2),ele.xyz(:,3), 'o', 'MarkerSize', 36, ...
    'MarkerFaceColor', 'none', 'Color', 'w', 'LineWidth', 2)

plot3(ele.xyz(:,1), ele.xyz(:,2),ele.xyz(:,3), 'o', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'Color', 'w', 'LineWidth', 2)

axis off;

% Viewpoint
if strcmpi(hemi, 'r'), set(gca, 'View', [-70 -10]);
else, set(gca, 'View', [110 -10]); end

% Lighting
if strcmpi(hemi, 'r'), pos = [-1 1 1]; else, pos = [1 1 1]; end
light('Position',100*pos,'Style','local')
lighting gouraud
material dull

end

