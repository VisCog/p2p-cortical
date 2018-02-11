function ebs_MakeFigure3()
% Reproduce Figure 3 for the following paper:
%  
%   Winawer and Parvizi (2016). Linking Electrical Stimulation of Human
%   Primary Visual Cortex, Size of Affected Cortical Area, Neuronal
%   Responses, and Subjective Experience Neuron. 92(6):1213?1219
%   http://dx.doi.org/10.1016/j.neuron.2016.11.008
%

site = 3;

% Plot radius, in degrees
rad = 8;  

% Set up figure
fH = figure; pos = get(fH, 'Position'); pos(3:4) = [400 1000];
set(fH, 'Color', 'w', 'NumberTitle', 'off', 'Position', pos, ...
    'name', 'Winawer and Parvizi Figure 3')

%% Plot ECoG PRF (panel A)

% Solve the ECoG pRF model
params = ebs_solvePRFmodels(site);
p = params.bb.params;
[x, y, s] = ecogFitPRFPix2Deg(site, p(1), p(2), p(3)/sqrt(p(5)));

% Plot settings
colors = getColors(site);

% Coordinates to draw circle for 1- and 2-sigma lines
angles = linspace(0,2*pi, 50);
ecc    = [angles*0+1*s; angles*0+2*s; ]; 

% Polar to cartesian
[x_c, y_c] = pol2cart([angles; angles], ecc);
to_plot.x = x_c+x;
to_plot.y = y_c+y;

% Plot the PRF
subplot(3,1,1); 
plotPolar(rad, fH);
fill(to_plot.x', to_plot.y', colors, 'EdgeColor', 'none');
plot(to_plot.x', to_plot.y', '--', 'LineWidth', 2, 'Color', 'k')
title('ECoG', 'FontSize', 20)

%% Plot Phosphene  (panel B)

% Use phosphene from condition 20 as an example
which_phosphene  = 20;

% Load the phosphenes from a file
pth = fullfile(ebsRootPath, 'data', 'ebs');
load(fullfile(pth, sprintf('phosphenes_site%d', site)));

% Set up plot
figure(fH); 
subplot(3,1,2)
plotPolar(rad, fH);
    
% Draw the phosphene
fill(phosphenes(which_phosphene).x,phosphenes(which_phosphene).y, colors, 'EdgeColor', 'none');
plot(phosphenes(which_phosphene).x,phosphenes(which_phosphene).y, '--', 'LineWidth', 2, 'Color', 'k')
title('EBS', 'FontSize', 20)

%% Plot fMRI pRF  (panel C)

% Load the fMRI pRF coverage plot for site 3
pth = fullfile(ebsRootPath, 'data', 'fmriPRF');
load(fullfile(pth, 'pRF_COV.mat'));

% Set up the plot
figure(fH); 
subplot(3,1,3); cla
plotPolar(rad, fH);

% 1 and 2 sigma thresholds
onesigma = exp(-1/2);
twosigmas = exp(-2/2);

% x-y coordinates of the image
rg = linspace(-rad,rad,size(pRF_COV,1));

% two sigma fill and line
C = contourc(rg, rg, double(pRF_COV), twosigmas*[1 1]);
fill(C(1,2:end), C(2,2:end), colors, 'EdgeColor', 'none');
plot(C(1,2:end), C(2,2:end), '--', 'LineWidth', 2, 'Color', 'k')

% one sigma line
C = contourc(rg, rg, double(pRF_COV), onesigma*[1 1]);
plot(C(1,2:end), C(2,2:end), '--', 'LineWidth', 2, 'Color', 'k')

title('fMRI', 'FontSize', 20)
end


function fH = plotPolar(rad, fH)

if ~exist('fH', 'var'); fH = figure; end

plot_prefs.linewidth = 2;
plot_prefs.gridcolor = .2*[1 1 1];

hold on
% radial grid lines
for ii = 0:3
    [xi, yi] = pol2cart(pi/4*[ii ii+4], rad *[ 1 1]);
    plot(xi, yi, 'color', plot_prefs.gridcolor, ...
        'LineWidth', plot_prefs.linewidth)
    
end

% circular grid lines
ang=0:0.01:2*pi;

r = rad * (2/3).^(0:7);

r = [1 2];
for ii = 1:6; r(end+1) = r(end-1)+r(end); end

for ii = 1:length(r)
    
    xp=r(ii)*cos(ang);
    yp=r(ii)*sin(ang);
    plot(xp, yp,  'color', plot_prefs.gridcolor, ...
        'LineWidth', plot_prefs.linewidth);
    %text(0, r(ii), num2str(round(r(ii))), 'FontSize', 14)
end

axis off

axis(rad*[-1 .25 -.25 1]), axis square

end

