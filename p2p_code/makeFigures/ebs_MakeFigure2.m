function ebs_MakeFigure2()
% Reproduce Figure 2 for the following paper:
%  
%   Winawer and Parvizi (2016). Linking Electrical Stimulation of Human
%   Primary Visual Cortex, Size of Affected Cortical Area, Neuronal
%   Responses, and Subjective Experience Neuron. 92(6):1213?1219
%   http://dx.doi.org/10.1016/j.neuron.2016.11.008
%

% Plot options
opts_bb.sl          = false; % stimulus locked pRF? (if not, then broadband)
opts_bb.rad         = 40;    % maximum radius, in degrees
opts_bb.srate       = 0.1;   % grid sampling rate, in degrees
opts_bb.whichsigmas = [1 2]; % how many pRF sigmas to plot
opts_bb.plottype    = 'contours'; %  'contours' or 'shaded';

opts_sl = opts_bb;
opts_sl.sl = true;

opts_phosp = opts_bb;
opts_phosp.plottype    = 'shaded';

% Which sites
sites = 1:5; 

% Get phosphene and pRF images
[pRFs_bb, r_bb]  = getPRFdata(opts_bb, sites);
[pRFs_sl, r_sl]  = getPRFdata(opts_sl, sites);
[phos, phos_all] = getPhosdata(opts_phosp, sites);

% Display R^2 of pRF models for sites 1-5
fprintf('bb r\t\t'), disp(r_bb)
fprintf('sl r\t\t'), disp(r_sl)

%% Figure 2A: All sites on one plot - broadband v phosphene
fH = figure; 
set(fH, 'Color','w', 'name','Winawer and Parvizi, Figure 2A', 'NumberTitle', 'off')
fH = plotPRFs(opts_phosp, sites, phos, fH);
fH = plotPRFs(opts_bb, sites, pRFs_bb, fH);
axis([-20 35 -20 35])

%% Figure 2B: All sites on one plot - stimulus locked v phosphene
fH = figure; 
set(fH, 'Color','w', 'name','Winawer and Parvizi, Figure 2B', 'NumberTitle', 'off')
fH = plotPRFs(opts_phosp, sites, phos, fH);
fH = plotPRFs(opts_sl, sites, pRFs_sl, fH);
axis([-20 35 -20 35])


%% Figure 2C: Zoom in to each channel
fH = figure;  pos = get(fH, 'Position');  pos(3) = 1000; pos(2) = 400;
set(fH, 'Position', pos, 'Color', 'w','name','Winawer and Parvizi, Figure 2C', 'NumberTitle', 'off')
for ii = 1:length(sites)
    subplot(2,length(sites),ii);
    fH   = plotPRFs(opts_phosp, sites(ii), phos(:,:,ii), fH);
    fH   = plotPRFs(opts_bb, sites(ii), pRFs_bb(:,:,ii), fH);
    axis(getZoom(sites(ii)))
    
    subplot(2,length(sites),ii+length(sites));
    fH   = plotPRFs(opts_phosp, sites(ii), phos(:,:,ii), fH);
    fH   = plotPRFs(opts_sl, sites(ii), pRFs_sl(:,:,ii), fH);
    axis(getZoom(sites(ii)))
end

%% Figure 2D and S3B: Plot overlap coefficients
fH = figure;  pos = get(fH(1), 'Position');  pos(3) = 1000; pos(4) = 200;
set(fH(1), 'Position', pos, 'Color', 'w', 'name','Winawer and Parvizi, Figure 2D', 'NumberTitle', 'off');

fH(2) = figure;  pos = get(fH(2), 'Position');  pos(3) = 1000; pos(4) = 200;
set(fH(2), 'Position', pos, 'Color', 'w', 'name', 'Winawer and Parvizi, Figure S3B', 'NumberTitle', 'off');

fH = plotDice(pRFs_bb, pRFs_sl, phos_all, phos, fH);


end

function fH = plotDice(pRFs_bb, pRFs_sl, phos_all, phos, fH)

numchannels = size(pRFs_bb,3);

% Compute overlap between phosphenes and pRFs
nsigmas = 2;
[o_bb, o_all_bb] = computeOverlap(pRFs_bb, phos, phos_all, nsigmas);
[o_sl, o_all_sl] = computeOverlap(pRFs_sl, phos, phos_all, nsigmas);

% ttests - dice coefficient of broadband/ebs v stimlocked/ebs
p_ttest = NaN(1,numchannels);
p_boot  = NaN(1,numchannels);
figure(fH(1));
set(gca, 'FontSize', 18); hold on
for ii = 1:numchannels 
    [~, p_ttest(ii)] = ttest(o_all_bb{ii} - o_all_sl{ii});
    d = o_all_bb{ii} - o_all_sl{ii};
    bootstat = bootstrp(10000, @mean, d);
    p_boot(ii) = mean(bootstat < 0);
    
    subplot(1,5, ii);
    set(gca, 'FontSize', 12); hold on
    plot([0 0], [0 400], 'k--')
    histogram(bootstat, 100); 
    xlim([-.5 .5]); set(gca, 'XTick', -5:.25:.5)
    
       title(sprintf('Site %d', ii));
        if ii == 3
            xlabel('Dice Coefficient between Phosphenes and pRFs, Broadband minus Stimulus-Locked')
        end
        if ii == 1
            ylabel('Number of bootstraps')
        end
    set(gca, 'XTickLabel', [], 'YTickLabel', [])
end


fprintf('bb o mean\t'), disp(o_bb)
fprintf('sl o mean\t'), disp(o_sl) 
fprintf('bb mean o\t'), disp(cellfun(@mean, o_all_bb)) 
fprintf('sl mean o\t'), disp(cellfun(@mean, o_all_sl))
fprintf('p-val ttest\t'), disp(p_ttest) 
fprintf('p-boot\t'), disp(p_boot) 
fprintf('p-val bootstrp\t'), disp(1 - (abs(.5 - p_boot) * 2)) 



ns = 0.5:0.5:4; idx = ns == 2;
o_bb = zeros(length(ns),numchannels); 
o_sl = zeros(length(ns),numchannels); 
o_all_bb = o_bb;
o_all_sl = o_bb;

for ii = 1:length(ns)
    [o_bb(ii,:), tmp_bb] = computeOverlap(pRFs_bb, phos, phos_all, ns(ii));
    [o_sl(ii,:), tmp_sl] = computeOverlap(pRFs_sl, phos, phos_all, ns(ii));
    o_all_bb(ii,:) = cellfun(@mean, tmp_bb);
    o_all_sl(ii,:) = cellfun(@mean, tmp_sl);
end


figure(fH(2));
for ii = 1:numchannels
    subplot(2,numchannels,ii), set(gca, 'FontSize', 12), hold on
    plot(ns, o_bb(:,ii), 'b',  ns, o_sl(:,ii), 'r', ...
         ns(idx), o_bb(idx,ii),'bo', ns(idx), o_sl(idx,ii),'ro', ...
         'LineWidth',4)
     %ylabel('Dice coefficient')
     set(gca, 'YLim', [0 .8], 'YTick', [0 .4 .8], 'FontSize', 18)

    subplot(2,numchannels,ii+numchannels), set(gca, 'FontSize', 12), hold on
    plot(ns, o_all_bb(:,ii), 'b',  ns, o_all_sl(:,ii), 'r', ...
         ns(idx), o_all_bb(idx,ii),'bo', ns(idx), o_all_sl(idx,ii),'ro', ...
         'LineWidth',4)
     %xlabel('PRF threshold (SDs)')
     %ylabel('Dice coefficient')
     set(gca, 'YLim', [0 .45], 'YTick', [0 .2 .4], 'FontSize', 18)
end

end

function fH = plotPRFs(opts, sites, pRFs, fH)

rnge = -opts.rad+opts.srate/2:opts.srate:opts.rad-opts.srate/2;

if ~exist('fH', 'var'), fH = figure; end
fH = plotPolar(opts.rad, fH);

colors = getColors(sites);

switch opts.plottype
    case 'shaded'
        im = reshape(pRFs, size(pRFs,1)*size(pRFs,2),[]);
        im = im*(1-colors);
        im = 1 - im;
        im = reshape(im, [size(pRFs,1), size(pRFs,2), 3]);
        imagesc(rnge, rnge, im);
        fH = plotPolar(opts.rad, fH);
        
    case 'contours'
        for ii = length(sites):-1:1
            
            contourHeights = exp(-opts.whichsigmas.^2/2);
            
            % If we are plotting a Gaussian of height 1, then we plot
            % contour lines height contourHeights. But if our data are not
            % Gaussian, then we we need to correct for the non-Gaussianity
            % and find the contour which would enclose the equivalent area            
            thisprf  = pRFs(:,:,ii);
            sortedRF = sort(thisprf(:));
            tmp      = cumsum(sortedRF); % cummulative sum of the sorted image is like the area                      
            for jj = 1:length(contourHeights)
                [~, idx(jj)] = min(abs(tmp-contourHeights(jj)*max(tmp)));
            end
            contourHeights = sortedRF(idx);
            
            
            contour(rnge,rnge, pRFs(:,:,ii),contourHeights, '-', 'LineWidth', 2, 'Color', colors(ii,:));
            contour(rnge,rnge, pRFs(:,:,ii),contourHeights, '--', 'LineWidth', 1, 'Color', 'k');
            
            %text(x(ii), y(ii), num2str(sites(ii)), 'FontSize', 18);
            
        end
end


end

function [pRFs, r] = getPRFdata(opts, sites)


rnge = -opts.rad+opts.srate/2:opts.srate:opts.rad-opts.srate/2;
[X,Y] = meshgrid(rnge);

for ii = 1:length(sites)
    params = ebs_solvePRFmodels(sites(ii));
    if opts.sl, params = params.sl; else, params = params.bb; end
            
    x0(ii) = params.params(1);
    y0(ii) = params.params(2);
    sigma(ii) = params.params(3) ./ sqrt( params.params(5));
    r(ii) = params.r;
end
[x, y, s] = ecogFitPRFPix2Deg(1:length(sites), x0, y0, sigma);

pRFs = rfGaussian2d(X(:), Y(:), s, s, s*0, x, y);
pRFs = reshape(pRFs, round(sqrt(size(pRFs,1))), round(sqrt(size(pRFs,1))), []);



end

function [phos, phos_all] = getPhosdata(opts, sites)

rnge = -opts.rad+opts.srate/2:opts.srate:opts.rad-opts.srate/2;
[X,Y] = meshgrid(rnge);

phos   = zeros(length(rnge), length(rnge), length(sites));
phos_all = cell(1,length(sites));

pth = fullfile(ebsRootPath, 'data', 'ebs');

for ii = 1:length(sites)
    load(fullfile(pth, sprintf('phosphenes_site%d', ii)));
    trial_data=load(fullfile(pth, sprintf('trial_data_site%d', ii)));
    
    valid_trials = isfinite(trial_data.which_drawing.val);
    these_drawings = trial_data.which_drawing.val(valid_trials);
    
    IN.grid_points = false(length(these_drawings), numel(X));
    
    for jj = 1:length(these_drawings)
        if isempty(phosphenes(these_drawings(jj)).x)
            % skip the drawing if it's empty
        else
            IN.grid_points(jj,:) = inpolygon(X(:),Y(:),...
                phosphenes(these_drawings(jj)).x,phosphenes(these_drawings(jj)).y);
        end
    end
    
    im = mean(IN.grid_points);
    im = im/max(im);
    
    phos(:,:,ii) = reshape(im, size(X));
    phos_all{ii} = IN.grid_points;
end

end

function [r, r_all] = computeOverlap(pRFs, phosphenes, phos_all, nsigmas)

dice = @(a,b) 2*sum(a(:)&b(:))/(sum(a(:))+ sum(b(:))); 
jaccard = @(a,b) sum(a(:) & b(:)) / sum(a(:) | b(:));

overlap = dice;

% pRfs is x by y by sites 
n = size(pRFs,3);
r = NaN(1,n);
r_all = cell(1,n);

for ii = 1:n
    
    contourHeights = exp(-nsigmas.^2/2);
    thisprf  = pRFs(:,:,ii);
    sortedRF = sort(thisprf(:));
    tmp      = cumsum(sortedRF); % cummulative sum of the sorted image is like the area
    [~, idx] = min(abs(tmp-contourHeights*max(tmp)));    
    contourHeights = sortedRF(idx);
    prf_thresh = thisprf(:) > contourHeights;
    
    contourHeights = exp(-nsigmas.^2/2);
    thisprf  = phosphenes(:,:,ii);
    sortedRF = sort(thisprf(:));
    tmp      = cumsum(sortedRF); % cummulative sum of the sorted image is like the area
    [~, idx] = min(abs(tmp-contourHeights*max(tmp)));    
    contourHeights = sortedRF(idx);
    phos_thresh = thisprf(:) > contourHeights;
    
    r(ii) = overlap(prf_thresh, phos_thresh);
    
    for jj = 1:size(phos_all{ii},1)                
        r_all{ii}(jj) = overlap(prf_thresh, phos_all{ii}(jj,:));
    end
    
end


end

function fH = plotPolar(rad, fH)

if ~exist('fH', 'var'); fH = figure; end

plot_prefs.linewidth = 1;
plot_prefs.gridcolor = .5*[1 1 1];

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

axis(rad*[-1 1 -1 1]), axis square

end

function saveFig(fH, fname, save_flag)

if ~save_flag, return; end

hgexport(fH, fullfile(ebsRootPath, 'figures', [fname '.eps']));
saveas  (fH, fullfile(ebsRootPath, 'figures', [fname '.fig']));

end