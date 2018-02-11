function ebs_MakeFigure1(site)
% Reproduce Figure 1 for the following paper:
%  
%   Winawer and Parvizi (2016). Linking Electrical Stimulation of Human
%   Primary Visual Cortex, Size of Affected Cortical Area, Neuronal
%   Responses, and Subjective Experience Neuron. 92(6):1213?1219
%   http://dx.doi.org/10.1016/j.neuron.2016.11.008
%
% For Figure 1, use site 2 from paper:
% ebs_MakeFigure1(1)
%
% For Supplementary Figure 1, loop through all 5 sites
% for site = 1:5; ebs_MakeFigure1(site); end

if notDefined('site'), site = 2; end

% Solve broadband and stimulus locked ECoG pRF models
[params, resp] = ebs_solvePRFmodels(site);

% Identify epochs with high versus low responses
epochs = ebs_identfyEpochs(site, resp);

% Visualize electrode location and V1 maps
ebs_Figure1AB();

% Set up plot
[fH, plotOpt] = ebs_setUpPlot();

% Plot time series averaged across epochs (Figure 1C, top)
ebs_Figure1Ctop(fH, plotOpt, site, epochs);

% Plot amplitude spectrum, averaged across epochs (Figure 1C, bottom)
ebs_Figure1Cbottom(fH, plotOpt, site, epochs, resp);

% Plot component broadband and stimulus locked time series (Figure 1D)
ebs_Figure1D(fH, plotOpt, epochs, resp, params);

end

function ebs_Figure1AB()

numSubjects = 4;

% Get two images per subject, one with the V1-V3 Benson ROIs, and one with
% the Hinds et al probabilistic atlas, and store the images in the variable
% <im>
im = cell(numSubjects, 2);

for ii = 1:numSubjects
    thisSubject = numSubjects - ii + 1;
    ebs_visualizeElectrodes(thisSubject, 'areas')
    axis off; grid off;   
    [im{ii,1}, ~] = getframe(gcf);
    close(gcf);
    
    ebs_visualizeElectrodes(thisSubject, 'v1')
    axis off; grid off; 
    [im{ii,2}, ~]  = getframe(gcf);
    close(gcf);
end

% Combine the images into one figure to show the panels for Figure1A and
% Figure1B
figure,set(gcf, 'Color', 'w', 'Name', 'Winawer & Parvizi (2016) Figure 1A&B', 'NumberTitle', 'off')

% These are the crop points to highlight just the occipital lobe
crop = [...
    200 600 550 875; ...
    200 600 550 875 ; ...
    150 550 550 875;...
    250 650 150 475;  ];

% Position of the first subplot in the figure [left, bottom, width, height]
pos = [0 0 .5 .7]/2;

% Mask to make soft white edge around the part of the brain that is masked
numrows = 401;
numcols = 326;
thresh  = 0.8;
mask = 255* ones(numrows,numcols);
[x, y] = meshgrid(linspace(-1,1,numcols), linspace(-1,1,numrows));
x = abs(x); y = abs(y);
weights = max(x,y);
weights(weights <= thresh) = 0;
weights(weights > thresh) = weights(weights > thresh) - thresh;
weights = weights / max(weights(:));

% Loop across subjects (4) and plot types (2 - V1-V3 templates and Hinds V1
% probabilistic map)
for col = 1:numSubjects
    for rows = 1:2
        % Benson atlas failed for subject 1 (col 4), so we don't plot it 
        if col == 4 && rows == 1, continue; end
        
        % New position for each subplot, dividing the figure into fourths
        pos(1) = (4-col)*.25;
        pos(2) = .05+(2-rows)*.5;
        
        % Crop the image and blur edges (near the crop) in each RGB channel
        subplot('Position', pos);
        thisim = double(im{col,rows});
        thisim = thisim(crop(col,1):crop(col,2), crop(col,3):crop(col,4),:);
        
        for ii=1:3
            thisim(:,:,ii) = weights.* mask + (1-weights) .* thisim(:,:,ii);
        end
        
        % Show it
        thisim = uint8(thisim);
        imshow(thisim)
        axis off
        title(sprintf('Subject %d', numSubjects - col+1));
    end
end

end

function [fH, plotOpt] = ebs_setUpPlot()

fH = figure; pos = get(fH, 'Position');
pos(3:4) = [1200 675]; set(fH, 'Position', pos, 'Color', 'w', ....
    'Name', 'Winawer & Parvizi (2016) Figure 1C&D', 'NumberTitle', 'off')

% Plot options
plotOpt.fontsz      = 20;
plotOpt.markersz    = 8;
plotOpt.fmax        = 200;
plotOpt.colors      = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
plotOpt.calcPower   = false; % plot amplitude, not power
plotOpt.useHann     = false; % rectangular window, not hann window
end

function epochs = ebs_identfyEpochs(site, resp)
% Identify epochs with high versus low responses

% Load the mapping stimulus (epochs x pixels; pixels have been converted
% from 101x101 matrix to 1021 x 1 vector for computational efficiency)
stimulus    = ecogGetPRFStimulus(site);

% Which epochs (if any) had blanks (ie no stimulus contrast)?
blanks = mean(stimulus, 2) < 10e-6;

% Find the 24 epochs with the largest response (averaged across broadband
% and stimulus locked, after normalizing each)
r1 = resp.bb/max(resp.bb);
r2 = resp.sl/max(resp.sl);
[~, idx] = sort(r1+r2, 'descend');
epochs.stimulus = false(size(resp.bb));
epochs.stimulus(idx(1:24)) = 1;

% Find the epochs with blanks (if any) or (if no blanks) find the 24 epochs
% with the lowest response (averaged across broadband and stimulus locked,
% after normalizing each)
if sum(blanks) == 0
    epochs.blanks = false(size(resp.bb));
    epochs.blanks(idx(end-23:end)) = 1; 
else  
    epochs.blanks = blanks; 
end

% Epochs with high and low responses
[~, ~, nruns] = ebsGetPRFTSeries(site);
epochs.stimulus_all = repmat(epochs.stimulus, [nruns, 1]);
epochs.blank_all    = repmat(epochs.blanks, [nruns, 1]);

end

function ebs_Figure1Ctop(fH, plotOpt, site, epochs)

figure(fH);

subplot(3,2,1)

% Return the epoched time series *without* averageing across runs
[ts, srate] = ebsGetPRFTSeries(site);

% Transpose so that ts is epochs x time
ts = ts';

% Time (in seconds)
t = (1:size(ts,2))/srate;

% Mean and SEM response to blank (or low-response epochs, if no blanks)
mn = mean(ts(epochs.blank_all,:));
se = std(ts(epochs.blank_all,:))/sqrt(sum(epochs.blank_all));

% Plot line and shading to indicate mean and standard error
fill([t flip(t)], [mn+se flip(mn-se)], [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5)
hold on,
plot(t, mn, '-', 'LineWidth', 3, 'Color', [0.7 0.7 0.7]);

% Mean and SEM response to stimulus (high responses epochs)
mn = mean(ts(epochs.stimulus_all,:));
se = std(ts(epochs.stimulus_all,:))/sqrt(sum(epochs.stimulus_all));

% Plot line and shading to indicate mean and standard error
fill([t flip(t)], [mn+se flip(mn-se)], 'w', 'EdgeColor', 'k' ,'FaceAlpha', 0.5)
plot(t, mn, 'k-', 'LineWidth', 3);
xtick_labels = {' ', ' ', ' ', '0.2', ...
    ' ', ' ', '0.4', ...
    ' ', ' ', '0.6', ...
    ' ', ' ', '0.8', ...
    ' ', ' ', '1.0'};
set(gca, 'XTick', (0:15)/15*max(t), 'XGrid', 'on', 'xlim', [4/15 1], ...
    'FontSize', plotOpt.fontsz , 'box', 'off', ... 'ylim' ,[-75 50],...
    'XTicklabel', xtick_labels );
ylabel('Amplitude (µV)'); xlabel('Time (s)')

end

function ebs_Figure1Cbottom(fH, plotOpt, site, epochs, resp)

figure(fH);

subplot(3,2,[3 5]);
set(gca, 'FontSize', plotOpt.fontsz); hold on

% Compute amplitude spectra for all epochs
[A, ~, f] = ecogGetSpectralData(site, plotOpt.calcPower, plotOpt.fmax, plotOpt.useHann);

% Reshape and average across repeated runs 
A = reshape(A, length(resp.bb),[], size(A,2));
A = squeeze(mean(A,2));

% Average across frequencies for stimulus epochs and for blanks (in each
% case, first log, then average, then exponentiate)
mn_Amp.stimulus = exp(mean(log(A(epochs.stimulus,2:length(f)))));
mn_Amp.blanks = exp(mean(log(A(epochs.blanks,2:length(f)))));

% Standard error across frequencies for stimulus epochs and for blanks
sem_Amp.stimulus = std(A(epochs.stimulus,2:length(f)))/sqrt(sum(epochs.stimulus));
sem_Amp.blanks = std(A(epochs.blanks,2:length(f)))/sqrt(sum(epochs.blanks));

% Omit DC for plotting
f = f(2:end);

% Plot error bars by filling mean +/- sems
fill([f flip(f)], [mn_Amp.stimulus + sem_Amp.stimulus flip(mn_Amp.stimulus - sem_Amp.stimulus)], 'w', 'EdgeColor', 'k' ,'FaceAlpha', 0.5, 'EdgeAlpha', .5)
plot(f, mn_Amp.stimulus, 'Color', 'k', 'LineWidth', 1)

fill([f flip(f)], [mn_Amp.blanks + sem_Amp.blanks flip(mn_Amp.blanks - sem_Amp.blanks)], [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeAlpha', .5)
plot(f, mn_Amp.blanks, 'Color', .7*[1 1 1], 'LineWidth', 1)

% Plot options
set(gca, 'XSCale', 'log', 'XLim', [10 max(f)]);
set(gca, 'YSCale', 'log');
yl = get(gca, 'YLim');
hold on

% Plot means
plot(f, mn_Amp.stimulus, 'Color', 'k', 'LineWidth', 2)
plot(f, mn_Amp.blanks, 'Color', .7*[1 1 1], 'LineWidth',2)

% Grid lines
for zz = 1:floor(plotOpt.fmax/15)
    plot(f(zz*15) * [1 1], yl, '-', 'Color', .8*[1 1 1], 'LineWidth', 1)
end

% Labels
xlabel('Frequency (Hz)'); ylabel('Amplitude (µV)')
set(gca, 'XTick', 15:15:200, 'XTickLabel', {'15', '30', '45', '60', '75', ' ', '105', ' ', ' ', '150', ' ', ' ', ' '})


%% Fit slope of power law (for blank condition)

% Exclude frequencies below 10 Hz or near multiples of 60
idx = abs(f-60) < 2 | abs(f-120) < 2 | abs(f-180) < 2 | f < 10;


x = log(f);         % log frequency
y = log(mn_Amp.blanks.^2);  % low power during blank epochs

x(idx) = NaN;

lm = fitlm(x, y); % linear model fit
m = lm.Coefficients.Estimate(2); % slope of power law

text(15, exp(mean(log(yl)/4)), sprintf('Power ~ f^%4.2f', m), 'FontSize', 18, 'Interpreter', 'none');

disp(m)
% site 1: -2.8988
% site 2: -2.2389
% site 3: -2.8670
% site 4: -2.4472
% site 5: -2.6433

end

function ebs_Figure1D(fH, plotOpt, epochs, resp, params)
figure(fH)
% Loop over 2 data types (broadband and stimulus locked)
for dt = 1:2
    if dt == 1
        thisresp = resp.bb; 
        thispar  = params.bb;
        plotnum  = 2;
    else
        thisresp = resp.sl; 
        thispar  = params.sl;
        plotnum = 6;
    end
    
    % num time points
    nt = length(thisresp);
    t = 1:nt;

    subplot(3,2,plotnum)
    set(gca, 'FontSize', plotOpt.fontsz, 'XTick', 0:12:nt, 'XLim', [0 nt], 'XGrid', 'on'); hold on
    y = thisresp-thispar.drift;
    plot(t, y , 'o', 'Color', plotOpt.colors(dt,:), 'LineWidth', 2, 'MarkerSize', plotOpt.markersz)
    plot(t(epochs.stimulus), y(epochs.stimulus) , 'ok', 'LineWidth', 2,'MarkerSize', plotOpt.markersz, 'MarkerFaceColor', 'w')
    plot(t(epochs.blanks), y(epochs.blanks) , 'o', 'Color', .7*[1 1 1], 'LineWidth', 2,'MarkerSize', plotOpt.markersz,'MarkerFaceColor', 'w')

    plot(t, thispar.signal , '-', 'Color', plotOpt.colors(dt,:), 'LineWidth', 2)
    xlabel('Time (s)');
    if dt == 1, ylabel('Power (µV^2)'); else, ylabel('Amplitude (µV)'); end

    axis tight;

end
end