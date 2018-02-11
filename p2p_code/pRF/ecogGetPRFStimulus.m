function stimulus = ecogGetPRFStimulus(site)
% Stimulus is a cell array of stimulus descriptions. Each stimulus
% description is pixels x time.

% path to stimulus
pth   = fullfile(ebsRootPath, 'stimuli');
fname = sprintf('bar_site%d', site);

% load it
load(fullfile(pth, fname));

% reshape to a matrix (pixels x time)
stimulus = reshape(images, size(images, 1)*size(images, 2), [])';

return