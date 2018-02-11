function [ts, srate, nruns] = ebsGetPRFTSeries(site)
% ts is an array time x epoch. srate is the sample rate in Hz. nruns is the
% number of repeated experiments with the same stimulus.

% path to stimulus
pth   = fullfile(ebsRootPath, 'data', 'ecogPRF');
fname = sprintf('prfDataEpoched_site%d', site);

% load it
tmp = load(fullfile(pth, fname));

ts      = tmp.ts;
srate   = tmp.srate;
nruns   = tmp.nruns;

return