function [A, P, f, T, ts] = ecogGetSpectralData(site, calcPower, fmax, useHann)
% Calculate the spectrogram via short-time fourier transform of an ECoG
% experiment.
%
%   [A, P, f, T, ts] = ecogGetSpectralData(site, calcPower, fmax, useHann)
%
% Inputs
%  site: site 1-5 (integer scalar)
%  calcPower: boolean. if true, then A = squared amplitude
%  fmax: maximum temporal frequency (scalar)
%  useHann: boolean. if true, then use a hanning window for spectra
%
% Return
%   A:      amplitude spectrum (trials x frequencies)
%   f:      row vector of frequencies
%   P:      phase spectrum (trials x frequencies)
%   T:      epoch length, in seconds
%   ts:     time series (time points x epoch)
%% Set up

if notDefined('calcPower'),   calcPower   = false; end
if notDefined('fmax'),        fmax        = 150;   end
if notDefined('useHann'),     useHann     = true;  end

%% Get the epoched data (at sampling rate srate)
[ts, srate, nruns] = ebsGetPRFTSeries(site);

% time series should be time points by trials 

nt = size(ts,1);    % time points per epoch
T  = nt/srate;      % epoch length (in seconds)

%% Short-time fourier transform

% window 
if useHann, H = hann(nt); else, H = ones(nt,1); end

% frequencies of interest (omit frequencies near line noise and harmonics)
f = 0:fmax;

% The frequencies need to be adjusted to account for the display refresh
% rate. Because the rate is close to, but not precisely 60 Hz, the time
% that a stimulus aperture is displayed is close to, but not precisely 1 s.
% the actual time window (T) is usually slightly more than 1. we adjust the
% frequencies down by the same amount. hence the steady state frequency is
% not 15 Hz, but rather 15/T Hz.
f = f/T;

% window the time series
ts = bsxfun(@times, ts, H);

% fourier transform
FT    = fft(ts);

% separate into amplitude and phase (scaling amplitude by 2/nt)
A = abs(FT)/nt*2;
P = angle(FT);

% if requested, square the amplitude 
if calcPower,   A = A.^2; end

% transpose for compatibility subsequent calculations
A = A';
P = P';

return
