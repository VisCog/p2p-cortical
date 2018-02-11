function [ss, bb, ssi, bbi] = ecogGetSSandBBFrequencies(f, T, stimF)
% [ss, bb, ssi, bbi] = ecogGetSSandBBFrequencies(f, T, stimF)

% stimulus reversal rate
if notDefined('stimF'), stimF = 15; end

% steady state frequencies 
%   (1st and 2nd harmonics of stimulus reversal rate)
[~, ssi]    = min(abs(f - stimF/T));
[~, ssi(2)] = min(abs(f - 2*stimF/T));

ss = f(ssi);

%% broadband frequencies - 
%  all frequencies except steady state harmonics, line noise harmonics, and
%  low frequncies (< 7 Hz)

% drop the steady state and harmonics (frequencies within 3 Hz of multiples
%   of the ss frequency, close to stimF Hz)
tmp     = (1:10)*stimF;
ssdrop  = sort([tmp-3 tmp-2 tmp-1 tmp tmp+1 tmp+2 tmp+3]);

% drop line noise(frequencies within 3 Hz of multiples of 60 Hz)
tmp = (1:5) * 60;
lndrop   = sort([tmp-3 tmp-2 tmp-1 tmp tmp+1 tmp+2 tmp+3]);

% low frequency drop
lfdrop = 0:7;

% remaining frequencies are used to compute broadband
[~, bbi]   = setdiff(round(f*T), [ssdrop lndrop lfdrop]);
bb = f(bbi);

%% Done
