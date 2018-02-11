function [tseries, freq, ph] = ecogGetSpectralTSeries(site, ...
    calcPower, datatype, useHann, subtractBBfromSS)
% [tseries, freq, phase] = ...
%   ecogGetSpectralTSeries(site, [calcPower=false], ...
%       datatype,[useHann=true],[subtractBBfromSS=true] )
%
% Example:
%   site = 1;
%   calcPower  = true; % return amplitude, not square amplitude (power)
%   datatype  = 'bb';
%   [ts freq]  = ecogGetSpectralTSeries(site, calcPower, datatype);


%% CHECK INPUTS

if notDefined('calcPower'), calcPower = false;  end % by default, we calculate signal amplitude, now power
if notDefined('fmax'),      fmax = 150;         end % by default, up to 150 Hz
if notDefined('useHann'),   useHann = true;     end % by default, use Hann window for FFT
if notDefined('subtractBBfromSS'), subtractBBfromSS = true; end % by default, discount bb fit when calculating ss amplitude


[S, Ph, f, T] = ecogGetSpectralData(site, calcPower, fmax, useHann);

[ss, bb, ssi, bbi] = ecogGetSSandBBFrequencies(f, T);

% we get up to three types of time series, and put each into the variable
% tseries (3 x nscans). to keep track of what we are doing, we name our
% three types of time series
%
% Models can be
%   'ss'
%   'bb line'
%   'bb pooled'




switch lower(datatype)
    case {'ss' 'ss2f'}
        % steady state (1st or 2nd harmonic of the reversal frequency, as requested)
        if strcmpi(datatype, 'ss'),       ssi = ssi(1); ss = ss(1);
        elseif strcmpi(datatype, 'ss2f'), ssi = ssi(2); ss = ss(2);
        end
        
        tseries = sum(S(:,ssi), 2);
        ph      = Ph(:,ssi(1));
        
        % if we want to subtract the BB from the SS
        if subtractBBfromSS
            tseriesBB  = ecogGetSpectralTSeries(site, calcPower, 'bb');
            tseries    = tseries - tseriesBB;
        end
        
    case {'bb'}
        % Broad band linear fit
        
        % fit the spectral data (after removing SS and harmonics) with a line in
        %   log-log space
        
        % Take the log of the signal and the frequencies to get log-log space
        Y = log(S(:,bbi));
        X = repmat(log(bb), size(Y, 1), 1);
        
        % fit ONE line to all the concatenated trials
        p = polyfit(X(:), Y(:), 1);
        
        % calculate the fitted broad band amplitude at the steady state
        % frequency
        ssa = polyval(p, log(ss(1)));
        
        % calculate the residuals from the linear prediction
        residual = bsxfun(@minus, Y, polyval(p, log(bb)));
        
        % take the mean of the residuals from each trial and add the interecept to
        % the single linear fit to get the  intercept for that trial
        intercepts = mean(residual,2) + ssa;
        
        % broad band based on line fit
        tseries = exp(intercepts);
        
        % because the broadband is comprised of many frequencies, the
        % phase has no specific meaning.
        ph   = NaN(size(tseries));
        
end


% return some information about the frequencies used for calculations
freq.f   = f;
freq.ss  = ss;
freq.bb  = bb;
freq.ssi = ssi;
freq.bbi = bbi;


return