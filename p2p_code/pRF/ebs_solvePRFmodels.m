function [params, resp] = ebs_solvePRFmodels(site)
% We compute  two pRF models per site, one for the broadband signal and one
% for the stimulus-locked signal. Notation:

% Windowing epochs for the Fourier spectra
useHann = true; % if true use a hann window, otherwise a rectangular window

% --- Broadband pRF ------------------------------------------------------
calcPower = true;   % power rather than amplitude
useExp    = true;   % compressive spatial summation model 
%                     (See Winawer et al, 2013, Current Biology for details).
[params.bb, resp.bb] = ecogFitPRF(site, useExp, calcPower, 'bb', useHann);

% --- Stimulus-locked pRF ------------------------------------------------
calcPower        = false;  % amplitude rather than power
useExp           = false;  % linear pRF model
subtractBBfromSS = true;   % for computing the stimulus locked resposne, 
%                            we subtract the modeled broadband signal at 
%                            the relevant frequency
[params.sl, resp.sl] = ecogFitPRF(site, useExp, calcPower, 'ss2f', useHann, subtractBBfromSS);

end