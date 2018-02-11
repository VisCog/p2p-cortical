function [params, tseries] = ecogFitPRF(site, useExp, ...
    calcPower, datatype,  useHann, subtractBBfromSS)
%[params, response] = ecogFitPRF(site, [useExp=true], ...
%    [calcPower=false], [datatype='bb'],...
%    [useHann=true], [subtractBBfromSS=true])

if notDefined('useExp'),            useExp = true;           end
if notDefined('calcPower'),         calcPower = false;       end
if notDefined('datatype'),          datatype = 'bb';         end
if notDefined('useHann'),           useHann = true;          end
if notDefined('subtractBBfromSS'),  subtractBBfromSS = true; end

%% Define input tseries

% stimulus is time (epoch #) x pixels
stimulus = ecogGetPRFStimulus(site);

% tseries is 1 x n time points, concatenated across repeated experiments 
% with the same stimulus. 
tseries  = ecogGetSpectralTSeries(site, calcPower, datatype, useHann, subtractBBfromSS);

% average tseries across repeated runs
tseries = reshape(tseries, 96, []);
tseries = mean(tseries,2);
[tseries, blanks] = ecogSubtractBlanks(tseries, stimulus);

%% Define all other fitprf inputs

% prfmodel
x = fpDefaultInputs(stimulus);
prfmodel = x.prfmodel; % 2-d gaussian (with static non-linearity)

% hrfmodel
hrfmodel = 1; % impulse response function

% flag
flag = 1; % NA

% maxpolydeg
if sum(blanks) > 0, maxpolydeg = NaN;    % no detrending
else,               maxpolydeg = 0;  end % model the mean

% ar
ar = []; % NA

% mode
if useExp,  mode = {1 1 [0; Inf]};
else,       mode = 0;  end

% maxiter
maxiter = 2000;

% wantresample
wantresample = 0; % for cross-validation

% tol
tol = x.tol;

% extraopt
extraopt = x.extraopt;
extraopt{4} = 'final';%'iter';

% extraregressors
extraregressors = [];

% hrfnormfun
hrfnormfun = [];

% derivemode
derivemode = [];

% metric
if sum(blanks) > 0, metric = @(x,y)calccod(x,y,[],[],0); % if we know the baseline
else,               metric = @(x,y)calccod(x,y,[],[],1); end

% outputfcn
outputfcn = [];

%% SOLVE!!

% run the prf fit!
params = fitprf(stimulus,tseries,prfmodel,hrfmodel,flag,maxpolydeg,ar,mode,maxiter, ...
    wantresample,tol,extraopt,extraregressors,hrfnormfun,derivemode,metric,outputfcn);
