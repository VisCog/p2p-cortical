function x = fpDefaultInputs(stimulus)

if ~exist('stimulus', 'var') 
    stimulus = []; 
    warning('No stimulus information. Will not be able to solve model') %#ok<*WNTAG>
end

% stimulus
x.stimulus = stimulus;

% function handle to pull time series out of data file
x.response        = [];

% 2d gaussian pRF model
if isempty(stimulus), x.prfmodel = [];
else
    x.prfmodel = constructprfmodel2Dgaussian(...
        cellfunfirst(@(x)sqrt(size(x,2)),stimulus)); 
end


% hrf: by default we use delta function for ECoG mdoels
x.hrfmodel        = 1;

% enforce separability between PRF and HRF
x.flag            = 1;    

% polynomial degree for detrending. knk thinks 2 is good based on the length of these scans (~3.5 min)
x.maxpolydeg      = 2;

% no longer used !?
x.ar              = [];

% allow exponent in prf model. {1 A B} means allow exponent.  A is the
% initial seed for the exponent (1 x 1). B are the bounds (2 x 1)
x.mode            = {1 0.5 [0; Inf]}; 

% how long to let matlab search for a good fit (using lsqnonlin)
x.maxiter         = 2000; % Inf

% here we deal with average data type (one scan) so we don't bootstrap
x.wantresample    = 0;

% tolerance for matlab lsqnonlin
x.tol             = 1e-06;

% parameters for lsqnonlin call
x.extraopt        = {'LargeScale' 'off', 'Display' 'final' 'Algorithm' 'levenberg-marquardt'};

% ok. 
x.extraregressors = [];

% huh?
x.hrfnormfun      = @(x)calchrfpeak(x,1);

% huh?
x.derivemode      = [];

% huh?
x.metric          = [];

