function [x, y, s] = ecogFitPRFPix2Deg(site, x0, y0, s0)
%Convert pixels to degrees. 
%
% [x, y, s] = ecogFitPRFPix2Deg(site, x0, y0, s0)
%
% We compute PRF models in pixels, but would like to report them in
% degrees.
%
% Inputs:   
%   site: site number (1-5)
%   x0:   x center of pRF in pixels
%   y0:   y center of pRF in pixels
%   s0:   pRF size in pixels (pRF size = sigma/sqrt(n))
%
% Outputs
%   x:  x center of pRF in degrees of visual angle (0 is fixation)
%   y:  y center of pRF in degrees of visual angle (0 is fixation)
%   s:  pRF size in degrees

% Make sure all input arguments are row vectors

nPoints = max([length(x0) length(y0) length(s0)]);
if isempty(x0), x0 = zeros(1,nPoints); end
if isempty(y0), y0 = zeros(1,nPoints); end
if isempty(s0), s0 = zeros(1,nPoints); end

if length(site) == 1, site = repmat(site, nPoints, 1); end
x0 = x0(:)'; y0 = y0(:)'; s0 = s0(:)'; site = site(:)';

viewing_distance = zeros(size(site));
screen_height    = zeros(size(site));

%  viewing distance, screen height in cm (different for each subject)
for ii = 1:length(site)
    switch site(ii)
        case {3, 5}, viewing_distance(ii) = 47;  screen_height(ii) = 20.7; % cm
        case 1,     viewing_distance(ii) = 30;  screen_height(ii) = 17.9; % cm
        case 4,     viewing_distance(ii) = 45;  screen_height(ii) = 17.9; % cm
        case 2,     viewing_distance(ii) = 50;  screen_height(ii) = 17.9; % cm
        otherwise
            error('Unknown viewing distance for subject %d', site);
    end
end

% height from center of screen to top of screen
radius_in_cm  = screen_height/2; % radius in cm
radius_in_deg = rad2deg(atan(radius_in_cm./viewing_distance)); % radius in degrees

numpix  = ones(size(site))*101; % number of pixels in stimulus description (one side) as input to fitprf
pix2deg = radius_in_deg*2./numpix;
xC      = (numpix+1)/2; % center location in pixels
yC      = xC;            
xC(site==1) = 0; 

% convert pixels to degrees
x =  pix2deg.*(x0-xC);
y = -pix2deg.*(y0-yC);
s =  pix2deg.*s0;

end