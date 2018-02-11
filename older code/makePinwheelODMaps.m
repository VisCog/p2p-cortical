function [pinwheel,OD] = makePinwheelODMaps(x,y,sig, ODsize)
% [pinwheel,OD] = makePinwheelODMaps(x,y,sig)
%
% sig determines the distribution of OD values. Default is 5.  The larger
% sig, the more the distribution tends toward 0 and 1.

%initial parameters:

if ~exist('sig','var')
    sig = .5; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex.
end
if ~exist('ODsize','var')
    ODsize = 0.863; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex.
end

FOV = [range(x(1,:)),range(y(:,1))];
sz = size(x);
pixpermm = sz/FOV;
% Rojer and Schwartz' method of bandpassing random noise:

% Rojer, A.S. and E.L. Schwartz, Cat and monkey cortical columnar patterns
%modeled by bandpass-filtered 2D white noise. Biol Cybern, 1990. 62(5): p. 381-91.

%Make random noise: complex numbers where the angle is the orientation
z = exp(sqrt(-1)*rand(sz)*pi*2);

%Make a gabor filter
filtSz = 3; % 3mm
%sig = 5; %  .8 mm
freq = 1/ODsize; %cycles/mm (Try zero for big columns)
filtPix = ceil(filtSz*pixpermm);
[x,y] = meshgrid(linspace(-filtSz/2,filtSz/2,filtPix),linspace(-filtSz/2,filtSz/2,filtPix));
r = sqrt(x.^2+y.^2);
filt = exp(-r.^2/sig.^2).*cos(2*pi*freq*r);  %Gabor

%Convolve z with the filter
w = conv2(z,filt,'same');
pinwheel = angle(w);


[WX,WY] = gradient(w);

Gx = angle(WX);
OD = normcdf(Gx*sig);

 
