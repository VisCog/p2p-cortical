function img = makePinwheelMap(voxSize,x,y)

%initial parameters:
FOV = [range(x(1,:)),range(y(:,1))];
sz = size(x);
pixpermm = sz/FOV;
%Rojer and Schwartz' method of bandpassing random noise:

%Rojer, A.S. and E.L. Schwartz, Cat and monkey cortical columnar patterns
%modeled by bandpass-filtered 2D white noise. Biol Cybern, 1990. 62(5): p. 381-91.

%Make random noise: complex numbers where the angle is the orientation
z = exp(sqrt(-1)*rand(sz)*pi*2);

%Make a gabor filter
filtSz = 3; % 3mm
sig = .8; %  .8 mm
freq = 1; %cycles/mm (Try zero for big columns)
filtPix = ceil(filtSz*pixpermm);
[x,y] = meshgrid(linspace(-filtSz/2,filtSz/2,filtPix),linspace(-filtSz/2,filtSz/2,filtPix));
r = sqrt(x.^2+y.^2);
filt = exp(-r.^2/sig.^2).*cos(2*pi*freq*r);  %Gabor

%Convolve z with the filter
img = angle(conv2(z,filt,'same'));

