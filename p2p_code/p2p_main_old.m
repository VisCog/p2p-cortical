function [c,trl] = p2p_main()
c = define_cortex();
v = define_visualmap();
c = generate_corticalmap(c, v);

c = define_electrodes(c, v);
c = c2v(c);
c = find_rf_in_c(c);
c = gridv2c(c,v); % move visual grid to cortex

plotcortgrid(64 * (c.ORmap+pi)/(pi*2), c, 'Orientation pinwheels', hsv(64), 1, 'subplot(1, 1, 1)');
plotcortgrid(64 * c.ODmap, c, 'Ocular dominance columns', gray(64), 2, 'subplot(1, 1, 1)');
plotcortgrid(64*c.rf.size/6, c, 'Receptive field size', hot(64), 3, 'subplot(1, 1, 1)');

c = generate_ef(c);
plotcortgrid(c.e.ef * 64, c, 'electric field', gray(64), 4, 'subplot(1,1,1)');

tp = define_temporalparameters();
trls = define_trials(tp);

v = generate_rfmap(c, v);
v = generate_phosphene(c, v, tp, trls);

plotretgrid(v.e.phos_noRF*64, c,v, 'phosphene no RF', gray(64), 10, 'subplot(1,3,1)');
plotretgrid(v.e.phos(:, :, 1)*64, c, v,'phosphene', gray(64), 10, 'subplot(1,3,2)');
plotretgrid(v.e.phos(:, :, 2)*64, c, v, 'phosphene', gray(64), 10, 'subplot(1,3,3)');

end
        
% definitions
function v = define_visualmap()
v.e.pos = [26, 30]; % position of the electode in visual co-ordinates
v.retinaSize = [70,70]; %  [height, width diameter in degrees]
v.retinaCenter = [0,0];
v.pixperdeg = 10;
v.x = linspace(.5/v.pixperdeg,v.retinaSize(2)-.5/v.pixperdeg,v.retinaSize(2)*v.pixperdeg)+v.retinaCenter(1) - v.retinaSize(2)/2; 
v.y = linspace(.5/v.pixperdeg,v.retinaSize(1)-.5/v.pixperdeg,v.retinaSize(1)*v.pixperdeg)+v.retinaCenter(2) - v.retinaSize(1)/2;
[v.X,v.Y] = meshgrid(v.x, v.y);

%Make the grid in retinal coordinates
v.angList = -90:20:90;
v.radList = [1 2 3 5 8 13 21 34];
v.n = 201;
v.zAng = linspace(0,max(v.radList),v.n)'*exp(sqrt(-1)*v.angList*pi/180);
v.zRad = v.radList'*exp(sqrt(-1)*linspace(-90,90,v.n)*pi/180);
end
function c = define_cortex()
% cortical magnification, typical log z transformation parameters (Based on Duncan and Boynton)
c.k = 20; %scale
c.a = 0.5; %fovea expansion

% receptive field parameters
c.ar = .1;  % aspect ratio of V1 RFs
c.sig = .5; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex. sig determines the distribution of OD values. Default is 5.  The larger sig, the more the distribution tends toward 0 and 1.
c.ODsize = 0.863; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
c.thr = 0.014; % the level of cortical response that is considered above perceptual threshold
c.filtSz = 3; % 3mm creates the initial OD and orientation maps

% define the size and resolution of cortical and visual space parameters
c.gridColor = [1,1,0];
c.cortexSize = [80,100];  % %[height, width] Size of cortical maps (mm)
c.cortexCenter = [30,0]; % center of electrode array (mm on cortex)
c.pixpermm = 5; % choose the resolution to sample in mm.
end
function c = define_electrodes(c,v)
% either needs an electrode position in cortical co-ordinates
% or needs to take in the position of the electrode in visual co-ordinates
c.e.size = 1150/1000;
c.e.hemi = 'rh'; 
if ~isfield(c.e, 'pos')
    c.e.pos = vpos2cpos(v.e.pos, c);
end
% current spread parameters, based on Ahuja
c.e.afac=1.69; % parameters for current spread
c.e.ct=14;
end

% generate functions
function c = generate_corticalmap(c, v)
% define cortex meshgrid
c.x = linspace(.5/c.pixpermm,c.cortexSize(2)-.5/c.pixpermm,c.cortexSize(2)*c.pixpermm) + c.cortexCenter(1) - c.cortexSize(2)/2;
c.y = linspace(.5/c.pixpermm,c.cortexSize(1)-.5/c.pixpermm,c.cortexSize(1)*c.pixpermm) + c.cortexCenter(2) - c.cortexSize(1)/2;
[c.X,c.Y] = meshgrid(c.x,c.y);

% Make the orientation and OD maps
c = generate_ORmapODmap(c);
end
function c = find_rf_in_c(c)
% finds the angle, eccentricity and size of rfs for every point in cortex
c.rf.size = max(.1652 .* c.c2v.ECC, 0.9);
end
function [c] = generate_ORmapODmap(c)
% [pinwheel,OD] = makePinwheelODMaps(x,y,sig)

FOV = [range(c.X(1,:)),range(c.Y(:,1))];
sz = size(c.X);
pixpermm = sz/FOV;
% Rojer and Schwartz' method of bandpassing random noise:

% Rojer, A.S. and E.L. Schwartz, Cat and monkey cortical columnar patterns
%modeled by bandpass-filtered 2D white noise. Biol Cybern, 1990. 62(5): c. 381-91.

%Make random noise: complex numbers where the angle is the orientation
Z = exp(sqrt(-1)*rand(sz)*pi*2);

% filter the noise to create initial columns
freq = 1/c.ODsize; %cycles/mm (Try zero for big columns)
filtPix = ceil(c.filtSz*pixpermm);
[X,Y] = meshgrid(linspace(-c.filtSz/2,c.filtSz/2,filtPix),linspace(-c.filtSz/2,c.filtSz/2,filtPix));
R = sqrt(X.^2+Y.^2);
FILT = exp(-R.^2/c.sig.^2).*cos(2*pi*freq*R);  %Gabor

%Convolve z with the filter
W = conv2(Z,FILT,'same');
c.ORmap = angle(W);
WX = gradient(W);
Gx = angle(WX);
c.ODmap = normcdf(Gx*c.sig);
end
function [c] = generate_ef(c)
if ~isfield(c.e, 'boxsize') % how far from each electrode to simulate
    c.e.boxsize = 3;
end
R=sqrt((c.X-c.e.pos(1)).^2+(c.Y-c.e.pos(2)).^2);
pt_ef = ones(size(c.X));
pt_ef(R>c.e.size)=2/pi*(asin(c.e.size./R(R>c.e.size)));
ss=c.x(1,2)-c.x(1,1);
[X,Y]= meshgrid(min(c.e.size*c.e.boxsize, -2):ss:max(ss:c.e.size*c.e.boxsize, 2));
d=sqrt(X.^2+Y.^2);
C_att=c.e.ct./(c.e.ct+d.^c.e.afac);
ef = conv2(pt_ef, C_att, 'same');
ef=ef(:)./max(ef(:));
c.e.ef = reshape(ef,size(c.X));
end
function v  = generate_rfmap(c, v)
cosMap = cos(c.ORmap);
sinMap = sin(c.ORmap);
v.e.phos_noRF = zeros(size(v.X)); % percept based on electric field
v.e.phos = zeros([size(v.X), 2]); % percept that includes a cortical model

for pixNum = 1:length(c.X(:))
    if mod(pixNum, 1000)==0
        disp([num2str(round((100*pixNum)/length(c.X(:)))),  '% complete' ]);
    end
    if c.e.ef(pixNum) > 0.5
        % created oriented mesh in visual space 
        xrot = (v.X-c.c2v.X(pixNum)) * cosMap(pixNum) +  ...
            (v.Y-c.c2v.Y(pixNum)) * sinMap(pixNum);
        yrot = -(v.X-c.c2v.X(pixNum)) * sinMap(pixNum) + ...
            (v.Y-c.c2v.Y(pixNum)) * cosMap(pixNum);
        
        % scoreboard version
        v.e.phos_noRF = v.e.phos_noRF + (c.e.ef(pixNum) * exp(-( xrot.^2/(0.1) + yrot.^2/(.1))));
        
        % Cortical Model version
        v.e.phos(:, :, 1)  =   v.e.phos(:, :, 1)  + (c.ODmap(pixNum) * c.e.ef(pixNum) * exp(-( xrot.^2/(2 * c.rf.size(pixNum) * c.ar) + yrot.^2/(2 * c.rf.size(pixNum)))));
        v.e.phos(:, :, 2)  =   v.e.phos(:, :, 2)  + ( (1-c.ODmap(pixNum)) * c.e.ef(pixNum) * ...
            exp(-( xrot.^2/(2 * c.rf.size(pixNum) * c.ar) + yrot.^2/(2 * c.rf.size(pixNum)))));      
    end
end
end
function v = generate_phosphene(c,v, tp, trls)
        % calculate the neural response over time for each trial
        for t=1:length(trls)
           v.trls(t).phos = max(p2p_finite_element(tp, trls(t).pt)); % the scaling due to current integration over time
        end
end
%         idx=find(trls(t).pt>0);
%         if ii==5; idx=setdiff(idx,3:5); end % drawings not obtained
%         
%         poly_area = D(ii).poly_area.val(idx);
%         poly_area(isnan(poly_area))=0;
%             
%         if ~c.smallflag
%             p.thr = 50;  p=UWfit('findDrawingThreshold', p, {'thr'}, Dsim(ii).EF.normalizedphosphene(:), Dsim(ii).CI.val(idx), ...
%                 (poly_area * c.pixperdeg^2));
%             Dsim(ii).EF.drawingthreshold = p.thr;
%             disp(['EF drawing threshold = ', num2str(Dsim(ii).EF.drawingthreshold)]);
%             
%             p.thr = 50; p=UWfit('findDrawingThreshold', p, {'thr'}, Dsim(ii).CM.normalizedphosphene(:), Dsim(ii).CI.val(idx), ...
%                 (poly_area * c.pixperdeg^2));
%             Dsim(ii).CM.drawingthreshold = p.thr;
%             disp(['CM drawing threshold = ', num2str(Dsim(ii).CM.drawingthreshold)]);
%         end
%         for cc=1:length(idx)
%             Dsim(ii).EF.poly_area.val(idx(cc)) = (1/(c.pixperdeg^2)) * ...
%                 length( find((Dsim(ii).EF.normalizedphosphene(:) * Dsim(ii).CI.val(idx(cc)) ) ...
%                 >  Dsim(ii).EF.drawingthreshold) );
%             Dsim(ii).CM.poly_area.val(idx(cc)) = (1/(c.pixperdeg^2)) * ...
%                 length( find((Dsim(ii).CM.normalizedphosphene(:) *  Dsim(ii).CI.val(idx(cc)) ) ...
%                 >  Dsim(ii).CM.drawingthreshold) );
%         end
% 
% end
% transformation functions
function c = gridv2c(c,v)
c.wAng = v2c(c, v.zAng);
c.wRad = v2c(c, v.zRad);
c.cropPix = ~(c.c2v.ANG<max(v.angList) & c.c2v.ANG>min(v.angList) & c.c2v.ECC<=max(v.radList));
end
function cpos = vpos2cpos(vpos, c)
% converts positions in visual coordinates into cortical coordinates
tmp =vpos(1).*exp(sqrt(-1)*vpos(2)*pi/180); % turn them into complex numbers, in degrees
tmp2=v2c(c, tmp); % turn into mm, imaginary
cpos = [real(tmp2), imag(tmp2)]; % turn into mm real
end
function v = v2c(c, z)
v = c.k*log(z + c.a);
end
function c = c2v(c)
tmp = c.X+sqrt(-1)*c.Y;
%v = map(z,p)
z = exp(tmp/c.k)-c.a;
c.c2v.X = real(z);
c.c2v.Y = imag(z);
c.c2v.ECC = abs(tmp);
c.c2v.ANG = angle(tmp) * 180/pi;
end

% time functions
function tp = define_temporalparameters()
tp.dt = 1/1000; % time sampling in ms
tp.tau1 =.2/1000; %.42/1000 for retina. Is this right???
tp.tau2_ca = 45.25/1000;  %38-57
tp.tau3 =  26.25/1000; % 24-33ex
tp.e = 8.73; %2.25;  %2-3 for threshold or 8-10 for suprathreshold

% leak out of charge accumulation
tp.flag_cl=0; % 1 if you want to charge to leak back out of the system
tp.tau2_cl = tp.tau2_ca * 8; % used for the conv model, fe uses p.tau2_ca

% nonlinearity parameters
tp.slope=.3; % larger number = shallower slope
tp.asymptote=14;
tp.shift=47; % shifts curve along x-axis
end
function pt = generate_pt(trl, tp)
if isnan(trl.freq) % if not using a temporal model at all
    pt = 1;
else
    p = cat(1, zeros(round(trl.lag./tp.dt), 1), ...
        trl.order * ones(round(trl.pw./tp.dt), 1), ...
        zeros(round(trl.ip./tp.dt), 1), ...
        -1 * trl.order * ones(round(trl.pw./tp.dt), 1));
    
    if trl.freq == -1 % single pulse
        tmp = zeros(trl.dur./tp.dt, 1);
        tmp(1:length(p)) = p;
        pt = tmp;
    else % not a single pulse
        tmp = zeros(round((1000/trl.freq)./tp.dt), 1);
        tmp(1:length(p)) = p;
        pt = repmat(tmp, floor(length(trl.t)./length(tmp)), 1);
    end
end
end
function trl = define_trials(tp)
trl.name = 'generic';
trl.dur = 150; % duration in ms
trl.t = 0:tp.dt:trl.dur-tp.dt;
trl.ip = 0; % interphase delay
trl.lag = 50/1000; % delay before the pulse train begins in ms
trl.pw = 200/1000; % pulse width in ms
trl.order = 1; % 1 = cathodic first, -1  = anodic first
trl.freq = 60; % -1 for a single pulse, NaN if not using a temporal model
trl.amp = 50; % current amplitude in microAmps
trl.pt = generate_pt(trl, tp);
end
function  out = p2p_finite_element(tp, tsform )
% [ out ] = p2p_finite_element(tp, tsform)
%
% Implements accumulation of current over time using a very simple finite element method
% written GMB 11/10/2017
% adapted for cortex IF 3/2/2018
% Finite difference method:
% Initial conditions

tmp.chargeacc = 0;
tmp.ca = 0;
tmp.cl = 0;
tmp.R1 = 0;


tmp.R3norm = 0;
tmp.R4a =  zeros(1,4);
tmp.R4ax =  zeros(1,4);

for i=1:length(tsform)-1
    % R1
    tmp.R1= tmp.R1 + tp.dt * (tsform(i)-tmp.R1)/tp.tau1;
    tmp.R4a(:,1) = max(tmp.R1, 0);
    for j=1:3
        tmp.R4a(:,j+1) = tmp.R4a(:,j+1) + tp.dt*(tmp.R4a(:,j) - tmp.R4a(:,j+1))/tp.tau3;
    end
    tmp.R4(i) = tmp.R4a(:,4);
end

out=tmp.R4;
end

% plotting functions
function plotcortgrid(img, c, titlestr, cmap,figNum, spstr)

% if isfield(c,'cropPix')
%     img(c.cropPix) = NaN;
%     img= img+2;
%     cmap = [0,0,0;cmap];
% end

fH=figure(figNum);set(fH, 'Name', titlestr);
eval(spstr); colormap(cmap);
image(c.x, c.y, flipud(img)); hold on
xlabel('mm'); ylabel('mm')
set(gca,'YDir','normal');
plot(c.wAng, '-', 'Color', c.gridColor);
plot(c.wRad', '-', 'Color', c.gridColor);

axis equal;  axis tight
set(gca,'XLim',[min(c.x(:)),max(c.x(:))]);
set(gca,'YLim',[min(c.y(:)),max(c.y(:))]);
end
function plotretgrid(img, c,v, titlestr, cmap, figNum, spstr)

fH=figure(figNum);set(fH, 'Name', titlestr);
eval(spstr);
 if strcmp(c.e.hemi, 'lh')
     image(v.x, v.y, img); hold on
 else
    image(v.x, v.y, fliplr(img)); hold on
 end

colormap(cmap);
set(gca,'YDir','normal');

plot(v.zAng,'-','Color',c.gridColor);
plot(v.zRad','-','Color',c.gridColor);
plot(-v.zAng,'-','Color',c.gridColor);
plot(-v.zRad','-','Color',c.gridColor);

axis equal;  axis tight
xlabel('degrees'); ylabel('degrees')
set(gca,'XLim',[min(v.x(:)),max(v.x(:))]);
set(gca,'YLim',[min(v.y(:)),max(v.y(:))]);

end



         

%% NOTES
% USING CORTICAL POSITION
% Winawer site one is LH, all the rest are RH
% original values based on the Benson template extracted from Winawer, are converted into x y co-ordinates using:
% Z=eData(s,1).*exp(sqrt(-1)*eData(s,2)*pi/180); xe=real(Z); ye=imag(Z);
% orig values in visual space co-ordinates were ecc, angle =
% ea_ret = [7.8774	8.8141; 18.2743	23.5752; 2.4848	48.3118 ; 7.8671	166.642 ;  3.337	49.077]; % pulled from Win paper, in degrees
% % convert so all projecting to the same hemisphere
% ea_ret = [7.8774	8.8141; 18.2743	23.5752; 2.4848	48.3118 ; 7.8671	180-166.642 ;  3.337	49.077]; % ecc, angle, pulled from Win paper, in degrees

% USING PHOSPHENE LOCATION
%xyloc = [25 9; -7.65 -1.85; -4.05 3.14; -1.35 1.35; -.75 .5]; % pulled using ginput from fig 2
% ea_ret = [ 26.6 19.8; 9  -166.4; 5.12 142.2; 1.9 135; 0.90 146.3]; % ecc, angle
% ea_ret = [ 26.6 19.8; 9  -166.4; 5.12 142.2; 1.9 135; 1 146.3]; % ecc, angle

%ea_ret = [ 26.6 19.8; 9  -(180-166.4); 5.12 (180-142.2); 1.9 (180-135); 1 (180-146.3)]; % ecc, angle, flipped 2:5 to right visual