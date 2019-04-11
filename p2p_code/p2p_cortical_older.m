function D=p2p_cortical()
% Reproduce Figure 4 for the following paper:
%
%   Winawer and Parvizi (2016). Linking Electrical Stimulation of Human
%   Primary Visual Cortex, Size of Affected Cortical Area, Neuronal
%   Responses, and Subjective Experience Neuron. 92(6):1213?1219
%   http://dx.doi.org/10.1016/j.neuron.2016.11.008

close all

% sites 1 - 5 are from Winawer
sites = 1;

addpath(pwd);
addpath('utilities')
opts = getPlotOpts();  % Plot options

c = getcorticalparams();

D = getData(sites);    % Data per trials

Dsim=simulateData(D, c); % generate simulated data to match
D = deriveData(D, opts); % Some additional derived parameters

opts.titlestr = 'winawer, 4a';
figure4a(D, opts);

opts.titlestr = 'phosphene 4a'; % cortical simulation
Dtmp=D; 
for ii = 1:length(D)
    Dtmp(ii).poly_area.val = Dsim(ii).A(1).chargedensity.val;
end
figure4a(Dtmp,  opts);

opts.titlestr = 'EF 4a'; % electric field
for ii = 1:length(D)
    Dtmp(ii).poly_area.val = Dsim(ii).A(2).chargedensity.val;
end
figure4a(Dtmp,  opts);
return


% Figure 4b
figure4b(D, opts);

% Figure 4c
figure4c(D, opts);

% Figure 4d
figure4d(D, opts);

% Figure 4e
figure4e(D, opts);

end

function Dsim = simulateData(D,c)

titlestr={'cortical model' 'electrical model'};
subplotstr=['subplot(2, 3, '];

%t = gettemporalmodelparams();

plotcortgrid(64 * (c.orientationMap+pi)/(pi*2), c, 'Orientation pinwheels', hsv(64), 10, 'subplot(1, 1, 1)');
plotcortgrid(64 * c.odMap, c, 'Ocular dominance columns', gray(64), 11, 'subplot( 1,1, 1)');

for ii=1:length(D) % for each site
    disp(['simulating electric field for site ', num2str(ii)])
    Dsim(ii).ERF = makeElectrodeERF(c, c.xy_cort(ii, :), c.e_size(ii));
    
    plotcortgrid(reshape(Dsim(ii).ERF * 64,size(c.x)), c, 'electric field', gray(64), 12, [subplotstr, num2str(ii), ')']);
    
    disp(['simulating phosphene for site ', num2str(ii)])
    
    for jj = 1:2 % 4 images, 1 for each eye, with and without neural modeling
        percept(jj).img = zeros(size(c.xr)); % left eye
    end
    %map from cortex to retinal coordinates
    vr = mapinv(c, c.x + sqrt(-1)*c.y);   
    RFsize = min(.9 + .1652*abs(vr), 0.9);  %based figure in Simoncelli and Freeman, Metamers ventral stream, NN 2011
            % assume inflexion at [6, .9] [20 3.35], 45 7.48], slope is
            % 0.1652
    RForientation = c.orientationMap;
    
    for pixNum = 1:length(c.x(:))
        if mod(pixNum, 1000)==0
            disp([num2str(round((100*pixNum)/length(c.x(:)))),  '% complete' ]);
        end
        if Dsim(ii).ERF(pixNum) > c.thr / 4
        
      % Oriented gaussian
            xrot = (c.xr-real(vr(pixNum)))*cos(RForientation(pixNum))+(c.yr-imag(vr(pixNum))) * sin(RForientation(pixNum));
            yrot = -(c.xr-real(vr(pixNum)))*sin(RForientation(pixNum)) + (c.yr-imag(vr(pixNum))) * cos(RForientation(pixNum));
            G = exp(-( xrot.^2/(2*RFsize(pixNum)*c.ar) + yrot.^2/(2*RFsize(pixNum))));
            percept(1).img =  percept(1).img  + Dsim(ii).ERF(pixNum) * G;   
            % percept(3).img =  percept(2).img  + (1-S.odMap(pixNum))*ERF(pixNum)*G; % right eye
            
            % scoreboard version
            G_SB = exp(-( xrot.^2/(0.01) + yrot.^2/(.01)));
            percept(2).img =  percept(2).img + Dsim(ii).ERF(pixNum) * G_SB;    
        end
    end
    for pp=1:length(percept)
        percept(pp).img = percept(pp).img./max(percept(1).img(:));
        plotretgrid(percept(pp).img * 64, c,['phosphene ', titlestr{pp}], gray(64), 12+pp, [subplotstr, num2str(ii), ')']);
        
        for cc=1:length(D(ii).condition.val)
            Dsim(ii).A(pp).chargedensity.val(cc) =         (1/c.pixperdeg^2) * length( find((percept(pp).img(:) * D(ii).chargedensity.val(cc))> c.thr) );
            Dsim(ii).A(pp).chargedensityPerPulse.val(cc) = (1/c.pixperdeg^2) * length( find((percept(pp).img(:) * D(ii).chargedensityPerPulse.val(cc)) > c.thr) );
            Dsim(ii).A(pp).totalCharge.val(cc) =           (1/c.pixperdeg^2) * length( find((percept(pp).img(:) * D(ii).totalCharge.val(cc)) > c.thr) );
            Dsim(ii).A(pp).chargedensityPerTime.val(cc) =  (1/c.pixperdeg^2) * length( find((percept(pp).img(:) * D(ii).chargedensityPerTime.val(cc)) > c.thr) );
        end
    end
end


end % end of nesting

function ERF = makeElectrodeERF(c, eloc, e_size)

disp('making electrode field');
R=sqrt((c.x-eloc(1)).^2+(c.y-eloc(2)).^2);
ERF = ones(size(c.x));
ERF(R>e_size)=2/pi*(asin(e_size./R(R>e_size)));
ERF=ERF(:);

end

function plotretgrid(img,c, titlestr, cmap,figNum, spstr)

fH=figure(figNum);set(fH, 'Name', titlestr);
eval(spstr);
image(c.xrax, c.yrax, img); hold on
 colormap(cmap);
set(gca,'YDir','normal');

plot(c.zAng,'-','Color',c.gridColor);
plot(c.zRad','-','Color',c.gridColor);
axis equal;  axis tight
xlabel('degrees'); ylabel('degrees')
set(gca,'XLim',[min(c.xrax(:)),max(c.xrax(:))]);
set(gca,'YLim',[min(c.yrax(:)),max(c.yrax(:))]);

end

function plotcortgrid(img, c, titlestr, cmap,figNum, spstr)

fH=figure(figNum);set(fH, 'Name', titlestr);
eval(spstr); colormap(cmap);
image(c.xax, c.yax, img); hold on
xlabel('mm'); ylabel('mm')
set(gca,'YDir','normal');
plot(c.wAng, '-', 'Color', c.gridColor);
plot(c.wRad', '-', 'Color', c.gridColor);

axis equal;  axis tight
set(gca,'XLim',[min(c.xax(:)),max(c.xax(:))]);
set(gca,'YLim',[min(c.yax(:)),max(c.yax(:))]);
end

function opts = getPlotOpts

fieldsToPlot = {'chargedensity' 'chargedensityPerPulse' 'totalCharge' 'chargedensityPerTime'};
opts.fieldToPlot = fieldsToPlot{3};

% Plot cortical area using individual maps ('area') or standard cm function
% ('cm_area')
opts.area = 'area';

% Plot colors
[opts.colors, opts.sites] = getColors;

% Legend text
opts.leg_txt = cellstr(num2str(opts.sites'));

% Axis scale and range
opts.yscale = 'log';
opts.xscale = 'log'; % 'linear'
opts.yl     = 10.^[-1 3]; % 'linear'

switch opts.xscale
    case 'log', opts.xl = 10.^[1.5 4.5];
    case 'linear', opts.xl = [0 1.2e4];
end

end

function D = getData(sites)

pth = fullfile(ebsRootPath, 'data', 'ebs');
for ii = sites
    fname = sprintf('trial_data_site%d', ii);
    D(ii) = load(fullfile(pth, fname));
end

%[~,~,~,units] = getConditionsFromFile('jt', pth.data);

end

function D = deriveData(D, opts)

for ii = 1:length(D)
    
    % Which stimulation parameter to plot?
    D(ii).stimulation_data.val = D(ii).(opts.fieldToPlot).val;
    
    % Which measure of surface area to plot? (derived from indivual
    % retinoptic map or from standard CM function)
    D(ii).surface_area.val = D(ii).(opts.area);
    
    % Derive number of pulses from frequency and duration
    D(ii).num_pulses.val            = D(ii).frequency.val .* D(ii).duration.val;
    
    % Derive one rating pre trial from separate ratings for color, motion,
    % brightness
    D(ii).subjective_rating.val     = nanmedian([...
        D(ii).motion.val;...
        D(ii).color.val; ...
        D(ii).brightness.val...
        ]);
    
end

% If we are plotting the surface area using indivual retinotopic maps, then
% we interpolate the area based on the standard CMF. This is because for
% small phosphenes, there may be no voxel whose center is insider the
% phosphene, yet we know the surface area cannot be 0.
for ii = 1:numel(D)
    idx = isfinite(D(ii).which_drawing.val);
    lm = fitlm(D(ii).cm_area.val(idx), D(ii).area.val(idx), 'Intercept', false);
    D(ii).area.val(idx) = lm.predict(D(ii).cm_area.val(idx)');
    D(ii).surface_area.val = D(ii).(opts.area).val';
end

end

function c = getcorticalparams()

experiment = 'winawer'
if strcmp(experiment, 'winawer')
% original values based on the Benson template extracted from Winawer
% didn't match the percept very well, so used values based on the center of
% mass  visual space co-ordinates, values are ecc, angle
ea_ret = [7.8774	8.8141; 18.2743	23.5752; 2.4848	48.3118 ; 7.8671	166.642 ;  3.337	49.077]; % pulled from Win paper, in degrees
% convert so all projecting to the same hemisphere
ea_ret = [7.8774	8.8141; 18.2743	23.5752; 2.4848	48.3118 ; 7.8671	180-166.642 ;  3.337	49.077]; % pulled from Win paper, in degrees
c.e_size=[1150/1000; 510/1000;  1150/1000;  1150/1000;  1150/1000]; %radius exposed electrode, in mm
elseif strcmp(experiment, 'hires')
    ea_ret = [ .2 0 ; .5 30; 1 50; 2 160; 4 50 ; 8 30; 16 0];
end
    
% visual to cortical mapping
% typical log z transformation parameters (Based on Duncan and Boynton)
c.k = 20; %scale
c.a = 0.5; %fovea expansion
c.ar = .1;  % aspect ratio of V1 RFs
c.gridColor = [1,1,0];

c.thr = 0.1; % the level of cortical response that is considered above perceptual threshold
% cortical and visual space parameters
c.cortexSize = [90,120];  % %[height, width] Size of cortical maps (mm)
c.cortexCenter = [0,0]; % center of electrode array (mm on cortex)

c2r = 3; % each degree is represented by about ~280 microns in the retina
c.retinaSize = [40,40]; %  [height, width in degrees]
c.retinaCenter = [.75,0];

% define the resolution of the cortical and visual maps
pixperelectrode = 5; % how fine resolution to sample, scaled to smallest electrode
c.pixpermm = pixperelectrode /min(c.e_size); 
c.pixperdeg = c.pixpermm*c2r;

%Make the grid in retinal coordinates
c.angList = [-90:20:90];
c.radList = [ 2.5, 5, 10, 15, 20, 25];
c.n = 201;

disp('generating visual/cortical template');
% define cortex meshgrid

c.xax = linspace(.5/c.pixpermm,c.cortexSize(2)-.5/c.pixpermm,c.cortexSize(2)*c.pixpermm) + c.cortexCenter(1) - c.cortexSize(2)/2;
c.yax = linspace(.5/c.pixpermm,c.cortexSize(1)-.5/c.pixpermm,c.cortexSize(1)*c.pixpermm) + c.cortexCenter(2) - c.cortexSize(1)/2;
[c.x,c.y] = meshgrid(c.xax,c.yax);

% Make the orientation and OD maps
[c.orientationMap,c.odMap] = makePinwheelODMaps(c.x,c.y, 0.863);

c.xrax = linspace(.5/c.pixperdeg,c.retinaSize(2)-.5/c.pixperdeg,c.retinaSize(2)*c.pixperdeg)+c.retinaCenter(1) - c.retinaSize(2)/2;
c.yrax = linspace(.5/c.pixperdeg,c.retinaSize(1)-.5/c.pixperdeg,c.retinaSize(1)*c.pixperdeg)+c.retinaCenter(2) - c.retinaSize(1)/2;
[c.xr,c.yr] = meshgrid(c.xrax,c.yrax);

ea_retZ =ea_ret(:,1).*exp(sqrt(-1)*ea_ret(:,2)*pi/180); % turn them into complex numbers, in degrees
ea_cortZ=map(c, ea_retZ); % turn into mm, imaginary
c.xy_cort = [real(ea_cortZ), imag(ea_cortZ)]; % turn into mm real

%% RF to V1 mapping
c.zAng = linspace(0,max(c.radList),c.n)'*exp(sqrt(-1)*c.angList*pi/180);
c.zRad = c.radList'*exp(sqrt(-1)*linspace(-90,90,c.n)*pi/180);

% project the grid to cortical coordinates (mapinv goes the other way)
c.wAng = map(c, c.zAng);
c.wRad = map(c, c.zRad);

end

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
%modeled by bandpass-filtered 2D white noise. Biol Cybern, 1990. 62(5): c. 381-91.

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

end

function v = map(c, z)
v = c.k*log(z + c.a);
end

function z = mapinv(c, v)
%v = map(z,p)
z = exp(v/c.k)-c.a;
end

function t = gettemporalmodelparams()
t.tau1 = .42/1000;  %.24 - .65
t.tau2_ca = 45.25/1000;  %38-57
t.tau3 =  26.25/1000; % 24-33ex
t.e = 8.73; %2.25;  %2-3 for threshold or 8-10 for suprathreshold\
% leak out of charge accumulation
t.flag_cl = 0; % 1 if you want to charge to leak back out of the system
t.tau2_cl = t.tau2_ca * 8; % used for the conv model, fe uses p.tau2_ca

t.slope = .3; % larger number = shallower slope
t.asymptote = 14;
t.shift = 47; % shifts curve along x-axis
t.tsample = 0.25/1000; % time sampling
end

function fH = figure4a(D, opts)
fH = figure;

set(fH, 'Color', 'w' ,'name', opts.titlestr, 'NumberTitle', 'off');

set(gca, 'FontSize', 10); hold on

fit_type = 'power';

x_all = [];
for ii = 1:length(D)
    idx = isfinite(D(ii).which_drawing.val);
    x_all = [x_all D(ii).stimulation_data.val(idx)];
end

x_all = sort(unique(x_all));

for ii = 1:length(D)
    idx = isfinite(D(ii).which_drawing.val);
    switch fit_type
        case 'linear'
            b = regress(D(ii).poly_area.val(idx)', ...
                D(ii).stimulation_data.val(idx)');
            pred_y = x_all * b;
            
            
        case 'power'
            [f, gof] = fit(D(ii).stimulation_data.val(idx)', ...
                D(ii).poly_area.val(idx)','b*x^m', ...
                'StartPoint',[mean(D(ii).poly_area.val(idx)) ...
                / mean(D(ii).stimulation_data.val(idx)) 1]);
            pred_y = f(x_all);
            disp(f)
            disp(gof)
            
    end
    plot(D(ii).stimulation_data.val(idx), D(ii).poly_area.val(idx),...
        'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
    plot(x_all,  pred_y, 'Color', opts.colors(ii,:), 'LineWidth', 2)
end

set(gca, 'YScale', opts.yscale, 'XScale', opts.xscale, 'XTick', 10.^[0 1 2 3], 'XLim', 10.^[0 3])
if strcmp(opts.yscale, 'log'), set(gca, 'YLim', 10.^[-3 3]); end

xlabel('Charge Deposited per Trial (µC)')
ylabel('Phosphene size (deg^2)')

end

function fH = figure4b(D, opts)

fH = figure;
set(fH, 'Color', 'w' ,'name', 'Winawer and Parvizi Figure 4B', 'NumberTitle', 'off');

set(gca, 'FontSize', 30); hold on

x_all = []; for ii = 1:length(D); x_all = [x_all D(ii).stimulation_data.val]; end
x_all = unique(x_all);

for ii = 1:length(D)
    plot(D(ii).stimulation_data.val, D(ii).surface_area.val,...
        'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
end

for ii = 1:length(D)
    
    idx = isfinite(D(ii).surface_area.val);
    x = D(ii).stimulation_data.val(idx)';
    y = D(ii).surface_area.val(idx);
    
    [f, gof] = fit(x,y,'b*x^m', 'StartPoint',[mean(x) / mean(y) 1]);
    pred_y = f(x_all);
    
    disp(f)
    disp(gof)
    
    plot(D(ii).stimulation_data.val, D(ii).surface_area.val, ...
        'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
    plot(x_all,  pred_y, 'Color', opts.colors(ii,:), 'LineWidth', 2)
end

xlabel('Charge Deposited Per Trial (µC)');
ylabel('Cortical Area (mm^2)');

set(gca, 'YScale', opts.yscale, 'XScale', opts.xscale, 'XLim', 10.^[0 3], 'XTick', 10.^[0 1 2 3])
if strcmp(opts.yscale, 'log'), set(gca, 'YLim', 10.^[-3 3]); end
% if strcmp(xscale, 'log'), set(gca, 'XLim', 10.^[1.5 4.5]); end

xl = get(gca, 'XLim');
plot(xl, 1.13 * [1 1], 'k--')
plot(xl, 4.15 * [1 1], 'k--')

end

function fH = figure4c(D, opts)
fH = figure;
set(gcf, 'Color', 'w' ,'name', 'Winawer and Parvizi Figure 4C');

% all channels, one plot
hold on, set(gca, 'FontSize', 30)
fit_type = 'power';

for ii = 1:numel(D)
    
    [x{ii}, inds] = sort(D(ii).poly_area.val');
    y{ii} = D(ii).surface_area.val(inds);
    idx = isfinite(x{ii});
    x{ii} = x{ii}(idx);
    y{ii} = y{ii}(idx);
    sz = D(ii).stimulation_data.val;
    sz = sz(inds) / max(sz)*200;
    
    
    plot(x{ii}, y{ii}, 'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
    %scatter(x{ii}, y{ii}, sz, 'ko','MarkerFaceColor', opts.colors(ii,:))
    
    switch fit_type
        case 'power'
            % power law fit
            [f, gof] = fit(x{ii},y{ii},'b*x^m', 'StartPoint', [mean(y{ii})/mean(x{ii}) 1]);
            xpred{ii} = [min(x{ii})/10; x{ii}; max(x{ii})*10];
            ypred{ii} = f(xpred{ii});
            disp(f)
            disp(gof)
        case 'linear'
            % linear fit
            [b{ii}, ~ ,~, ~, stats] = regress(y{ii}, x{ii});
            ypred{ii} = x{ii}* b{ii};
            r2(ii)=(stats(1));
    end
end
for ii = 1:length(D)
    plot(xpred{ii}, ypred{ii}, '-', 'Color', opts.colors(ii,:), 'LineWidth', 2);
end
xl =  10.^[-3 3];
yl = 10.^[0 3];
set(gca, 'YScale', opts.yscale,  'XScale', 'log', 'YLim', yl, 'XLim',xl, 'XTick', 10.^[-2 0 2])
ylabel('Cortical area (mm^2)'), xlabel('Phosphene size (deg^2)')
legend(opts.leg_txt, 'Location', 'Best')

end

function fH = figure4d(D, opts)

fH = figure; set(fH, 'Color', 'w' ,'name', 'Winawer and Parvizi Figure 4D');

x_all = [];
for ii = 1:length(D)
    idx = isfinite(D(ii).which_drawing.val);
    x_all = [x_all D(ii).stimulation_data.val(idx)];
end

x_all = sort(unique(x_all));


set(gca, 'FontSize', 30); hold on
% title('Subjective intensity rating')

fit_type = 'power';

% % if a trial has the same electrode and same condition number as the trial
% % before, then it must indicate a second drawing on the same trial (subject
% % drew two phosphenes for one stimulation and we only use the first)
% isrepeat =  [1  diff(electrode)] == 0 & [1 diff(condition)] == 0;

for ii = 1:length(D)
    
    if all(isnan(D(ii).subjective_rating.val))
    else
        
        idx = isfinite(D(ii).subjective_rating.val);
        switch fit_type
            case 'linear'
                b = regress(D(ii).subjective_rating.val(idx)', ...
                    D(ii).stimulation_data.val(idx)');
                pred_y = x_all * b;
                
                
            case 'power'
                [f, gof] = fit(D(ii).stimulation_data.val(idx)', ...
                    D(ii).subjective_rating.val(idx)','b*x^m', ...
                    'StartPoint',[mean(D(ii).subjective_rating.val(idx)) ...
                    / mean(D(ii).stimulation_data.val(idx)) 1]);
                pred_y = f(x_all);
                disp(f)
                disp(gof)
        end
        
        plot(D(ii).stimulation_data.val(idx), D(ii).subjective_rating.val(idx),...
            'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
        plot(x_all,  pred_y, 'Color', opts.colors(ii,:), 'LineWidth', 2)
        
    end
end

yscale = 'linear';
set(gca, 'YScale', yscale, 'XScale', opts.xscale, 'YLim', [0 11], ...
    'YTick', 0:2:10, 'XLim', 10.^[0 3], 'XTick', 10.^[0 1 2 3])
if strcmp(yscale, 'log'), set(gca, 'YLim', 10.^[0 1]); end

xlabel('Charge Deposited per Trial (µC)')
ylabel('Subjective rating')

end

function fH = figure4e(D, opts)

x1 = 'chargePerPulse';
x2 = 'frequency'; % 'num_pulses'
% dvs = {'surface_area_indiv' 'surface_area_cm' 'subjective_rating'}; % 'phosphene_area'
dvs = { 'surface_area' 'subjective_rating'};
x1all = []; x2all = [];
for ii = 1:length(D)
    idx = isfinite(D(ii).which_drawing.val);
    x1all = [x1all D(ii).(x1).val(idx)];
    x2all = [x2all D(ii).(x2).val(idx)];
end

fH = figure; pos = get(fH, 'Position'); pos(3:4) = [800 600];
set(gcf, 'Color', 'w', 'name', 'Winawer & Parvizi, Figure S4',...
    'NumberTitle', 'off', 'Position', pos);

for z = 1:length(dvs)
    dv = dvs{z};
    
    for ii = 1:length(D)
        
        if all(isnan(D(ii).(dv).val)), skip = true; else, skip = false; end
        
        disp(dv)
        subplot(length(D),length(dvs),(ii-1)*length(dvs)+z)
        set(gca, 'FontSize', 15); hold on;
        if ~skip
            idx = isfinite(D(ii).which_drawing.val);
            
            lm = fitlm([D(ii).(x1).val(idx)', D(ii).(x2).val(idx)'],...
                D(ii).(dv).val(idx)', 'linear',  'varnames', {x1, x2, dv});
            
            lm.plotEffects; xl = get(gca, 'XLim'); set(gca, 'XLim', [-1 1] * xl(2));
        else
            axis off;
        end
        
        if ii == 1, title(dv, 'interpreter', 'none'); end
        
    end
    
    %     subplot(length(D)+1,length(dvs),length(D)*length(dvs)+z)
    %     dvall = []; for ii = 1:5; dvall = [dvall D(ii).(dv)']; end
    %     lm = fitlm([x1all; x2all]', dvall', 'linear');
    %     lm.plotEffects; xl = get(gca, 'XLim'); set(gca, 'XLim', [-1 1] * xl(2));
    
end


fH(2) = figure; pos = get(fH(2), 'Position');

set(fH(2), 'Position', [pos(1) pos(2) 400 800], 'Color', 'w', ...
    'name', 'Winawer and Parvizi Figure 4E', 'NumberTitle', 'off');

for z = 1:length(dvs)
    
    xmx = 10.^ceil(log10(max(x1all)));
    xmn = 10.^floor(log10(min(x1all)));
    
    %xl = [8 400];  yl = [4  120];
    xl = [xmn xmx];  yl = [4  120];
    xl(1) = .1;
    xt = 10.^(log10(xmn):log10(xmx)); yt = [10 100];
    %xt = [10 100]; yt = [10 100];
    dv = dvs{z};
    
    for ii = 1:length(D)
        
        idx = isfinite(D(ii).which_drawing.val);
        if all(isnan(idx)), skip = true; else, skip = false; end
        
        subplot(length(D),length(dvs),(ii-1)*length(dvs)+z)
        set(gca, 'FontSize', 20); hold on;
        if ii == 1, title(dv, 'interpreter', 'none'); end
        
        if ~skip
            sz = D(ii).(dv).val(idx);
            sz = sz * 300 / max(sz);
            sz(sz == 0) = eps;
            
            scatter(D(ii).(x1).val(idx), D(ii).(x2).val(idx), sz, ...
                'MarkerFaceColor', opts.colors(ii,:),...
                'MarkerEdgeColor', 'k', 'LineWidth', 1)
            scatter(D(ii).(x1).val(idx), D(ii).(x2).val(idx), sz,...
                'MarkerEdgeColor', 'k', 'LineWidth', 2)
            if ii == length(D)
                xlabel(sprintf('Charge per pulse\n%s', D(ii).chargePerPulse.units));
            end
            if z == 1 && ii == 3, ylabel('Frequency (Hz)'); end
            
            axis([xl yl])
            set(gca, 'YScale', opts.yscale, 'XScale', opts.xscale, 'XTick', xt,'YTick', yt)
            axis square
            plot(xl, [15 15], 'k--', 0.7*[1 1], yl, 'k--')
        end
    end
    
    
end

end




