% clear all
close all

rand('seed',0);

%% parameters
cortexSize = [4,4];  % %[height, width] Size of cortical maps (mm)
cortexCenter = [0,0]; % center of electrode array (mm on cortex)

c2r=2; % how many degrees of visual angle do you want to represent per mm of cortex
retinaSize = [cortexSize/c2r,cortexSize/c2r]; %  [height, width]
retinaCenter = [.75,0]; % don't change

eSize  = .5;   % radius of electrode 'RF' radius in mm, 1000 microns/mm
pixpermm=round(1/(eSize/3)); % choose the resolution on the cortex to sample in mm.
pixperdeg = pixpermm*c2r; % use a similar sampling for visual space

%eList = [1:64];
eList = [10,15,50, 55]; % which electrodes are going to be stimulated
% eList = [1,8,57,64];

%% CORTEX

% Define the cortex meshgrid
xax = linspace(.5/pixpermm,cortexSize(2)-.5/pixpermm,cortexSize(2)*pixpermm)+cortexCenter(1)-cortexSize(2)/2;
yax = linspace(.5/pixpermm,cortexSize(1)-.5/pixpermm,cortexSize(1)*pixpermm)+cortexCenter(2)-cortexSize(1)/2;
[x,y] = meshgrid(xax,yax);

% Define the electrode centers (mm cortex)
nElectrodes = 8^2;
[xe,ye] = meshgrid(linspace(eSize*2,cortexSize(2)-eSize*2,sqrt(nElectrodes)),linspace(eSize*2,cortexSize(1)-eSize*2,sqrt(nElectrodes)));
xe = xe+cortexCenter(1)-cortexSize(2)/2;
ye = ye+cortexCenter(1)-cortexSize(1)/2;

% Make the electrode array electric fields
EF = makeElectrodeEFs(xe,ye,eSize,x,y);

% Make the orientation and OD maps on cortex
[orientationMap,odMap] = makePinwheelODMaps(x,y);
RF_ar = .1;  % aspect ratio of orientation tuned cortical receptive fields


%% VISUAL SPACE (RETINOTOPIC SPACE)

% Define the meshgrid for visual space
xrax = linspace(.5/pixperdeg,retinaSize(2)-.5/pixperdeg,retinaSize(2)*pixperdeg)+retinaCenter(1)-retinaSize(2)/2;
yrax = linspace(.5/pixperdeg,retinaSize(1)-.5/pixperdeg,retinaSize(1)*pixperdeg)+retinaCenter(2)-retinaSize(1)/2;
[xr,yr] = meshgrid(xrax,yrax);

%% Show the orientation (pinwheel) map on the cortical surface
figure(1); clf
image(xax,yax,64*(orientationMap+pi)/(pi*2));
set(gca,'YDir','normal');
colormap(hsv(64));
axis equal
axis tight
hold on
title('Orientation map');

%% Show the ocular dominance map on the cortical surface
figure(2); clf
imagesc(xax,yax,64*odMap);
set(gca,'YDir','normal');
colormap(gray(64));
axis equal
axis tight
title('Ocular dominance map');

%% Show the sums of the electrode electric fields on the cortical surface
figure(3); clf
img = reshape(sum(EF(:, eList),2),size(x));
image(xax,yax,256*img/max(img(:)));
set(gca,'YDir','normal');
colormap(gray(256));
axis equal
axis tight
title('Electrode EFs');

%% VISUAL SPACE to V1 mapping

% typical log z transformation parameters (Based on Duncan and Boynton)
p.k = -20; %scale
p.a = .5; %fovea expansion

gridColor = .5*[1,1,1];

% Make the grid in visual space coordinates
angList = [-90:20:90];
radList = [.5,1,2,4,8];
n = 201;

zAng = linspace(0,max(radList),n)'*exp(sqrt(-1)*angList*pi/180);
zRad = radList'*exp(sqrt(-1)*linspace(-90,90,n)*pi/180);

% project the grid to cortical coordinates (mapinv goes the other way)
wAng = map(p,zAng);
wRad =map(p,zRad);

% Draw the grid on the cortical surface figures
for i=1:3
    figure(i)
    hold on
    plot(wAng,'-','Color',gridColor);
    plot(wRad','-','Color',gridColor);
    xlabel('mm');
end

%% BUILD THE EXPECTED PERCEPTS

perceptL = zeros(size(xr)); % left eye
perceptR = zeros(size(xr)); % right eye
perceptL_scoreboard = zeros(size(xr)); % left eye
perceptR_scoreboard = zeros(size(xr)); % right eye


vr = mapinv(p,x(:)+sqrt(-1)*y(:));
sizeMap = min(.1+.1*abs(vr),1);  %based figure in Simoncelli and Freeman
tic


for pixNum = 1:length(x(:))  
    % Oriented gaussian
    xrot = (xr-real(vr(pixNum)))*cos(orientationMap(pixNum))+(yr-imag(vr(pixNum)))*sin(orientationMap(pixNum));
    yrot = -(xr-real(vr(pixNum)))*sin(orientationMap(pixNum)) + (yr-imag(vr(pixNum)))*cos(orientationMap(pixNum));
    G_scoreboard = exp(-( xrot.^2/(0.001) + yrot.^2/(.001))); % make the receptive fields really tiny
    G = exp(-( xrot.^2/(2*sizeMap(pixNum)*RF_ar) + yrot.^2/(2*sizeMap(pixNum))));
    for eNum = eList
        if EF(pixNum,eNum)>.1
            perceptL =   perceptL + reshape(odMap(pixNum)*EF(pixNum,eNum)*G, size(perceptL)); % right eye
            perceptR =   perceptR +  reshape((1-odMap(pixNum))*EF(pixNum,eNum)*G, size(perceptR)); % right eye
            % scoreboard version
            perceptL_scoreboard = perceptL_scoreboard + reshape(EF(pixNum,eNum)*G_scoreboard, size(perceptL_scoreboard));
            
        end
    end
end
disp(toc)
perceptR_scoreboard =perceptL_scoreboard;
% both eyes the same


%%
titleStr = {'Left Eye','Right Eye'};
figure(4)
clf;
for i=1:4
    subplot(2,2,i)
    hold on
    if i<3
        if mod(i,2)==0
            img=256*3*(perceptL_scoreboard./max(perceptL_scoreboard(:)));
        else
            img=256*3*(perceptR_scoreboard./max(perceptR_scoreboard(:)));
        end
        img(img>256)=256;
        image(xrax,yrax,img);
        title(titleStr(i));
    else
        if mod(i,2)==0
            img=256*2*(perceptL./max(perceptL(:)));
        else
            img=256*2*(perceptR./max(perceptR(:)));
        end
        img(img>256)=256;
        image(xrax,yrax,img);
    end
    colormap(gray(256))
    plot(zAng,'-','Color',gridColor);
    plot(zRad','-','Color',gridColor);
    xlabel('deg');
    axis equal
    axis tight
    set(gca,'XLim',[min(xrax(:)),max(xrax(:))]);
    set(gca,'YLim',[min(yrax(:)),max(yrax(:))]);
end

margin = 2; %mm
for i=1:3
    figure(i);
    set(gca,'XLim',[min(xax(:))-margin,max(xax(:))+margin]);
    set(gca,'YLim',[min(yax(:))-margin,max(yax(:))+margin]);
end

%% save files
cd('C:\Users\Ione Fine\Documents\Work\Science\Projects\Ione Fine\p2p-corticalstimulation\images');
for f=1:4
    gcf=figure(f);
    filename=['corticalpros_', num2str(eSize*1000), '_', num2str(f), '.jpg'];
    saveas(gcf,filename)
end
