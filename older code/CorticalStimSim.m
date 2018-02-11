clear all
rand('seed',0);

cd('C:\Users\Ione Fine\Documents\Work\Science\Projects\Ione Fine\p2p-corticalstimulation');

%% parameters

cortexSize = [2,2];  % %[height, width] Size of cortical maps (mm)
cortexCenter = [0,0]; % center of electrode array (mm on cortex)

c2r=2; % each degree is represented about ~280 microns, so a scaling factor of 3 to 4, reduce from there so enough visual field is shown
retinaSize = [cortexSize/c2r,cortexSize/c2r]; %  [height, width]
retinaCenter = [.75,0];

ECONFIG='Winawer';
if strcmp(ECONFIG, '2sight')
    eSize = 500 / 1000;   % size of electrode 'RF' radius in mm, second sight
    eSep = 520 / 1000;
    nElectrodes = 8^2;
    eLoc=[0 0]; % center of the electrode array in cortical co-ordinates (mm)
elseif strcmp(ECONFIG, 'Troyk')
    eSize  = 15/1000;   % size of electrode 'RF' radius in mm, second sight
    eSep = 400 / 1000;
    nElectrodes = 8^2;
    eLoc=[0 0]; % center of the electrode array in cortical co-ordinates (mm)
    elseif strcmp(ECONFIG, 'Winawer')
    eSize  = 2300/1000;   % size of electrode 'RF' radius in mm, second sight
    eSep = 4000 / 1000;
    nElectrodes = 8^2;
    eLoc=[0 0]; % center of the electrode array in cortical co-ordinates (mm)
    
elseif strcmp(ECONFIG, 'Default')
    eSize=500/1000; % unrealistically large but fast
    eSep = 2000 / 1000;
    nElectrodes = 2^2;
    eLoc=[0 0]; % center of the electrode array in cortical co-ordinates (mm)
else
    error('ECONFIG incorrectly defined, line 11');
end

% define the resolution of the cortical and visual maps
pixpermm=80 ; max(round(1/(eSize/2)), 30); % choose the resolution to sample in mm.
pixperdeg = pixpermm*c2r;

%% really basic current integration model
% Model Parameters
% core parameters

%% DEFINE CORTEX MESHGRID

xax = linspace(.5/pixpermm,cortexSize(2)-.5/pixpermm,cortexSize(2)*pixpermm)+cortexCenter(1)-cortexSize(2)/2;
yax = linspace(.5/pixpermm,cortexSize(1)-.5/pixpermm,cortexSize(1)*pixpermm)+cortexCenter(2)-cortexSize(1)/2;
[x,y] = meshgrid(xax,yax); 

%% Define the electrode centers (mm cortex)

s=(eSep * (sqrt(nElectrodes)-1)) / 2;
[xe, ye] = meshgrid(linspace(-s, s, sqrt(nElectrodes)), linspace(-s, s, sqrt(nElectrodes)));
xe=xe+eLoc(1); ye=ye+eLoc(2);

% Make the electrode array RF's
ERF = makeElectrodeRFs(xe,ye,eSize,x,y);

RF.ar = .1;  % aspect ratio of pixel RFs
eList = [1]; % corner electrodes
%eList = [10,15,50, 55]; % just inside the corners

% Make the orientation and OD maps
[orientationMap,odMap] = makePinwheelODMaps(x,y, 0.863);


%% VISUAL SPACE (RETINOTOPIC SPACE)

%% Define the retina meshgrid
xrax = linspace(.5/pixperdeg,retinaSize(2)-.5/pixperdeg,retinaSize(2)*pixperdeg)+retinaCenter(1)-retinaSize(2)/2;
yrax = linspace(.5/pixperdeg,retinaSize(1)-.5/pixperdeg,retinaSize(1)*pixperdeg)+retinaCenter(2)-retinaSize(1)/2;
[xr,yr] = meshgrid(xrax,yrax);

%% Show the orientation (pinwheel) map, cortical surface
figure(1)
clf
image(xax,yax,64*(orientationMap+pi)/(pi*2));
set(gca,'YDir','normal');
colormap(hsv(64));
axis equal
axis tight
hold on
title('Orientation map');

%% Show the occular dominance map, cortical surface
figure(2)
clf
imagesc(xax,yax,64*odMap);
set(gca,'YDir','normal');
colormap(gray(64));
axis equal
axis tight
title('Ocular dominance map');

%% Show the sums of the electrode RFs, cortical surface
figure(3)
clf
img = reshape(sum(ERF(:, eList),2),size(x));
image(xax,yax,256*img/max(img(:)));
set(gca,'YDir','normal');
colormap(gray(256));
axis equal
axis tight
title('Electrode RFs');

return
%%
% RF to V1 mapping

% typical log z transformation parameters (Based on Duncan and Boynton)
p.k = 20; %scale
p.a = 0.5; %fovea expansion

gridColor = .5*[1,1,1];

%Make the grid in retinal coordinates
angList = [-90:20:90];
radList = [1.25,2.5,5,10];
n = 201;

zAng = linspace(0,max(radList),n)'*exp(sqrt(-1)*angList*pi/180);
zRad = radList'*exp(sqrt(-1)*linspace(-90,90,n)*pi/180);

% project the grid to cortical coordinates (mapinv goes the other way)
wAng = map(p,zAng);
wRad =map(p,zRad);

% Draw the grid on the figures
for i=1:3
    figure(i)
    hold on
    plot(wAng,'-','Color',gridColor);
    plot(wRad','-','Color',gridColor);
    xlabel('mm');
end

figure(10)
hold on
plot(wAng,'-','Color',gridColor);
plot(wRad','-','Color',gridColor);
xlabel('mm');

%%

perceptL = zeros(size(xr)); % left eye
perceptR = zeros(size(xr)); % right eye
perceptL_SB = zeros(size(xr)); % left eye
perceptR_SB = zeros(size(xr)); % right eye

% load_windata;
% for s=1:5
%     ind=find(windata(:, 1) ==s);
%     
% for w=1:length(ind)
%     STIM.dur = windata(ind(w),6);  % seconds
%     STIM.pulsedur = windata(ind(w),5)/1000; %.45/1000;  %sec
% 
%     STIM.amp = windata(ind(w),3);
%     STIM.freq = windata(ind(w),4);


for eNum = eList % length(xe(:));
    for pixNum = 1:length(x(:))
        if ERF(pixNum,eNum)>.1
            
            figure(3)
            %    plot(x(pixNum),y(pixNum),'w.');
            
            %map from cortex to retinal coordinates
            vr = mapinv(p,x(pixNum)+sqrt(-1)*y(pixNum));
            
            RFsize = min(.1+.1*abs(vr),1);  %based figure in Simoncelli and Freeman
            RForientation = orientationMap(pixNum);
            
            % Oriented gaussian
            xrot = (xr-real(vr))*cos(RForientation)+(yr-imag(vr))*sin(RForientation);
            yrot = -(xr-real(vr))*sin(RForientation) + (yr-imag(vr))*cos(RForientation);
            G = exp(-( xrot.^2/(2*RFsize*RF.ar) + yrot.^2/(2*RFsize)));
            perceptL =  perceptL + odMap(pixNum)*ERF(pixNum,eNum)*G;    % left eye
            perceptR =  perceptR + (1-odMap(pixNum))*ERF(pixNum,eNum)*G; % right eye
            
            % scoreboard version
            G_SB = exp(-( xrot.^2/(0.0001) + yrot.^2/(.0001)));
            perceptL_SB =  perceptL_SB+ ERF(pixNum,eNum)*G_SB;    % left eye
            perceptR_SB =  perceptR_SB+ ERF(pixNum,eNum)*G_SB; % right eye
        end
    end
end


%%
titleStr = {'Left Eye','Right Eye'};
figure(4)
clf;
for i=1:4
    subplot(2,2,i)
    hold on
    if i<3
        if mod(i, 2)==0
            img=256*3*(perceptL_SB./max(perceptL_SB(:)));
        else
            img=256*3*(perceptR_SB./max(perceptR_SB(:)));
        end
        img(img>256)=256;
        image(xrax,yrax,img);
        title(titleStr(i));
    else
        if mod(i, 2)==0
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
    %set(gca,'XLim',[-.5,max(xrax(:))]);
    %set(gca,'YLim',[min(yrax(:)),max(yrax(:))]);
    set(gca,'XLim',[-.5,1.25]);
    set(gca,'YLim',[-1.25,1.25]);
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
    filename=[ECONFIG, num2str(eSize*1000), '_', num2str(f), '.jpg'];
    saveas(gcf,filename)
end
