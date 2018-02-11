clear all
rand('seed',0);

cd('C:\Users\Ione Fine\Documents\Work\Science\Projects\Ione Fine\p2p-corticalstimulation');
addpath(genpath(C:\Users\Ione Fine\Documents\Work\Science\Projects\Ione Fine\p2p-corticalstimulation'));
%% cortical and visual space parameters

cortexSize = [30,30];  % %[height, width] Size of cortical maps (mm)
cortexCenter = [0,0]; % center of electrode array (mm on cortex)

c2r=2; % each degree is represented about ~280 microns
retinaSize = [cortexSize/c2r,cortexSize/c2r]; %  [height, width]
retinaCenter = [.75,0];

% define the resolution of the cortical and visual maps
pixpermm=20 ;% choose the resolution to sample in mm.
pixperdeg = pixpermm*c2r;

%% DEFINE CORTEX MESHGRID AND ORIENTATION/OD COLUMNS

xax = linspace(.5/pixpermm,cortexSize(2)-.5/pixpermm,cortexSize(2)*pixpermm)+cortexCenter(1)-cortexSize(2)/2;
yax = linspace(.5/pixpermm,cortexSize(1)-.5/pixpermm,cortexSize(1)*pixpermm)+cortexCenter(2)-cortexSize(1)/2;
[x,y] = meshgrid(xax,yax);
% Make the orientation and OD maps
[orientationMap,odMap] = makePinwheelODMaps(x,y, 0.863);
%% Show the orientation (pinwheel) map, cortical surface
figure(1); clf; subplot(1, 3,1)
image(xax,yax,64*(orientationMap+pi)/(pi*2));
set(gca,'YDir','normal');
colormap(hsv(64));
axis equal
axis tight
hold on
title('Orientation map');

%% Show the occular dominance map, cortical surface
subplot(1, 3,2)
imagesc(xax,yax,64*odMap);
set(gca,'YDir','normal');
colormap(gray(64));
axis equal
axis tight
title('Ocular dominance map');

%% VISUAL SPACE (RETINOTOPIC SPACE)

xrax = linspace(.5/pixperdeg,retinaSize(2)-.5/pixperdeg,retinaSize(2)*pixperdeg)+retinaCenter(1)-retinaSize(2)/2;
yrax = linspace(.5/pixperdeg,retinaSize(1)-.5/pixperdeg,retinaSize(1)*pixperdeg)+retinaCenter(2)-retinaSize(1)/2;
[xr,yr] = meshgrid(xrax,yrax);

%% RF to V1 mapping

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


%% lay out our electrode grid
load_windata;

D=getData(1);
eData=[ 7.8774	8.8141  1150/1000; ... % eccentricity, angle and size, extracted from Winaurer
    18.2743	23.5752 510/1000; ...
    2.4848	48.3118 1150/1000; ...
    7.8671	166.642 1150/1000; ...
    3.337	49.077 1150/1000];

s=3;
Z=eData(s,1).*exp(sqrt(-1)*eData(s,2)*pi/180); % convert to cortical co-ordinates (mm)
xe=real(Z); ye=imag(Z);
ERF = makeElectrodeRFs(xe,ye,eData(s,3),x,y);
subplot(1,3,3)
img = reshape(sum(ERF(:),2),size(x));
image(xax,yax,256*img/max(img(:)));
set(gca,'YDir','normal');
colormap(gray(256));
axis equal
axis tight
title('Electrode RFs');

% Draw the cortical grid on the figures
margin = 2; %mm
for i=1:3
    subplot(1, 3,i);   hold on
    plot(wAng,'-','Color',gridColor);
    plot(wRad','-','Color',gridColor);
    xlabel('mm');
    set(gca,'XLim',[min(xax(:))-margin,max(xax(:))+margin]);
    set(gca,'YLim',[min(yax(:))-margin,max(yax(:))+margin]);
end

ind=find(windata(:, 1) == s);
for w=1; %length(ind)
    duration(w) = windata(ind(w),6);  % seconds
    pulsewidth(w)= windata(ind(w),5)/1000; %.45/1000;  %sec
    current(w)= windata(ind(w),3);
    frequency(w) = windata(ind(w),4);
    totalCharge(w) = pulsewidth(w) * current(w);
    
    for i=1:4
        percept(i).img = zeros(size(xr)); % left eye
    end
    
    for pixNum = 1:length(x(:))
        if ERF(pixNum)>.02
            %map from cortex to retinal coordinates`
            vr = mapinv(p,x(pixNum)+sqrt(-1)*y(pixNum));
            
            RFsize = min(.1+.1*abs(vr),1);  %based figure in Simoncelli and Freeman
            RForientation = orientationMap(pixNum);
            
            % Oriented gaussian
            xrot = (xr-real(vr))*cos(RForientation)+(yr-imag(vr))*sin(RForientation);
            yrot = -(xr-real(vr))*sin(RForientation) + (yr-imag(vr))*cos(RForientation);
            RF.ar = .1;  % aspect ratio of V1 RFs
            G = exp(-( xrot.^2/(2*RFsize*RF.ar) + yrot.^2/(2*RFsize)));
            percept(1).img =  percept(1).img  + odMap(pixNum)*ERF(pixNum)*G;    % left eye
            percept(2).img =  percept(2).img  + (1-odMap(pixNum))*ERF(pixNum)*G; % right eye
            
            % scoreboard version
            G_SB = exp(-( xrot.^2/(0.0001) + yrot.^2/(.0001)));
            percept(3).img =  percept(3).img + ERF(pixNum)*G_SB;    % left eye
            percept(4).img =  percept(4).img + ERF(pixNum)*G_SB; % right eye
        end
    end
    
    titleStr = {'Left Eye','Right Eye'};
    figure(2); clf;
    for i=1:4
        subplot(2,2,i)
        img=256*3*(percept(i).img./max(percept(i).img(:)));
        img(img>256)=256;
        image(xrax,yrax,img);
        title(titleStr(mod(i, 2)+1));
        colormap(gray(256)); hold on
        plot(zAng,'-','Color',gridColor);
        plot(zRad','-','Color',gridColor);
        xlabel('deg');
        axis equal; axis tight
        set(gca,'XLim',[-.5,1.25]);
        set(gca,'YLim',[-1.25,1.25]);
    end
end

%% save files
cd('C:\Users\Ione Fine\Documents\Work\Science\Projects\Ione Fine\p2p-corticalstimulation\images');
for f=1:4
    gcf=figure(f);
    filename=[ECONFIG, num2str(eSize*1000), '_', num2str(f), '.jpg'];
    saveas(gcf,filename)
end


