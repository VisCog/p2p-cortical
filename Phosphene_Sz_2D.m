% Phosphene_Sz_2D.m
%
% Packs phosphenes in visual space that vary in size as a parametric
% function of eccentricity.  Phosphene centers are then projected into
% Schwartz cortical space to show how spacing of electrodes are less
% densely packed near the fovea.

% Data for phosphene size as a function of eccentricity:

% cortical map parameters
c.a = .1;
c.k = 15;
c.shift = c.k*log(c.a);

option = 4;

switch option
    case 1
        % option 1: use the model's simulated phosphene sizes
        T = readtable('datasets/PhospheneSz_vs_Ecc.xlsx');
        Tid = T.e_size==T.e_size(1);
        
        % hack - interpolate back to zero degrees eccentricity
        e.ecc = [0;T(Tid,:).eccentricity]; % eccentricity
        e.sz = [2*T(Tid(1),:).size_sigma;2*T(Tid,:).size_sigma];   % '2*' because sigma is the radius, 'sz' is the diameter of the phosphene
        n=400; % number of electrodes/phosphenes
        
        % You should fiddle with these physics parameters to bring phosphenes
        % together but not overlapping.
        p.force = .125;
        p.friction = .25;
        p.g = .001;
        nIter = 10000;
        % e.sz =  .007*e.ecc.^2+1; % diameter of phosphene
        
        
    case 2
        % option 2: use the winawer straight line in figure D:
        
        e.ecc = [0,30];
        e.sz =  0.2620*e.ecc +  0.1787;  % straight line from Bosking figure D
       % e.sz = 0.2922*e.ecc+.2862;  % polyfit from GrabData
        n = 265;
        % You should fiddle with these physics parameters to bring phosphenes
        % together but not overlapping.
        p.force = .11;
        p.friction = 0;
        p.g = .0017;
        nIter = 50000;
    case 3
        
        sep0 = 4;  % sep0=4 gives a similar slope as option 2
        e.ecc = [0,40];
        e.sz =  sep0*(e.ecc+c.a)/c.k;  % straight line that predicts constant separation
        %e.sz = 0.2667*e.ecc+0.1333;  % slopes and intercepts for sep0 = 4
        n = 275;
        % You should fiddle with these physics parameters to bring phosphenes
        % together but not overlapping.
        p.force = .11;
        p.friction = 0;
        p.g = .0014;
        nIter = 20000;
    case 4
        % option 2: use the winawer straight line in figure D:
        
        e.ecc = [0,30];
        %e.sz =  .5*(0.2620*e.ecc +  0.1787);  % straight line from Bosking figure D
        % e.sz = 0.2922*e.ecc+.2862;  % polyfit from GrabData
        e.sz = 0.1305*e.ecc*0.524; % bosking simulation
        n = 1000;
        % You should fiddle with these physics parameters to bring phosphenes
        % together but not overlapping.
        p.force = .17;
        p.friction = 0;
        p.g = .0014;
        nIter = 1000;
end


p.electrodeSize =1.25;

% For polar axes on plots
p.eccList = [0,.5,1,2,4,8,16,32];
p.angList = [90.1,135,180,225,270]*pi/180;

%%
% Plot size vs eccentricity
ecc = linspace(0,max(e.ecc),101);
figure(1)
sz = ecc2sz(e,ecc);
plot(ecc,sz,'ko-','MarkerSize',3,'MarkerFaceColor',[.5,.5,.5]);
set(gca,'YLim',[0,max(sz)*1.1]);
grid
xlabel('Eccentricity (deg)');
ylabel('Phosphene size (deg)');

%%
% Pack phosphenes

% initial phosphene positions (spiraling out)
ecc = (linspace(0,9,n))'.^2;
ang = linspace(0,n*15*pi,n)';
p.x = ecc.*cos(ang);
p.y = ecc.*sin(ang);

% initial velocities for migrating phosphenes
dx = zeros(n,1);
dy = zeros(n,1);

% iterate though 'time' to draw phosphenes toward fovea but keep them from
% overlapping

for t=1:nIter
    % matrix of distances
    distx = repmat(p.x,1,n)-repmat(p.x',n,1);
    disty = repmat(p.y,1,n)-repmat(p.y',n,1);
    D = sqrt(distx.^2+disty.^2);
    % forces
    sz = ecc2sz(e,sqrt(p.x.^2+p.y.^2));
    sz = min(sz,10);  % limit max size while packing
    
    r2 = repmat(sz/2,1,n)+repmat(sz'/2,n,1);  % radius = sz/2
    id = D < r2- eye(n);
    
    % push phosphenes away from each other but toward the center
    dx = p.friction*dx+-sum(p.force*distx.*id)'-p.g*p.x;
    dy = p.friction*dy+-sum(p.force*disty.*id)'-p.g*p.y;
    
    p.x = p.x+dx;
    p.y = p.y+dy;
end

%%
% Draw phosphenes in visual space

% define colors for each electrode/phosphene

% Plot left visual field electrodes:
id = p.x<0;

cmap = hsv(100);
col = 0.9*ones(n,3);  % gray
col(id,:) = cmap(ceil(sz(id)*10),:);

% draw the phosphenes in visual coordinates
figure(2)
clf
ecc = repmat(p.eccList,101,1);
ang = repmat(linspace(pi/2,3*pi/2,101)',1,length(p.eccList));
gridx = ecc.*cos(ang);
gridy = ecc.*sin(ang);
plot(gridx,gridy,'k-','Color',[.5,.5,.5]);

hold on
ang = repmat(p.angList,101,1);
ecc = repmat(exp(linspace(-20,log(max(p.eccList)),101))',1,length(p.angList));
gridx = ecc.*cos(ang);
gridy = ecc.*sin(ang);
plot(gridx,gridy,'k-','Color',[.5,.5,.5]);

% unit circle
cx = cos(linspace(-pi,pi,61));
cy = sin(linspace(-pi,pi,61));



for i=1:length(p.x)
    patch(p.x(i)+sz(i)/2*cx,p.y(i)+sz(i)/2*cy,col(i,:),'EdgeColor','w');
end
axis equal

m = max(max(abs(p.x)+sz),max(abs(p.y)+sz));
m=38;
set(gca,'xLim',[-m,m]);
set(gca,'YLim',[-m,m]);

%%
% Draw the corresponding 'electrodes' in cortical space

figure(3)
clf

ecc = repmat(p.eccList,101,1);
ang = repmat(linspace(90.1,270,101)'*pi/180,1,length(p.eccList));
x = ecc.*cos(ang);
y = ecc.*sin(ang);
[gridx,gridy] = p2p_c.v2c_real(c,x,y);
plot(gridx,gridy,'k-','Color',[.5,.5,.5]);

hold on
ang = repmat(p.angList,101,1);
ecc = repmat(exp(linspace(-20,log(max(p.eccList)),101))',1,length(p.angList));
x = ecc.*cos(ang);
y = ecc.*sin(ang);

[gridx,gridy] = p2p_c.v2c_real(c,x,y);
plot(gridx,gridy,'k-','Color',[.5,.5,.5]);

id = p.x<0;
p.z = p.x+sqrt(-1)*p.y;
[vx,vy] = p2p_c.v2c_real(c,p.x,p.y);


% for circles:
% cx = cos(linspace(-pi,pi,61));
% cy = sin(linspace(-pi,pi,61));

% for squares:
cx = cos(-pi/4:pi/2:5*pi/4);
cy = sin(-pi/4:pi/2:5*pi/4);
sz = ecc2sz(e,sqrt(p.x.^2+p.y.^2));

hold on
for i=1:length(p.x)
    if id(i)
       h= patch(vx(i)+p.electrodeSize/2*cx,vy(i)+p.electrodeSize/2*cy,col(i,:));
    end
end
axis equal

% find axis limits

[limx,limy] = p2p_c.v2c_real(c,[-.01,-.01],max(p.eccList)*1.2*[1,-1]);
set(gca,'XLim',[-3,limx(1)]);
set(gca,'YLim',limy*1.2);

%%
% Colorbar

figure(4)
clf

maxSz = 8; %deg

img = repmat((1:(ceil(maxSz*10)))',1,2);

image(1,linspace(0,maxSz,size(col,1)),img);
colormap(cmap)
set(gca,'YDir','normal');
set(gca,'XTick',[]);
ylabel('Phosphene size (deg)');
set(gca,'FontSize',18);
axis equal
axis tight
set(gcf,'PaperPosition',[1,1,1,1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Support function  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sz  = ecc2sz(e,ecc)
% interpolate data to get phosphene size
ecc = min(ecc,max(e.ecc));
sz = interp1(e.ecc,e.sz,ecc);
end


