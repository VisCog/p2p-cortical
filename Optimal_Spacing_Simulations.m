% Optimal_Spacing_Simulations.m
%
% Packs phosphenes in visual space that vary in size as a parametric
% function of eccentricity.  Phosphene centers are then projected into
% Schwartz cortical space to show how spacing of electrodes are less
% densely packed near the fovea.

% Data for phosphene size as a function of eccentricity:

clear

maxEcc = 40;
maxEcc2 = 2.^ceil(log2(maxEcc));
colScale = 30; % scaling so using the full colormap
e.ecc = [0,maxEcc];
%e.sz = 0.4218*e.ecc*0.537; % keliris simulation, ignore this

slope = 0.0524;
intercept = 0.1305;
e.sz = (intercept+e.ecc*slope); % bosking simulation


% cortical map parameters
c.a = .5;
c.k = 15;
c.shift = c.k*log(c.a);
c.squish = 1;
p.electrodeSize =1;

% For polar axes on plots
p.eccList = [0,2.^(1:log2(maxEcc2))];
p.angList = [90.1,135,180,225,270]*pi/180;

%%
% Plot size vs eccentricity
ecc = linspace(0,maxEcc2,101);
figure(1);
sz = ecc2sz(e,ecc);
idx = find(ecc<=maxEcc);
plot(ecc(idx),sz(idx),'-','Color',[.5,.5,.5],"LineWidth",2);
set(gca,'YLim',[0,max(sz)*1.1]);
grid
xlabel('Eccentricity (deg)');
ylabel('Phosphene size (deg)');


%%
% Pack phosphenes

% generate optimal spacing along the horizontal meridian by adding up
% phosphene sizes
xi = 0;
si = intercept;
i = 1;
while(xi(end)+si(end)<maxEcc)
   si(i) = xi(i)*slope+intercept;
   xi(i+1) = (xi(i)+(si(i)+intercept)/2)/(1-slope/2);
   i=i+1;
end
si(i) = xi(i)*slope+intercept;
   
% make concentric rings of phosphene at each spacing distance   
ni = floor(2*pi*xi./si)';
p.x = [0];
p.y = [0];
for i=1:length(xi)
    a = linspace(0,2*pi,ni(i)+1)';
    a = a(1:end-1);
    p.x = [p.x;xi(i).*cos(a)];
    p.y = [p.y;xi(i).*sin(a)];    
end

rad = sqrt(p.x.^2+p.y.^2);
sz = slope*rad + intercept;
n= length(p.x);

% Hack to get the foveal phosphene inside the plotting range
p.x = p.x-.0001;

%%
% Draw phosphenes in visual space

% define colors for each electrode/phosphene

% Plot left visual field electrodes:
id = p.x<=0;

ncol = 100;
cmap = hsv(100);
col = 0.9*ones(n,3);  % gray
col(id,:) = cmap(min(ceil(sz(id)*colScale),100),:);

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
    patch(p.x(i)+sz(i)/2*cx,p.y(i)+sz(i)/2*cy,col(i,:));
end
axis equal

m = 42; %max(max(abs(p.x)+sz),max(abs(p.y)+sz));
%m=maxEcc2*1.1;
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
        patch(vx(i)+p.electrodeSize/2*cx,vy(i)+p.electrodeSize/2*cy,col(i,:));
    end
end
axis equal
% weirdly vx and vy are the locations on the cortical surface

% find axis limits

[limx,limy] = p2p_c.v2c_real(c,[-.01,-.01],max(p.eccList)*1.2*[1,-1]);
set(gca,'XLim',[-3,limx(1)]);
set(gca,'YLim',limy*1.2);


%%
% Colorbar

figure(4); clf

maxSz = ceil(max(sz)); %deg
img = repmat((1:(colScale*ceil(maxSz)))',1,7);
image(1,linspace(0,maxSz*colScale,size(col,1)),img);
colormap(cmap)
set(gca,'YDir','normal');
set(gca,'YTick',colScale*linspace(min(sz), max(sz), 4));
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


