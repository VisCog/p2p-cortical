% BeauchampFig4.m
% 
% Finds best cortical mapping (a,k,squish) and electrode array location
% (ang,xc,yc) for an array on our V1 map to match the projections of
% phosphenes from Beauchamp (2020).

% Note: electrode spacing dx is fixed at a 2mm for both x and y.

% load in phosphene parameters
T = readtable('datasets/Beauchamp_2020_data.xlsx','Sheet','ElectrodeLocations');
eid =strcmp(T.experiment,'Figure 4');
T = T(eid,:)
% Phosphene locations in visual space
id = strcmp(T(:,:).patient,'YBN');
vx = T(id,:).x';
vy = T(id,:).y';

% model parameters
% kmapping parameters
p.k = 15;
p.a = 0.5; %.1474;
p.squish = 1; %.63;
% electrode location and spacing on cortex
p.ang = -30*pi/180;
p.xc = -44;
p.yc = -8;
p.dx = 2;  %mm

% Axis range for plotting:
% zoomed in:
p.eccList = [4,8,16];
p.angList = 0:pi/8:pi/2;

% zoomed out:
% p.eccList = [0.01,0.25,.5,1,2,4,8,16];
% p.angList = -pi/2:pi/8:pi/2;

% Draw the phosphene locations in visual space
figure(1)
clf
hold on

% Draw the grid
ecc = repmat(p.eccList,101,1);
ang = repmat(linspace(min(p.angList),max(p.angList),101)',1,length(p.eccList));
gridx = ecc.*cos(ang);
gridy = ecc.*sin(ang);
plot(gridx,gridy,'k-','Color',[.5,.5,.5]);

ang = repmat(p.angList,101,1);
ecc = repmat(exp(linspace(-20,log(max(p.eccList)),101))',1,length(p.angList));
gridx = ecc.*cos(ang);
gridy = ecc.*sin(ang);
plot(gridx,gridy,'k-','Color',[.5,.5,.5]);

% Plot the phosphenes in visual space (blue)
drawLocations(vx,vy,[0,0,1]);
axis equal
title('Visual Space');
xlabel('x (deg)');
ylabel('y (deg)');

% Fit the array
p = fit('p2p_c.fitElectrodeGrid',p,{'ang','a','k','squish','xc','yc'},vx,vy);
p.shift = p.k*log(p.a);  % shouldn't really be a parameter

% Get the predictions
[err,predcx,predcy,cx,cy] = p2p_c.fitElectrodeGrid(p,vx,vy);

% Draw the projected phosphenes (blue) and the best fitting array (green)
figure(2)
clf
hold on
% Plot the grid
ecc = repmat(p.eccList,101,1);
ang = repmat(linspace(min(p.angList),max(p.angList),101)',1,length(p.eccList));
x = ecc.*cos(ang);
y = ecc.*sin(ang);
[gridx,gridy] = p2p_c.v2c_real(p,x,y);
plot(gridx,gridy,'k-','Color',[.5,.5,.5]);

ang = repmat(p.angList,101,1);
ecc = repmat(exp(linspace(log(min(p.eccList)),log(max(p.eccList)),101))',1,length(p.angList));
x = ecc.*cos(ang);
y = ecc.*sin(ang);
[gridx,gridy] = p2p_c.v2c_real(p,x,y);
plot(gridx,gridy,'k-','Color',[.5,.5,.5]);

title('Cortical Space');
xlabel('x (mm)');
ylabel('y (mm)');
axis equal

% Projected phosphenes
drawLocations(cx,cy,[0,0,1]);

% Best fitting array
drawLocations(predcx,predcy,[0,1,0]);

% project best fitting electrode array back to visual space
[predvx,predvy] = p2p_c.c2v_real(p,predcx,predcy);

% Plot the projected best fitting array in visual space (green)
figure(1)
drawLocations(predvx,predvy,[0,1,0]);

function drawLocations(x,y,col)

for i=1:length(x)
    plot(x(i),y(i),'ko','MarkerSize',14,'Color',col,'MarkerFaceColor','none');
    text(x(i),y(i),num2str(i),'HorizontalAlignment','center','VerticalAlignment','middle','Color',col);
end

end
