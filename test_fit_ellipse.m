% test_fit_ellipse.m

[v.X,v.Y] = meshgrid(linspace(-10,10,1000));

start_p.x0 = 1;
start_p.y0 =-1;
start_p.theta = 2*pi*rand(1);
start_p.sigma_x = randn(1)*4;
start_p.sigma_y = randn(1)*4;

img = zeros(size(v.X));
[err,G] = fit_ellipse_to_phosphene(start_p,img,v);
img = double(G>exp(-1));

figure(1)
clf
imagesc(v.X(1,:),v.Y(:,1),img);
axis equal
axis tight
grid
colormap(hot)
% 
% clear p
% p.x0 = sum(img(:).*v.X(:))/sum(img(:)); % COM for x
% p.y0 = sum(img(:).*v.Y(:))/sum(img(:)); % COM for y
% p.theta = 0;
% p.sigma_x = 1;
% p.sigma_y = 1;
% p.shutup = 1;
% 
% p = fitcon('fit_ellipse_to_phosphene',p,{'theta','sigma_x>0','sigma_y>0'},img,v);
% p = fitcon('fit_ellipse_to_phosphene',p,{'theta','sigma_x','sigma_y','x0','y0'},img,v);
% p.theta = mod(p.theta,pi);
% 
% if p.sigma_x < p.sigma_y
%     tmp = p.sigma_x;
%     p.sigma_x = p.sigma_y;
%     p.sigma_y = tmp;
%     p.theta = pi-p.theta;
% end


% [err,G] = fit_ellipse_to_phosphene(p,img,v);

% figure(2)
% clf
% image(v.X(1,:),v.Y(:,1),(G+img)*128);
% colormap([gray(128);hot(128)]);
% axis equal
% axis tight
% grid

% figure(3);
% clf
% imagesc(v.X(1,:),v.Y(:,1),(G>exp(-1))+img);
% colormap([0,0,1;1,0,0;1,1,0]);
% axis equal
% axis tight
% grid

%% using moments

M00 = sum(sum(img));

M10 = sum(sum(v.X.*img));
M01 = sum(sum(v.Y.*img));

M11 = sum(sum(v.X.*v.Y.*img));

M20 = sum(sum(v.X.^2.*img));
M02 = sum(sum(v.Y.^2.*img));

x0 = M10/M00;
y0 = M01/M00;

mu20 = M20/M00 - x0^2;
mu02 = M02/M00 - y0^2;
mu11 = M11/M00 - x0*y0;


a = (mu20+mu02)/2;
b = .5*sqrt(4*mu11^2+(mu20-mu02)^2);

lambda_1 = a+b;
lambda_2 = a-b;

p.x0 = x0;
p.y0 = y0;
p.theta = -.5*atan2(2*mu11,mu20-mu02)
p.sigma_x = sqrt(2*lambda_1);
p.sigma_y = sqrt(2*lambda_2);

[err,G] = fit_ellipse_to_phosphene(p,img,v);




figure(2);
clf
imagesc(v.X(1,:),v.Y(:,1),G>exp(-1));
axis equal
axis tight
grid

figure(3);
clf
imagesc(v.X(1,:),v.Y(:,1),(G>exp(-1))+img);
colormap([0,0,1;1,0,0;1,1,0]);
axis equal
axis tight
grid
colorbar
