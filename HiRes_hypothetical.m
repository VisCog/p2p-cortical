% HiRes_hypothetical.m

clear all; close all; clear mex

rng(1)  % fix the random number generator. This affects the ocular dominance/orientation maps
c.efthr = 0.15; % what magnitude of electric field goes through the model, just a speed thing not important
v.drawthr = .015; % what phosphene strength is visible, threshold of visibility/oval drawing, not important for you
nelect = 70;


%% load image
%
% in_filename = 'TeddySmiling.jpg';
% in_img = imread(in_filename); in_img = mean(double(in_img), 3);
% sz = min(size(in_img));
% in_img = in_img(1:sz, 1:sz);
% in_img_resize = imresize(in_img, [size(v.X, 1), size(v.X, 2)]);
% out_filename = ['TeddySmiling', num2str(length(v.e))];

%% define cortex

c.cortexSize = [60,100]; % degrees top to bottom, degrees LR, divide by 2 to get the actual mm that are useful
c.cortexCenter = [30,0];
c.pixpermm = 8; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased
fov = 30;
v.retinaSize = [fov, fov]*2.5; v.pixperdeg = 20;  %visual field map size and samping
c = p2p_c.define_cortex(c); % define the properties of the cortical map

%% set up the electrode locations

xList = linspace(0,fov, nelect/2);
yList = linspace(-fov,fov,nelect);

[X, Y] = meshgrid(xList, yList) ; X = X(:); Y = Y(:);
for i = 1:length(X)
    v.e(i).x = X(i);
    v.e(i).y = Y(i);
    [a, e]= cart2pol(v.e(i).x, v.e(i).y);
    v.e(i).ang = a*180/pi;
    v.e(i).ecc = e;
    c.e(i).radius = .001;% units are cm
end

%% define electrodes
c = p2p_c.define_electrodes(c, v); % defines properties for each electrode in retinal space
v = p2p_c.c2v_define_electrodes(c,v); % convert electrode locations from cortex to visual space

%% set up model

c = p2p_c.define_electrodes(c, v); % defines properties for each electrode in retinal space
v.eccList = [1 2 4 8 12];
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
c = p2p_c.generate_ef(c); % generate map of the electric field on cortical surface

%% display electric fields on cortex
ef = c.e(1).ef;
for e=2:length(c.e)
    ef = ef + c.e(e).ef;
end
figure(1); clf
p2p_c.plotcortgrid(256.*(ef./max(ef(:))), c, gray(256), 1, ['title(''electric field'')']); % shows electrodes on the cortical surface

%% load in image
in_filename = 'FruitNinjaIone.jpg';
in_img = imread(in_filename); in_img = mean(double(in_img), 3);
big_img = 127+zeros(max(size(in_img)));
in_img = insertImg(big_img, in_img);
in_img_resize = imresize(in_img, [size(v.X, 1), size(v.X, 2)]);

%% create rf maps

% generates the sum of weighted receptive fields activated by an electrode
% normalized so the max is 1
frameImg = zeros(size(v.X));frameImgNorm = frameImg;
idx = 1:length(v.e); ne = length(v.e);
for ii = 1:ne
    disp([num2str(round((100*ii)/length(idx))),  '% electrodes complete' ]);
    rfmap1 = zeros([size(v.X)]); rfmap2 = rfmap1; amp1 = 0; amp2 = 0;
    
    for pixNum = 1:length(c.X(:))
        if c.e(idx(ii)).ef(pixNum) > c.efthr
            x0 = c.v2c.X(pixNum); % x center
            y0 = c.v2c.Y(pixNum); % y center
            if strcmp(c.e(idx(ii)).hemi, 'lh')
                x0=-x0; y0 = -y0;
            end
            theta = pi-c.ORmap(pixNum);  %orientation
            sigma_x = c.RFmap(pixNum) * c.ar; % major axis sd
            sigma_y = c.RFmap(pixNum); % minor axis sd
            
            % things you need to go from sigmas and theta to G
            aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
            bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
            cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);
            
            % oriented 2D Gaussian
            
            G = c.e(idx(ii)).ef(pixNum) * exp( - (aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
            rfmap1 = rfmap1 + c.ODmap(pixNum)*G;
        end
    end
    if sum(rfmap1(:)>0)<20
        disp('WARNING! Too few pixels passed ef threshold.');
        disp(' try lowering c.efthr, checking location of electrodes relative to cortical sheet & ');
        disp('checking the sampling resolution of cortex');
    end
    ind = find(rfmap1>.5);
    amp1 = mean(in_img_resize(ind));
    frameImg = frameImg + rfmap1.*amp1;
    frameImgNorm = frameImgNorm + rfmap1;
    
    % do the other side
    rfmap2 = fliplr(rfmap1);
    ind = find(rfmap2>.5);
    amp2 = mean(in_img_resize(ind));
    frameImg = frameImg + rfmap2.*amp2;
    frameImgNorm = frameImgNorm + rfmap2;
    figure(5); subplot(1, 2,1)
    imagesc(rfmap1); axis off; axis equal;
    subplot(1, 2,2)
    imagesc(rfmap2); axis off; axis equal;colormap(gray); drawnow
end

figure(5);
imagesc(frameImg./frameImgNorm); axis off; axis equal; colormap(gray); drawnow

function img = insertImg(img,testImg)
% insertImg(img,testImg)
%
% Inserts testImg into the center of img.
% if testImg is larger than img, testImg is cropped and centered.

if size(testImg,1)>size(img,1)pwd
    x0 = ceil((size(testImg,1)-size(img,1))/2)+1;
    testImg = testImg(x0:(x0+size(img,1)-1),:);
end

if size(testImg,2)>size(img,2)
    y0 = ceil((size(testImg,2)-size(img,2))/2)+1;
    testImg = testImg(:,y0:(y0+size(img,2)-1),:);
end

x0 = ceil((size(img,2)-size(testImg,2))/2)+1;
y0 = ceil((size(img,1)-size(testImg,1))/2)+1;
img(y0:(y0+size(testImg,1)-1),x0:(x0+size(testImg,2)-1)) = testImg;
end

