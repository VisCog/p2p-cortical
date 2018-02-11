function ERF = makeElectrodeRFs(xe,ye,sd,x,y)

ne = length(xe(:));
npix = length(x(:));

ERF = zeros(npix,ne);
for i=1:length(xe(:))
    R=sqrt((x-xe(i)).^2+(y-ye(i)).^2);
    cspread=ones(size(x));
    cspread(R>sd)=2/pi*(asin(sd./R(R>sd)));
   % G = exp(-( (x-xe(i)).^2 + (y-ye(i)).^2)/(2*sd^2));
   % ERF(:,i) = G(:);
   ERF(:,i) = cspread(:);
end

    


% 
% %make cspread
% [x,y] = meshgrid(-1000:25:1000);  %microns
% rad = sqrt(x.^2+y.^2);
% cspread = ones(size(x));
% cspread(rad>STIM.electrodeRAD) = 2/pi*(asin(STIM.electrodeRAD./rad(rad>STIM.electrodeRAD)));