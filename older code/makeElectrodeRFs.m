function ERF = makeElectrodeRFs(xe,ye,sd,x,y)

ne = length(xe(:));
npix = length(x(:));

ERF = zeros(npix,ne);
for i=1:length(xe(:))
    G = exp(-( (x-xe(i)).^2 + (y-ye(i)).^2)/(2*sd^2));
    ERF(:,i) = G(:);
end

    



