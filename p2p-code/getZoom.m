function zoom = getZoom(sites)

zoom = zeros(length(sites),4);

for ii = 1:length(sites)
   switch sites(ii)
       case 1, zoom(ii,:) = [5 35 -5 25];
       case 2, zoom(ii,:) = [-16 -4 -6 6];
       case 3, zoom(ii,:) = [-8 0 -1 7];
       case 4, zoom(ii,:) = [-3.5 1.5 -2 3];
       case 5, zoom(ii,:) = [-2 1 -1.5 1.5];
    
   end
end
 