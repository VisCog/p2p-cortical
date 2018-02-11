function [colors, sites] = getColors(sites)
% RFB colors for plotting 5 sites
%
% Example: 
%   colors = getColors(1)
%   plot(rand(1,5), rand(1,5), 'Color', colors);

all_sites = 1:5;

all_colors = [...
    0.5000         0         0    
    0.5000    1.0000    0.5000;
    1.0000    0.8125         0
         0    0.8750    1.0000
         0         0    1.0000 %        0.5625
];
    
if ~exist('sites', 'var'),  sites = all_sites; end

inds = zeros(size(sites));

for ii = 1:length(sites)
    inds(ii) = find(all_sites == sites(ii));
end
colors = all_colors(inds,:);

return 