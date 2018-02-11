function [nearest, inds] = findnearest(a, b)
% find the nearest value in b for each element in a. probably some built-in
% matlab function already does this. but who cares...
%
%   [nearest, inds] = findnearest(a, b)
%
% Example:
%   a = rand(1,100);
%   b = 0:.001:1;
%   [nearest inds] = findnearest(a, b)

% 
% % let's make A a row vector and B a column vector
% a = a(:)';
% b = b(:);
% 
% [~, inds] = min(abs(bsxfun(@minus, b, a)));
% 
% 
% nearest = b(inds);

nearest = interp1(b,b,a,'nearest');
inds    = interp1(b,1:length(b),a,'nearest');


