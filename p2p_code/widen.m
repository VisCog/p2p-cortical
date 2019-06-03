function widen(p)
% widen(p)
%
% Speads x-axis limits by a proportion p > 1. Designed to pull extreme data
% points off of the horizontal axes.
% 
% Input:
%   p          Proportion that scales the current figure's x-axis, numeric
%              (default: 1.1)

% Written GMB 2/27/97 - Wrote the three lines of code.
% Edited by Kelly Chang - March 15, 2017

%% Input Control

if ~exist('p', 'var')
	p = 1.1;
end

if p < 1
    error('proportion ''p'' must be greater than 1');
end

%% Widen X-Axis

xlim = get(gca,'XLim');
set(gca, 'XLim', ([-0.5 0.5]*diff(xlim)*p)+mean(xlim));