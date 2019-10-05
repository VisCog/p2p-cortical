function [params,err] = UWfit(funName, params, freeList, varargin)
% [params,err] = UWfit(funName, params, freeList, var1, var2, var3,...)
%
% Helpful interface to matlab's 'fminsearch' function.
%
% Inputs:
%   funName        Function to be optimized. Must have form
%                  [err] = <funName>(params, var1, var2, ...)
%
%   params         A structure of with field names with corresponding
%                  parameter values for fitted function. 'freeList'
%                  parameters must only have a singular value per parameter
%       options    A structure with pptions for fminsearch program
%                  (see OPTIMSET)
%
%   freeList       Cell array containing list of parameter names (strings)
%                  to be free in fitting
%
%   var<n>         Extra variables to be sent into fitted function
%
% Outputs:
%   params         A structure with best fitting parameters as fieldnames
%
%   err            Error value at minimum, numeric
%
% Notes:
% - Dependencies: params2var.m, var2params.m, fitFunction.m

% Written by Geoffrey M. Boynton, Summer of '00
% Edited by Kelly Chang, February 10, 2017

%% Input Control

options = []; % options for fminsearch (see OPTIMSET)
if isfield(params, 'options')
    options = params.options;
end

if isempty(freeList)
    freeList = fieldnames(params);
end

%% Fit Function and Calculate Final Error

% turn free parameters in to 'var'
vars = params2var(params, freeList);

% calling fminsearch
vars = fminsearch('fitFunction', vars, options, funName, params, freeList, varargin);

% assign final parameters into 'params'
params = var2params(vars, params, freeList);

% evaluate the function 'funName' for error at minimum
err = fitFunction(vars, funName, params, freeList, varargin);