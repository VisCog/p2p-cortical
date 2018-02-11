function rootPath=ebsRootPath()
%
%        rootPath =ebsRootPath;
%
% Determine path to root of the ebs directory
%
% This function MUST reside in the directory at the base of the ebs directory structure
%

rootPath=which('ebsRootPath');

rootPath=fileparts(rootPath);

return
