function f = constructprfmodel2Dgaussian(res)

% function f = constructprfmodel2Dgaussian(res)
%
% <res> is the number of pixels along each side of the image
%
% return a PRF model of the form {A B C} suitable for use with fitprf.m.
% the model is an isotropic 2D Gaussian function parameterized by three
% parameters:
%   x is the column index of the peak of the Gaussian
%   y is the row index of the peak of the Gaussian
%   sigma is the standard deviation of the Gaussian (in matrix units)
% the Gaussian function is unit-length-normalized and reshaped into
% a column vector.  (notice that gain flexibility is not included
% in the model.)
%
% technical details:
% - we use a fairly uninformative initial seed of [(1+res)/2 (1+res)/2 res].
% - we enforce lower and upper bounds of [1-res 1-res 0;
%                                         2*res 2*res Inf].
%
% example:
% f = constructprfmodel2Dgaussian(100);
% model = feval(f{3},[30 10 5]);
% figure; imagesc(reshape(model,[100 100])); axis equal tight;

% pre-compute xx and yy for speeding up execution
[d,xx,yy] = makegaussian2d(res,[],[],1,1);

% do it
seed = [(1+res)/2 (1+res)/2 res];
bounds = [1-res 1-res 0; 2*res 2*res Inf];
f = {seed bounds @(params) vflatten(unitlength(makegaussian2d(res,params(2),params(1),params(3),params(3),xx,yy,0),[],[],0))};
