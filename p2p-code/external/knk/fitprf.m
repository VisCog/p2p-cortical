function [params,paramsse,r,rrun,polyparams,polymeans,numiters,hrf,betas,signal,drift] = ...
  fitprf(stimulus,response,prfmodel,hrfmodel,flag,maxpolydeg,ar,mode,maxiter, ...
         wantresample,tol,extraopt,extraregressors,hrfnormfun,derivemode,metric,outputfcn)

% function f = fitprf(...)
%
% by specifying a single output argument, you will receive in this argument
% a struct that encapsulates all the usual outputs as detailed below.
%
%   OR
%
% function [params,paramsse,r,rrun,polyparams,polymeans,numiters,hrf,betas,signal,drift] = ...
%   fitprf(stimulus,response,prfmodel,hrfmodel,flag,maxpolydeg,ar,mode,maxiter, ...
%          wantresample,tol,extraopt,extraregressors,hrfnormfun,derivemode,metric,outputfcn)
%
% <stimulus> is the stimulus with dimensions time x stim.
%   can be a cell vector that indicates different runs.  different runs can 
%   have different numbers of time points.  <stimulus> can also be a function
%   that returns the stimulus when evaluated.  
% <response> is the data with dimensions time x 1.  can be a cell vector and 
%   should mirror <stimulus>.  all values should be finite and double-format.
% <prfmodel> is
%   {A B C} where
%     A is the initial seed (1 x p)
%     B are the bounds (2 x p)
%     C is a function handle that transforms parameters (1 x p) into
%       a PRF (stim x 1).  note that allowing gain flexibility here
%       is not necessary since gain is explicitly handled 
%       (see evaluateprfmodel.m for details).
%  OR
%   a matrix of dimensions stim x n with different basis functions
%   along the second dimension.  n must be >= 1.
%  OR
%   0 is a special case that means act as if we have delta basis functions,
%   e.g., eye(stim).  the point of this is to enable us to speed-up the code.
% <hrfmodel> is
%   {A B C D} where
%     A is the initial seed (1 x q)
%     B are the bounds (2 x q)
%     C is a function handle that transforms parameters (1 x q) into
%       an HRF (t x 1).
%     D is whether the model includes gain flexibility.  can be omitted
%       or [], in which case we default to 0.
%  OR
%   a matrix of dimensions t x b with different basis functions
%   along the second dimension.  b must be >= 1.
% <flag> takes effect in the case that PRF and HRF are both basis functions.
%   0 or {0} means do not enforce separability between the PRF and HRF, so we can go
%     ahead and use ordinary least-squares (OLS) estimation (which is fast!).  (this assumes
%     that <mode> is 0 or {0}; if not, then we have to resort to the usual lsqnonlin.m
%     optimization.)
%   1 or {1} means enforce separability between the PRF and HRF.  this forces us to
%     use the usual slow lsqnonlin.m optimization.  the main advantage of separability
%     is that it may vastly reduce the number of free parameters in the model.  of
%     course, the benefit of enforcing separability depends on whether that assumption
%     is reasonable for your specific paradigm.
%   {2 X Y} means do not enforce separability between the PRF and HRF, and go ahead and
%     use regularized regression as implemented in gradientdescent_wrapper.m.  (this
%     assumes that <mode> is 0 or {0}; if not, then we have to resort to the usual
%     lsqnonlin.m optimization).  X and Y are the <fittype> and <randomtype> inputs 
%     to gradientdescent_wrapper.m.  we use defaults for the other inputs.
%     note that it is tricky to combine resampling with regularized regression 
%     due to the issue that the number of iterations can change a lot across resamples,
%     causing large variance in the magnitude of the parameter estimates.  note that
%     we ignore the DC estimate returned by gradientdescent_wrapper.m; since we already
%     project out the polynomial nuisance functions (when <maxpolydeg> is finite) or
%     omit polynomials (when <maxpolydeg> is NaN), we are in effect already 
%     assuming that DC is 0.
% <maxpolydeg> (optional) is a non-negative integer with the maximum polynomial degree
%   to use for polynomial nuisance functions.  perhaps use round(L/2) where L is the
%   number of minutes long your run is, but no less than 2?  special case is NaN which
%   means do not use any polynomials at all.  default: 1.
% <ar> (optional) is no longer used.  omit or pass in [].
% <mode> (optional) allows for a divisive-normalization nonlinearity to be applied
%     after <stimulus> is projected onto <prfmodel> (see divnorm.m).  this is probably
%     useful in only certain special cases.  note that if <prfmodel> is the basis 
%     function case, the nonlinearity is applied to the projection of <stimulus> 
%     onto each basis function (and this occurs before any weighting of the various
%     basis functions).
%   0 or {0} means normal
%   {1 A B} means allow exponent.  A is the initial seed for the exponent (1 x 1).
%     B are the bounds (2 x 1).
%   {2 A B} means allow full divisive normalization.  A is the initial seed for
%     the exponent and normalization term (1 x 2).  B are the bounds (2 x 2).
%   when <mode> is specified (i.e. not 0 nor {0}), then we perform stepwise fitting
%   in which we optimize all parameters except the nonlinear parameters, with the
%   nonlinear parameters fixed at the initial seed as specified by <mode>.  then,
%   we optimize all parameters simultaneously.  the point of this is to try to
%   avoid local minima problems.  the <numiters> output reflects the number of
%   iterations in the second (and final) optimization.
%   default: 0.
% <maxiter> (optional) is the maximum number of iterations.  this input is used
%   only in the lsqnonlin.m case.  default: Inf.  can also be [A B C].  
%   in this case, we first try to fit with a max number of iterations of A.  
%   then we check if we stopped because we reached the max number of iterations; 
%   if so, then if <r> >= B, then we re-do the whole process with a max number 
%   of iterations of C.  (if there are multiple resamplings, we check if we 
%   stopped for any of them because of reaching the max number of iterations.)
% <wantresample> (optional) is
%   0 means to not bootstrap nor cross-validate.
%   N or {N SEED GROUP} means number of bootstraps to take, the random number seed to use,
%       and the grouping to use (a vector of positive integers).  for example, a grouping
%       of [1 1 1 2 2 2] means to draw six bootstrap samples in total, with three
%       bootstrap samples using only the first three cases and three bootstrap samples
%       using only the second three cases.
%     positive values for N mean to bootstrap individual data points.  this is not
%       recommended for various reasons (see note below).
%     negative values for N mean to bootstrap entire runs (thus, <stimulus> should be a
%       cell vector with multiple elements).
%     if the seed is [] or not provided, we use sum(100*clock).
%     if the grouping is [] or not provided, we use ones(1,D) where D is the
%       number of data points or runs.
%   [A B C ...] specifies the cross-validation scheme.  this should be a vector of
%     the same length as <stimulus> in the cell vector case, and should have at least
%     two elements.  elements of the vector should be non-negative integers 
%     and should cover either the range 0 to J for some positive J or the range 1 to J
%     for some positive J.  on the first cross-validation run, we train on all data not
%     labeled 1 and then validate on data labeled 1.  on the second cross-validation
%     run, we train on all data not labeled 2 and then validate on data labeled
%     2.  and so forth.  you can also specify 'n-fold' for <wantresample>; this is
%     a shortcut that means [1 2 3 ...].
%  -M specifies the cross-validation scheme.  (this scheme subsumes the previous scheme
%     and provides additional functionality.)  M should be a matrix with dimensions
%     A x B where A >= 1 indicates different cross-validation cases and B >= 2 matches
%     the number of elements in the cell vector case of <stimulus>.  elements of M should
%     be consecutive integers, possibly skipping 0.  positive integers indicate runs to use 
%     for training; negative integers indicate runs to use for testing.  runs that share the 
%     same value are grouped together, which affects the derivation of beta weights 
%     (see <derivemode>).  every cross-validation case must involve at least one training 
%     run (but can include zero testing runs).  for example, -[1 1 1 1 -1 -1 -2 -2 0] 
%     specifies a cross-validation scheme in which one cross-validation case is performed.  
%     in this cross-validation case, the first four runs are used for training and then the 
%     5th through 8th runs are used for testing.
%   default: 0.
% <tol> (optional) is the tolerance to use for TolFun and TolX in the lsqnonlin.m case.
%   default: 1e-10.  can be [A B] where A is the tolerance on TolFun and B is the tolerance on TolX.
% <extraopt> (optional) is a cell vector of extra parameter/value pairs to use
%   in the lsqnonlin.m optimization options (e.g. {'Display' 'final'}).  default: {}.
% <extraregressors> (optional) is one or more extra regressors with dimensions time x r.
%   can be a cell vector and should mirror <stimulus>.  the number of extra regressors
%   does not have to be the same across different runs.  we orthogonalize the 
%   <extraregressors> with respect to the polynomials up front.  <extraregressors> can 
%   also be a function that returns the extra regressors when evaluated.  default: [].
% <hrfnormfun> (optional) is a function that takes an HRF (time x 1) and outputs a scalar
%   that the HRF should be divided by in order to normalize the HRF.  <hrfnormfun> is 
%   useful for applying a scale normalization to the HRF before proceeding with derivation
%   of beta weights (see <derivemode>).  please note that in the inseparable case 
%   (see <flag>), we consider only the first PRF component's HRF.  default: @alwaysone.
% <derivemode> (optional) indicates the way in which beta weights are to be derived.  derivation of
%   beta weights occurs only when <prfmodel> is the basis function case and <mode> is 0 or {0}.
%   if these conditions are not met, then the derivation of beta weights is skipped.  when 
%   <wantresample> is 0 or N or {N SEED GROUP}, <derivemode> is either 0 which means skip the 
%   beta weight derivation or 1 which means perform the derivation.  when <wantresample> is 
%   [A B C ...] or -M, the format of <derivemode> is [A B] where A and B are
%     0 means include no beta weights
%     1 means include beta weights from fitting all groups of runs at once
%     2 means include beta weights from fitting individual groups of runs
%     3 means include both types of beta weights
%   A controls beta weights from the training data; B controls beta weights from the testing data.
%   default is [] which means to skip the beta weight derivation stuff.
% <metric> (optional) is a function that accepts two column vectors (the first is the model,
%   the second is the data) and outputs a value.  default is @calccod, which computes the 
%   percent variance of the data explained by the model.  note that polynomial nuisance functions
%   are projected out from both the data and the model before computation of the <metric>
%   (for more details, see below), except, of course, when <maxpolydeg> is NaN.
% <outputfcn> (optional) is a function suitable for use as an 'OutputFcn'.
%   <outputfcn> matters only when the lsqnonlin.m optimization is used.
%   default: @alwayszero.
%
% this function fits a general class of models suitable for application to fMRI time-series.
% we generally refer to the models as PRF models (population receptive field models) only
% for the reason that one type of model that can be fit is the PRF models described in
% Dumoulin and Wandell, 2008.  but more generally, we note that the kinds of models that can
% be fit by this function include what is commonly referred to as GLM models in the fMRI
% literature.  some functionality that this function provides that may be unique include
% allowing for flexibility in the shape of the HRF (hemodynamic response function), enforcing
% separability between the HRF and the PRF, allowing for bootstrapping entire time-series
% runs (in order to obtain good estimates of standard error), and allowing for cross-validation
% (in order to obtain unbiased measures of model accuracy).  in particular, bootstrapping
% and cross-validation are essential to those interested in exploring and comparing different
% models of the data.
%
% both the PRF and HRF can be modeled using two different approaches.  one approach is to
% use a parametric function to generate the PRF or HRF.  for example, a gamma function
% might be used to generate the HRF.  the advantage of this approach is that the parameters
% of the parametric function might be easily interpretable.  also, parametric functions
% may help constrain the space of potential solutions to ones that are likely to be accurate,
% while only expending a few free parameters.  disadvantages of the parametric approach
% include the fact that the parametric function may be too restrictive (i.e. the space of
% solutions may be too restrictive), the fact that fitting the parameters of a parametric
% function in general requires nonlinear fitting techniques (e.g. MATLAB's lsqnonlin.m),
% and the fact that fitting the parameters may suffer from local minima issues.
% the second approach is to use a set of basis functions to generate the PRF or HRF (whereby
% the final model consists of a linear combination of the basis functions).  for example, 
% a set of sinusoids might be used to generate the HRF.  the main advantage of this approach
% is that fitting the parameters becomes relatively easy (i.e. in certain cases, linear
% fitting techniques can be used; the potential for local minima issues becomes less likely).
% the disadvantage of this approach is that a large number of basis functions may be 
% necessary to achieve the desired amount of flexibility.  in the end, which specific
% combination of PRF and HRF models is best can be decided objectively through the use
% of cross-validation.
%
% this function leaves it up to the caller of the function to specify the specific types 
% of PRF and HRF models to be used.  this generality provides an immense amount of power 
% and flexibility, but potentially increases complexity of use.
%
% the basic form of the kinds of models fitted by this function is given by:
%   conv(divnorm(stimulus*prfmodel) * prfweights?,hrfmodel * hrfweights?).
% this can be decomposed as follows:
%   1. stimulus*prfmodel => project stimulus onto prfmodel
%   2. divnorm(...) => apply divisive-normalization nonlinearity (if applicable)
%   3. ... * prfweights? => perform weighted sum if there are multiple prfmodel basis functions
%   4. hrfmodel * hrfweights? => perform weighted sum if there are multiple hrfmodel basis functions
%   5. conv(...,...) => convolve the prf component with the hrf component
% also, the PRF model includes a set of polynomials of increasing degree (e.g. 0, 1, and 2)
% for each individual run (except when <maxpolydeg> is NaN).  these polynomials allow 
% the capturing of low-frequency drift which is common in fMRI time-series data.
%
% return:
%  <params> is 1 x p with the parameters (see below for details).
%    when <wantresample>, we return <params> as b x p where the rows correspond to different resamples.
%  <paramsse> is 1 x p with the standard deviation across resamples.
%    NaNs are returned if isequal(<wantresample>,0).  be careful in the interpretation of <paramsse>, as
%    certain variants of <params> may have scaling issues that should be dealt with before calculating
%    the standard deviation across resamples.  we leave it up to the user to do this.
%  <r> is the <metric> between <modelfit> and <responseB>.
%    when <wantresample> is 0 or N or {N SEED GROUP}, this is straightforward:
%     <modelfit> is time x 1 with the model fit after projecting out both the polynomial nuisance functions
%       and the <extraregressors>.  if <wantresample> is N or {N SEED GROUP}, we calculate the mean model 
%       fit across bootstraps.
%     <responseB> is time x 1 with the data after projecting out both the polynomial nuisance functions
%       and the <extraregressors>.
%    when <wantresample> is [A B C ...] or -M, this is a bit tricky:
%     <modelfit> is timeV x 1 with the cross-validated model fits after projecting out the polynomial
%       nuisance functions.
%     <responseB> is timeV x 1 with the cross-validation data concatenated together.  these data have had
%       the polynomial nuisance functions projected out.
%     note that if there are no testing runs, then it is possible for <modelfit> and <responseB> to be
%       both empty.  in this case, we return <r> as NaN.
%  <rrun> is a vector of the correlations between each pair of <modelfitsep> and <responseBsep>.
%    <modelfitsep>,<responseBsep> is just like <modelfit> and <responseB> except that the data are
%      segregated by runs and embedded in a cell vector.  note that it is possible for <modelfitsep>
%      and <responseBsep> to be both empty like {}.  in this case, we return <rrun> as NaN.
%  <polyparams> is length(stimulus) x <maxpolydeg>+1 with the weights on the polynomial nuisance functions.
%    note that there is a separate set of weights for each stimulus run.  note that the cross-validation
%    cases of <wantresample> presents a tricky situation.  we handle this situation by overwriting any previously
%    determined sets of weights.  for example, suppose <wantresample> is [0 1 1 2 2].  on the first cross-
%    validation run, we fill in rows 1, 4, and 5 of <polyparams>.  on the second cross-validation run,
%    we fill in rows 1, 2, and 3 of <polyparams> (overwriting the previous entries for the first row).
%    any rows that are not filled in are left being NaNs.  when <wantresample> is N or {N SEED GROUP} 
%    where abs(N) >= 1, then we return <polyparams> as length(stimulus) x <maxpolydeg>+1 x abs(N).
%    note that when <maxpolydeg> is NaN, <polyparams> will be returned as [].
%  <polymeans> is 1 x length(stimulus) with the mean of the time-series of the fitted nuisance functions.
%    for the cross-validation cases of <wantresample>, we handle this in a similar way as we do for <polyparams>.
%    when <wantresample> is N or {N SEED GROUP} where abs(N) >= 1, then we return <polymeans> as abs(N) x length(stimulus).
%    note that when <maxpolydeg> is NaN, <polymeans> will be returned as [].
%  <numiters> is 1 x 1 with the number of iterations taken.  when <wantresample>, we return <numiters>
%    as 1 x b where each element corresponds to a different resampling.  in the special OLS case,
%    <numiters> is returned as [].  in the special regularized-regression case, <numiters> is returned
%    with the <numiters> output from gradientdescent_wrapper.m.
%  <hrf> is time x 1 is the estimated HRF after dividing by the output of <hrfnormfun>.  when <wantresample>, we
%    return <hrf> as time x b where the columns correspond to different resamples.
%  <betas> is n x z with derived beta weights.  this output is assigned only when specific conditions are met
%    (see <derivemode>).  if these conditions are not met, <betas> is returned as [].  the number of rows in
%    <betas> is equal to the number of basis functions in <prfmodel> (with one rare exception, see below).
%    when <wantresample> is 0, the number of columns is 1.  when <wantresample> is N or {N SEED GROUP}, 
%    the number of columns is the number of bootstraps.  when <wantresample> is [A B C ...] or -M, the 
%    number of columns depends on the number of cross-validation cases and the specific value of <derivemode>.  
%    specifically, columns are ordered by cross-validation case and then, within each cross-validation case, 
%    by training beta weights and testing beta weights.
%  <signal> is time x 1 with the model fit.  when <wantresample> is N or {N SEED GROUP}, we concatenate model fits from 
%    different bootstraps.  when <wantresample> is [A B C ...] or -M, we calculate model fits for only
%    the training runs and concatenate model fits from different resamplings.  
%  <drift> is time x 1 with the time-series of the fitted nuisance functions in a format that is 
%    analogous to <signal>.
%
% history:
% 2011/10/09 - add input <outputfcn>.  also, we automatically enforce a sanity check via outputfcnsanitycheck.m.
% 2011/09/21 - allow <extraregressors> to be a function
% 2011/08/30 - now use stepwise fitting in the nonlinear case (i.e. <mode> is not 0 nor {0}).
% 2011/08/04 - add a speed up
% 2011/06/27 - previously, we would crash when the data had no variance (e.g. all zeros).  now we don't crash.
% 2011/04/09 - implement NaN case for <maxpolydeg>
% 2010/12/06 - introduce nanreplace into the lsqnonlin objective to prevent crazy parametric model return values from causing chaos.
% 2010/12/06 - implement encapsulation (single output case); implement 'n-fold' for <wantresample>; change initial seeds for basis function cases (in particular, we no longer force weights to 1.  it seems that nonlinear fitting can operate successfully even if there is a gain ambiguity); fix minor bug (it would have crashed)
% 2010/12/03 - change the form of the input <hrfnormfun> (this included changes to <betas>)!  allow <derivemode> to operate in all resampling circumstances.
% 2010/12/02 - implement {N SEED GROUP} for <wantresample>; fix bug --- in the re-do case the signal and drift were not getting assigned;
%              add the <metric> input.
% 2010/12/01 - the <ar> input has been deprecated (but still works); the 
%              <meanint>,<driftstd>,<signalrms>,<noiserms> outputs have been removed.
% 2010/11/15 - allow derivation of beta weights even in the inseparable case now.  (in evaluateprfmodel, 
%              hrfactual is now returned for the internal hrfactual variable.)
% 2010/10/22 - fix bug concerning exitflag and maxiters.  if there were multiple resamplings, then
%              we were checking only the last exitflag when deciding to re-do the fitting process.
%              now, we do the right thing --- we check to see if any of the resamplings stopped because
%              they reached the max number of iterations.
% 2010/10/22 - change initial seed for prfmodel and hrfmodel basis functions to 0!
% 2010/10/21 - implement gain flexibility for hrfmodel
% 2010/10/09 - add meanint,driftstd,signalrms,noiserms outputs.
% 2010/10/06 - add signal and drift outputs.
% 2010/09/08 - major update: implement bootstrapping of runs; implement -M case for <wantresample>; add
%              new input <hrfnormfun> and return new output <hrf>; add new input <derivemode> and 
%              return new output <betas>.
% 2010/07/17 - add input <extraregressors>
% 2010/07/09 - use sparse matrices for speedup
% 2010/06/15 - use calccod.m instead of calccorrelation.m.  also, remove rgain and rrungain output arguments.
%
% the meaning of the output <params> depends on the specific combination
% of <prfmodel> and <hrfmodel> used:
%
%   case = PRF: parametric, HRF: parametric, no gain in HRF
%     [a1 a2 ... b1 b2 ... c n sigma] where
%     a1,a2,... are parameters for the PRF model.
%     b1,b2,... are parameters for the HRF model.
%     c is the final gain.
%     n and sigma are parameters for divisive-normalization nonlinearity.
%
%   case = PRF: parametric, HRF: parametric, gain in HRF
%     [a1 a2 ... b1 b2 ... n sigma] where
%     a1,a2,... are parameters for the PRF model.
%     b1,b2,... are parameters for the HRF model.
%     n and sigma are parameters for divisive-normalization nonlinearity.
%
%   case = PRF: parametric, HRF: basis
%     [a1 a2 ... b1 b2 ... n sigma] where
%     a1,a2,... are parameters for the PRF model.
%     b1,b2,... are weights on the HRF basis functions.
%     n and sigma are parameters for divisive-normalization nonlinearity.
%
%   case = PRF: basis, HRF: parametric, no gain in HRF
%     [a1 a2 ... b1 b2 ... n sigma] where
%     a1,a2,... are weights on the PRF basis functions.
%     b1,b2,... are parameters for the HRF model.
%     n and sigma are parameters for divisive-normalization nonlinearity.
%
%   case = PRF: basis, HRF: parametric, gain in HRF
%     [a1 a2 ... b1 b2 ... n sigma] where
%     a1,a2,... are weights on the PRF basis functions.
%     b1,b2,... are parameters for the HRF model.
%     n and sigma are parameters for divisive-normalization nonlinearity.
%
%   case = PRF: basis, HRF: basis, separable
%     [a1 a2 ... b1 b2 ... n sigma] where
%     a1,a2,... are weights on the PRF basis functions.
%     b1,b2,... are weights on the HRF basis functions.
%     n and sigma are parameters for divisive-normalization nonlinearity.
%
%   case = PRF: basis, HRF: basis, inseparable
%     [1b1 1b2 ... 2b1 2b2 ... n sigma] where
%     1b1,1b2,... are the PRF weights for the first HRF basis function,
%     2b1,2b2,... are the PRF weights for the second HRF basis function,
%       and so forth.
%     n and sigma are parameters for divisive-normalization nonlinearity.
%
% a note on bootstrapping data points:
%   we implement bootstrapping by controlling for data points at the very last possible
% moment --- for example, in the call to lsqnonlin.m.  strictly speaking, this does
% not ensure perfect resampling / independence between data points, since we do things
% like projecting out the polynomial nuisance functions based on all the data.  moreover,
% noise may be correlated from time point to time point, making the bootstrap samples 
% non-independent.  of course, you can go all the way back to the pre-processing and 
% claim that no consecutive time points can ever be independent.  it is not clear how 
% much of a stickler we need to be.  to be completely independent, the best strategy 
% is to resample entire runs, which is the other option.
%
% a note on beta weight derivation:
% - in the case where <wantresample> is 0 or N or {N GROUP SEED}, beta weight derivation
%   is easy --- the only thing we need to do is to take the parameters in <params> and
%   multiply them by the output of <hrfnormfun>.  (however, there is one tricky rare case:
%   in the inseparable case, there is the potential to have multiple PRF components that each
%   have their own HRF.  how we handle this case is to use the first PRF component's HRF 
%   as our HRF, divide it by the output of <hrfnormfun>, and return in <betas> the adjusted
%   beta weight for ONLY the first PRF component.  the reason is that it is not straightforward
%   to compare different PRF components that have different HRFs.)
% - in the case where <wantresample> is [A B C ...] or -M, things are a bit more tricky,
%   since the user might want to see derived beta weights for individual runs (and this
%   is not available from the main model fit).  in this case, each column of <betas> 
%   represents a set of beta weights obtained using OLS regression.  The OLS regression 
%   uses the estimated HRF after dividing by the output of <hrfnormfun> and attempts to 
%   determine beta weights corresponding to the various basis functions in <prfmodel>.  
%   The regression includes one regressor for each <prfmodel> basis function as well as 
%   regressors for the polynomial nuisance functions and the <extraregressors>.  however, 
%   only the estimated weights for the <prfmodel> basis functions are recorded (and are 
%   returned in <betas>); the estimated weights for the polynomials and <extraregressors> 
%   are not recorded.  
%
% other notes:
% - in the OLS case, different time-series could be fit in one single matrix multiplication.
%   unfortunately, given how the code is designed, this is actually very difficult to implement.
%   so, the current mode of operation is to process each time-series separately, which is
%   fine, except that it is relatively slow.
% - we use certain initial seeds for the basis function cases.  hopefully, these initial seeds
%   are generic and non-informative and will work in all cases.

% TODO:
% - to save space at the fitprfmulti level, we could postpone loading of extraregressors.

% input
if ~exist('maxpolydeg','var') || isempty(maxpolydeg)
  maxpolydeg = 1;
end
if ~exist('ar','var') || isempty(ar)
  ar = 0;
end
if ~exist('mode','var') || isempty(mode)
  mode = 0;
end
if ~exist('maxiter','var') || isempty(maxiter)
  maxiter = Inf;
end
if ~exist('wantresample','var') || isempty(wantresample)
  wantresample = 0;
end
if ~exist('tol','var') || isempty(tol)
  tol = 1e-10;
end
if ~exist('extraopt','var') || isempty(extraopt)
  extraopt = {};
end
if ~exist('extraregressors','var') || isempty(extraregressors)
  extraregressors = [];
end
if ~exist('hrfnormfun','var') || isempty(hrfnormfun)
  hrfnormfun = @alwaysone;
end
if ~exist('derivemode','var') || isempty(derivemode)
  derivemode = [];
end
if ~exist('metric','var') || isempty(metric)
  metric = @calccod;
end
if ~exist('outputfcn','var') || isempty(outputfcn)
  outputfcn = @alwayszero;
end
if isa(stimulus,'function_handle')
  stimulus = feval(stimulus);
end
if isa(extraregressors,'function_handle')
  extraregressors = feval(extraregressors);
end
if ~iscell(stimulus)
  stimulus = {stimulus};
end
if ~iscell(response)
  response = {response};
end
if ~iscell(flag)
  flag = {flag};
end
if ~iscell(mode)
  mode = {mode};
end
if length(tol) == 1
  tol = repmat(tol,[1 2]);
end
if ~iscell(extraregressors)
  extraregressors = {extraregressors};
end
if iscell(hrfmodel) && length(hrfmodel)==3
  hrfmodel{4} = [];
end
if iscell(hrfmodel) && isempty(hrfmodel{4})
  hrfmodel{4} = 0;
end

% calc stuff
stimulus = cellfun(@(x) full(x),stimulus,'UniformOutput',0);
numsets = length(stimulus);
time = cellfun(@(x) size(x,1),stimulus);  % number of time points
numstim = size(stimulus{1},2);  % number of stimuli
exitflag = NaN;  % just need to make sure it's defined

% figure out whether we want beta weights [TRICKY]
wantbetaweights = ...
  ~isempty(derivemode) && ...
  ~isequal(derivemode,0) && ...
  ~isequal(derivemode,[0 0]) && ...  % can't be the [0 0] case
  ~iscell(prfmodel) && ...  % must be the prfmodel basis function case
  mode{1}==0;  % must be the basic linear case

% convert wantresample case 'n-fold'
if isequal(wantresample,'n-fold')
  assert(numsets > 1,'you must have at least two runs to perform n-fold cross-validation');
  wantresample = 1:numsets;
end

% convert wantresample case [A B C ...] to -M format
wantresampleORIG = wantresample;  % keep a record of the original version due to the possibility of re-calling fitprf (see end of function)
if ~iscell(wantresample) && numel(wantresample) > 1 && all(wantresample(:) >= 0)
  wantresampleNEW = [];
  for p=1:max(wantresample)
    wantresampleNEW(p,:) = zeros(1,numsets);
    wantresampleNEW(p,wantresample~=p) = 1:count(wantresample~=p);
    wantresampleNEW(p,wantresample==p) = -(1:count(wantresample==p));
  end
  wantresample = -wantresampleNEW; clear wantresampleNEW;
end

% construct polynomial regressors matrix and a matrix that removes these polynomials from time-series
polyregressors = {}; polymatrix = {};
for p=1:numsets
  if isnan(maxpolydeg)
    polyregressors{p} = zeros(time(p),0);
    polymatrix{p} = eye(time(p));
  else
    polyregressors{p} = constructpolynomialmatrix(time(p),0:maxpolydeg);
    polymatrix{p} = projectionmatrix(polyregressors{p});
  end
end

% construct matrix that corrects for autocorrelated noise
armatrix = {};
for p=1:numsets
  armatrix{p} = toeplitz([1 -ar zeros(1,time(p)-2)]',[1 zeros(1,time(p)-1)]);
end

% construct matrix that removes the extra regressors
if length(extraregressors)==1
  extraregressors = repmat(extraregressors,[1 numsets]);
end
exmatrix = {};
for p=1:numsets
  if isempty(extraregressors{p})
    exmatrix{p} = eye(size(armatrix{p}));
  else
    extraregressors{p} = polymatrix{p}*extraregressors{p};  % orthogonalize wrt polynomials
    exmatrix{p} = projectionmatrix(extraregressors{p});
  end
end

% construct total transformation matrix
tmatrix = {};
for p=1:numsets
  tmatrix{p} = exmatrix{p}*armatrix{p}*polymatrix{p};
end

% define options
options = optimset('Display','iter','MaxFunEvals',Inf,'MaxIter',maxiter(1),'TolFun',tol(1),'TolX',tol(2),extraopt{:}, ...
                   'OutputFcn',@(a,b,c) outputfcnsanitycheck(a,b,c,tol(1),10) | feval(outputfcn,a,b,c));

% prep bootstrapping
if iscell(wantresample) || numel(wantresample)==1
  if iscell(wantresample)
    wantresample0 = wantresample{1};
    if length(wantresample) < 2 || isempty(wantresample{2})
      setrandstate;
    else
      setrandstate({wantresample{2}});
    end
    if length(wantresample) < 3 || isempty(wantresample{3})
      bgroup = [];
    else
      bgroup = wantresample{3};
    end
  else
    wantresample0 = wantresample;
    setrandstate;
    bgroup = [];
  end
  if wantresample0 > 0
    if isempty(bgroup)
      bgroup = ones(1,sum(time));
    end
    bootfun = {};
    for xx=1:wantresample0
      ix = [];
      for bbb=1:max(bgroup)
        ix = [ix subscript(find(bgroup==bbb),ceil(rand(1,count(bgroup==bbb))*count(bgroup==bbb)))];
      end
      bootfun{xx} = @(x) x(ix,:);  % each function picks a random subset of the time points
    end
  elseif wantresample0 < 0
    if isempty(bgroup)
      bgroup = ones(1,numsets);
    end
    bootfun = {};
    for xx=1:-wantresample0
      ix = [];
      for bbb=1:max(bgroup)
        ix = [ix subscript(find(bgroup==bbb),ceil(rand(1,count(bgroup==bbb))*count(bgroup==bbb)))];
      end
      bootfun{xx} = @(x) x(ix);  % each function picks a random subset of the runs
    end
  else
    bootfun = {@identity};  % this does nothing
  end
end

% deal with PRF and HRF seed and bounds (see documentation in evaluateprfmodel.m).
paramsA = [];  % extra parameters to tack on
if iscell(prfmodel) & iscell(hrfmodel) & hrfmodel{4}==0
  mainseed = [prfmodel{1} hrfmodel{1} 1];
  mainlb = [prfmodel{2}(1,:) hrfmodel{2}(1,:) -Inf];
  mainub = [prfmodel{2}(2,:) hrfmodel{2}(2,:) Inf];
elseif iscell(prfmodel) & iscell(hrfmodel) & hrfmodel{4}==1
  mainseed = [prfmodel{1} hrfmodel{1}];
  mainlb = [prfmodel{2}(1,:) hrfmodel{2}(1,:)];
  mainub = [prfmodel{2}(2,:) hrfmodel{2}(2,:)];
elseif iscell(prfmodel) & ~iscell(hrfmodel)
  mainseed = [prfmodel{1} zeros(1,size(hrfmodel,2))];
  mainlb = [prfmodel{2}(1,:) repmat(-Inf,[1 size(hrfmodel,2)])];
  mainub = [prfmodel{2}(2,:) repmat(Inf,[1 size(hrfmodel,2)])];
elseif ~iscell(prfmodel) & iscell(hrfmodel) & hrfmodel{4}==0  % this case used for spatial summation paper...
  temp = choose(isequal(prfmodel,0),zeros(1,numstim),zeros(1,size(prfmodel,2)));
  mainseed = [temp hrfmodel{1}];
  mainlb = [repmat(-Inf,[1 length(temp)]) hrfmodel{2}(1,:)];
  mainub = [repmat(Inf,[1 length(temp)]) hrfmodel{2}(2,:)];
elseif ~iscell(prfmodel) & iscell(hrfmodel) & hrfmodel{4}==1
  paramsA = [];
  temp = choose(isequal(prfmodel,0),zeros(1,numstim),zeros(1,size(prfmodel,2)));
  mainseed = [temp hrfmodel{1}];
  mainlb = [repmat(-Inf,[1 length(temp)]) hrfmodel{2}(1,:)];
  mainub = [repmat(Inf,[1 length(temp)]) hrfmodel{2}(2,:)];
elseif ~iscell(prfmodel) & ~iscell(hrfmodel) & flag{1}==1
  paramsA = [];
  temp = choose(isequal(prfmodel,0),ones(1,numstim),ones(1,size(prfmodel,2)));
  mainseed = [temp zeros(1,size(hrfmodel,2))];
  mainlb = [repmat(-Inf,[1 length(temp)]) repmat(-Inf,[1 size(hrfmodel,2)])];
  mainub = [repmat(Inf,[1 length(temp)]) repmat(Inf,[1 size(hrfmodel,2)])];
elseif ~iscell(prfmodel) & ~iscell(hrfmodel) & (flag{1}==0 | flag{1}==2)
  temp = choose(isequal(prfmodel,0),numstim,size(prfmodel,2));
  mainseed = [ones(1,temp*size(hrfmodel,2))];
  mainlb = [repmat(-Inf,[1 temp*size(hrfmodel,2)])];
  mainub = [repmat(Inf,[1 temp*size(hrfmodel,2)])];
else
  die;
end

% deal with nonlinearity seed and bounds
switch mode{1}
case 0
  paramsB = [1 Inf];  % extra parameters to tack on
  nlseed = [];
  nllb = [];
  nlub = [];
case 1
  paramsB = [Inf];
  nlseed = mode{2};  %[1];
  nllb = mode{3}(1);  %[0];
  nlub = mode{3}(2);  %[Inf];
case 2
  paramsB = [];
  nlseed = mode{2};  %[1 10];
  nllb = mode{3}(1,:);  %[0 0];
  nlub = mode{3}(2,:);  %[Inf Inf];
end

% handle fitting in the cross-validation cases
if ~iscell(wantresample) && numel(wantresample) > 1

  % do the fitting
  params = []; numiters = []; responseB = []; rawdata = []; signal = [];
  drift = []; driftchunks = []; modelfit = []; responseBsep = {}; modelfitsep = {};
  if isnan(maxpolydeg)
    polyparams = []; polymeans = [];
  else
    polyparams = NaN*zeros(numsets,maxpolydeg+1); polymeans = NaN*zeros(1,numsets);
  end
  for cnt=1:size(wantresample,1)
  
    % figure out train and test indices
    test = find(-wantresample(cnt,:) < 0);  % NO LONGER REQUIRED: assert(~isempty(test));
    train = find(-wantresample(cnt,:) > 0); assert(~isempty(train));
    
    % precompute
    tmatrix0 = sparse(blkdiag(tmatrix{train}));
    data = tmatrix0*cat(1,response{train}); datastd = choose(std(data)==0,1,std(data));

    % fit the training data
    if ~iscell(prfmodel) & ~iscell(hrfmodel) & (flag{1}==0 | flag{1}==2) & mode{1}==0  % in this very special case, we can do OLS or regularized-regression fitting
      [d,X] = evaluateprfmodel(prfmodel,hrfmodel,flag,[mainseed 1 Inf],stimulus(train));
      if flag{1}==0
        params0 = olsmatrix(tmatrix0*X) * data;
        params(cnt,:) = [params0' 1 Inf];
      else
        [h,dc,numiters0] = gradientdescent_wrapper(data,tmatrix0*X,flag{2},[],flag{3});
        params(cnt,:) = [full(h') 1 Inf];
        numiters(cnt) = numiters0;
      end
    else

      % if nonlinear case, optimize everything except the nonlinear stuff to figure out params0
      if ~isempty(nlseed)
        [params0,d,d,exitflag(cnt),output] = ...  % THIS TAKES THE MOST EXECUTION TIME
          lsqnonlin(@(params) (data - nanreplace(tmatrix0*evaluateprfmodel(prfmodel,hrfmodel,flag, ...
                                                                           [paramsA params nlseed paramsB],stimulus(train)),0,2)) / datastd, ...
            [mainseed],[mainlb],[mainub],options);
      % if linear case, just use the default seed
      else
        params0 = mainseed;
      end
        
      % now, optimize everything
      [params0,d,d,exitflag(cnt),output] = ...  % THIS TAKES THE MOST EXECUTION TIME
        lsqnonlin(@(params) (data - nanreplace(tmatrix0*evaluateprfmodel(prfmodel,hrfmodel,flag, ...
                                                                         [paramsA params paramsB],stimulus(train)),0,2)) / datastd, ...
          [params0 nlseed],[mainlb nllb],[mainub nlub],options);
      params(cnt,:) = [paramsA params0 paramsB];
      numiters(cnt) = output.iterations;
      assert(exitflag(cnt) >= -1);
    end

    % deal with polyparams and polymeans for the training data    
    signal0 = evaluateprfmodel(prfmodel,hrfmodel,flag,params(cnt,:),stimulus(train));
    rawdata0 = cat(1,response{train});
    if isnan(maxpolydeg)
      polyfit0 = zeros(size(signal0));  % hack it such that the effective drift is just zeros
    else
      [polyparams0,polyfit0] = olscontrol(blkdiag(polyregressors{train}),rawdata0, ...
        [signal0 blkdiag(extraregressors{train})]);
      polyparams(train,:) = reshape(polyparams0,maxpolydeg+1,[])';
      polymeans(train) = subsetfun(@mean,polyfit0',time(train));
    end
    rawdata = cat(1,rawdata,rawdata0);
    signal = cat(1,signal,signal0);
    drift = cat(1,drift,polyfit0);
    driftchunks = [driftchunks time(train)];

    % construct the data and model fit for the testing data (polynomials removed)
    temp = sparse(passmulti(@blkdiag,polymatrix(test),1));
    responseB = cat(1,responseB,temp*cat(1,response{test}));
    modelfit = cat(1,modelfit,temp*evaluateprfmodel(prfmodel,hrfmodel,flag,params(cnt,:),stimulus(test)));
    for zz=1:length(test)  % separated by cross-validation run
      temp = sparse(passmulti(@blkdiag,polymatrix(test(zz)),1));
      responseBsep{end+1} = temp*response{test(zz)};
      modelfitsep{end+1} = temp*evaluateprfmodel(prfmodel,hrfmodel,flag,params(cnt,:),stimulus(test(zz)));
    end

  end

% handle fitting in the other cases
else

  % precompute if possible (when wantresample < 0, we can't precompute here)
  if (iscell(wantresample) && wantresample{1} >= 0) || (~iscell(wantresample) && wantresample >= 0)
    tmatrix0 = sparse(blkdiag(tmatrix{:}));
    data = tmatrix0*cat(1,response{:}); datastd = choose(std(data)==0,1,std(data));
  end

  % do the fitting
  params = []; numiters = []; polyparams = []; polymeans = []; rawdata = []; signal = []; drift = []; driftchunks = [];
  for cnt=1:length(bootfun)
    
    % handle the bootstrap data points case [THIS IS STARTING TO GET UGLY!]
    if (iscell(wantresample) && wantresample{1} >= 0) || (~iscell(wantresample) && wantresample >= 0)
      if ~iscell(prfmodel) & ~iscell(hrfmodel) & (flag{1}==0 | flag{1}==2) & mode{1}==0  % in this very special case, we can do OLS or regularized-regression fitting
        [d,X] = evaluateprfmodel(prfmodel,hrfmodel,flag,[mainseed 1 Inf],stimulus);
        if flag{1}==0
          params0 = feval(bootfun{cnt},olsmatrix(tmatrix0*X)')' * feval(bootfun{cnt},data);
          params(cnt,:) = [params0' 1 Inf];
        else
          [h,dc,numiters0] = gradientdescent_wrapper(feval(bootfun{cnt},data),feval(bootfun{cnt},tmatrix0*X),flag{2},[],flag{3});
          params(cnt,:) = [full(h') 1 Inf];
          numiters(cnt) = numiters0;
        end
      else
      
        % if nonlinear case, optimize everything except the nonlinear stuff to figure out params0
        if ~isempty(nlseed)
          [params0,d,d,exitflag(cnt),output] = ...  % THIS TAKES THE MOST EXECUTION TIME
            lsqnonlin(@(params) feval(bootfun{cnt},(data - nanreplace(tmatrix0*evaluateprfmodel(prfmodel,hrfmodel,flag, ...
                                                                                                [paramsA params nlseed paramsB],stimulus),0,2)) / datastd), ...
              [mainseed],[mainlb],[mainub],options);
        else
          params0 = mainseed;
        end
        
        % now, optimize everything
        [params0,d,d,exitflag(cnt),output] = ...  % THIS TAKES THE MOST EXECUTION TIME
          lsqnonlin(@(params) feval(bootfun{cnt},(data - nanreplace(tmatrix0*evaluateprfmodel(prfmodel,hrfmodel,flag, ...
                                                                                              [paramsA params paramsB],stimulus),0,2)) / datastd), ...
            [params0 nlseed],[mainlb nllb],[mainub nlub],options);
        params(cnt,:) = [paramsA params0 paramsB];
        numiters(cnt) = output.iterations;
        assert(exitflag(cnt) >= -1);
      end

    % handle the bootstrap runs case
    else

      % precompute 
      temp = feval(bootfun{cnt},tmatrix);
      tmatrix0 = sparse(blkdiag(temp{:}));
      data = tmatrix0*catcell(1,feval(bootfun{cnt},response));
      datastd = choose(std(data)==0,1,std(data));

      % do it
      if ~iscell(prfmodel) & ~iscell(hrfmodel) & (flag{1}==0 | flag{1}==2) & mode{1}==0  % in this very special case, we can do OLS or regularized-regression fitting
        [d,X] = evaluateprfmodel(prfmodel,hrfmodel,flag,[mainseed 1 Inf],feval(bootfun{cnt},stimulus));
        if flag{1}==0
          params0 = olsmatrix(tmatrix0*X) * data;
          params(cnt,:) = [params0' 1 Inf];
        else
          [h,dc,numiters0] = gradientdescent_wrapper(data,tmatrix0*X,flag{2},[],flag{3});
          params(cnt,:) = [full(h') 1 Inf];
          numiters(cnt) = numiters0;
        end
      else
        
        % if nonlinear case, optimize everything except the nonlinear stuff to figure out params0
        if ~isempty(nlseed)
          [params0,d,d,exitflag(cnt),output] = ...  % THIS TAKES THE MOST EXECUTION TIME
            lsqnonlin(@(params) (data - nanreplace(tmatrix0*evaluateprfmodel(prfmodel,hrfmodel,flag, ...
                                                                             [paramsA params nlseed paramsB],feval(bootfun{cnt},stimulus)),0,2)) / datastd, ...
              [mainseed],[mainlb],[mainub],options);
        % if linear case, just use the default seed
        else
          params0 = mainseed;
        end
        
        % now, optimize everything
        [params0,d,d,exitflag(cnt),output] = ...  % THIS TAKES THE MOST EXECUTION TIME
          lsqnonlin(@(params) (data - nanreplace(tmatrix0*evaluateprfmodel(prfmodel,hrfmodel,flag, ...
                                                                           [paramsA params paramsB],feval(bootfun{cnt},stimulus)),0,2)) / datastd, ...
            [params0 nlseed],[mainlb nllb],[mainub nlub],options);
        params(cnt,:) = [paramsA params0 paramsB];
        numiters(cnt) = output.iterations;
        assert(exitflag(cnt) >= -1);
      end
    end

    % deal with polyparams and polymeans for the training data    
    signal0 = evaluateprfmodel(prfmodel,hrfmodel,flag,params(cnt,:),stimulus);
    rawdata0 = cat(1,response{:});
    if isnan(maxpolydeg)
      polyfit0 = zeros(size(signal0));  % hack it such that the effective drift is just zeros
    else
      [polyparams0,polyfit0] = olscontrol(blkdiag(polyregressors{:}),rawdata0, ...
        [signal0 blkdiag(extraregressors{:})]);
      polyparams(:,:,cnt) = reshape(polyparams0,maxpolydeg+1,[])';
      polymeans(cnt,:) = subsetfun(@mean,polyfit0',time);
    end
    rawdata = cat(1,rawdata,rawdata0);
    signal = cat(1,signal,signal0);
    drift = cat(1,drift,polyfit0);
    driftchunks = [driftchunks time];

  end

  % construct the data and model fit (polynomials removed)
  temp = sparse(blkdiag(exmatrix{:})*blkdiag(polymatrix{:}));
  responseB = temp*cat(1,response{:});
  modelfit = temp*mean(rowfun(params,@(x) evaluateprfmodel(prfmodel,hrfmodel,flag,x,stimulus),2),2);
  responseBsep = {}; modelfitsep = {};
  for zz=1:numsets
    temp = sparse(blkdiag(exmatrix{zz})*blkdiag(polymatrix{zz}));
    responseBsep{end+1} = temp*response{zz};
    modelfitsep{end+1} = temp*mean(rowfun(params,@(x) evaluateprfmodel(prfmodel,hrfmodel,flag,x,stimulus(zz)),2),2);
  end

end

% calculate hrf and hrfsc
hrf = []; hrfsc = [];
for xx=1:size(params,1)
  [d,d,d,hrf0] = evaluateprfmodel(prfmodel,hrfmodel,flag,params(xx,:),stimulus);
  hrfsc(xx) = feval(hrfnormfun,hrf0);
  hrf(:,xx) = hrf0 ./ hrfsc(xx);
end

% calculate beta weights [THIS IS POTENTIALLY THE UGLIEST CODE EVER.  SHOULD CLEAN UP.]
if wantbetaweights

  % handle the regular and bootstrap cases
  if isscalar(wantresample) || iscell(wantresample)

    % the inseparable case is weird, handle it here
    if ~iscell(hrfmodel) & (flag{1}==0 | flag{1}==2)
      betas = hrfsc;
      
    % handle other cases
    else
      numprf = choose(isequal(prfmodel,0),numstim,size(prfmodel,2));    % calc number of PRF things
      betas = bsxfun(@times,params(:,1:numprf)',hrfsc);
    end
  
  % handle the cross-validation cases
  else

    % figure out indices
    indices = {}; hrfindex = [];
    for xx=1:size(wantresample,1)
      vals = -wantresample(xx,:);
  
      % do training runs
      if ismember(derivemode(1),[1 3])
        indices{end+1} = find(vals > 0);
        hrfindex(end+1) = xx;
      end
      if ismember(derivemode(1),[2 3])
        for zz=1:max(vals)
          indices{end+1} = find(vals==zz);
          hrfindex(end+1) = xx;
        end
      end
      
      % do testing runs
      if ismember(derivemode(2),[1 3])
        indices{end+1} = find(vals < 0);
        hrfindex(end+1) = xx;
      end
      if ismember(derivemode(2),[2 3])
        for zz=-1:-1:min(vals)
          indices{end+1} = find(vals==zz);
          hrfindex(end+1) = xx;
        end
      end
  
    end
    
    % then, derive the beta weights
    betas = [];
    for p=1:length(indices)
    
      if ~isempty(indices{p})  % blkdiag gets mad if it gets no arguments.  also, the assignment involving "end+1" will fail with the empty case.
                               % so, let's explicitly detect to avoid this degenerate case.
    
        % calc the matrix that removes the polynomials and extra regressors
        tmatrix0 = sparse(blkdiag(tmatrix{indices{p}}));
    
        % calc the design matrix (note that we simulate the inseparable mode)
        [d,X] = evaluateprfmodel(prfmodel,hrf(:,hrfindex(p)),0,[ones(1,numstim) 1 Inf],stimulus(indices{p}));
        
        % fit the model
        betas(:,end+1) = olsmatrix(tmatrix0*X) * tmatrix0*cat(1,response{indices{p}});
  
      end
      
    end
    
  end

else
  betas = [];
end

% how well did we do?
r = feval(metric,modelfit,responseB);
rrun = cellfun(metric,modelfitsep,responseBsep);

% deal with empty cases
r = choose(isempty(r),NaN,r);
rrun = choose(isempty(rrun),NaN,rrun);

% deal with bootstrapping
if ~isequal(wantresample,0)
  paramsse = std(params,[],1);
else
  paramsse = NaN*zeros(1,size(params,2));
end

% should we re-do?
if any(exitflag==0) && length(maxiter)==3
  if r >= maxiter(2)
      % notice that we use wantresampleORIG here
    [params,paramsse,r,rrun,polyparams,polymeans,numiters,hrf,betas,signal,drift] = ...
      fitprf(stimulus,response,prfmodel,hrfmodel,flag,maxpolydeg,ar,mode,maxiter(3),wantresampleORIG,tol,extraopt,extraregressors,hrfnormfun,derivemode);
  end
end

% encapsulate if necessary
if nargout==1
  params = cell2struct({params,paramsse,r,rrun,polyparams,polymeans,numiters,hrf,betas,signal,drift}, ...
                       {'params' 'paramsse' 'r' 'rrun' 'polyparams' 'polymeans' 'numiters' 'hrf' 'betas' 'signal' 'drift'},2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTION (FOR INTERNAL DEVELOPMENT PURPOSES ONLY)

function [response,ts,intermediate,hrfactual] = evaluateprfmodel(prfmodel,hrfmodel,flag,params,stimulus)

% function [response,ts,intermediate] = evaluateprfmodel(prfmodel,hrfmodel,flag,params,stimulus)
%
% <prfmodel>,<hrfmodel>,<flag>,<params>,<stimulus> are as in fitprf.m.  all optional parameters should be supplied.
%
% given the PRF and HRF models, parameter values, and the stimulus,
% return <response> as the predicted response timecourse with dimensions
% time x 1.  also return <ts>, which is meaningful only in case 5.
% in case 5, <ts> is time x n*b with the model right before the linear
% weighting step.
%
% also return <intermediate> which is a cell vector with
% various contents that represent the intermediate calculations that
% we make.  this cell vector is incrementally built up
% as we process each <stimulus> case.  for a given <stimulus> case,
% if case 1, 2, or 5, we add to <intermediate> these quantities:
% stimulus*prf, divnorm(stimulus*prf), and conv2 with hrf.  if case 3,
% we add to <intermediate> these quantities: stimulus*prf, 
% divnorm(stimulus*prf), and weighting with prf weights.  if case 4,
% we add to <intermediate> these quantities: stimulus*prf, 
% divnorm(stimulus*prf), weighting with prf weights, and conv2 with hrf.
%
% (note that there is redundancy between <ts> and 
% <intermediate>, but this is intentional, since we do not calculate
% <intermediate> unless the user requests it.)
%
% for cases 1 and 2, the predicted response is schematically
% given by conv2(divnorm(stimulus*prf),hrf) * weights, where
% stimulus is the contrast mask with dimensions time x stim;
% prf is the spatial PRF model with dimensions stim x 1 and is
% obtained by evaluating a parametric function; divnorm is a point nonlinearity
% (see divnorm.m); hrf is one or more timecourses with dimensions t x b;
% and weights are linear weights with dimensions b x 1.  note that in
% the parametric HRF case, b==1 and the HRF timecourse is obtained by 
% evaluating a function, whereas in the basis HRF case, HRF basis functions
% are given directly by the user.
%
% for case 3, the predicted response is schematically
% given by conv2(divnorm(stimulus*prf) * weights,hrf), where
% stimulus is the contrast mask with dimensions time x stim;
% prf is one or more basis functions with dimensions stim x n;
% divnorm is a point nonlinearity (see divnorm.m); weights are
% linear weights with dimensions n x 1; and hrf is an HRF timecourse
% with dimensions t x 1 and is obtained by evaluating a function.
%
% for case 4, the predicted response is schematically
% given by conv2(divnorm(stimulus*prf) * a_weights,hrf * b_weights), where
% stimulus is the contrast mask with dimensions time x stim;
% prf is one or more basis functions with dimensions stim x n;
% a_weights is [1 a2 a3 ...]', which are linear weights with dimensions
% n x 1; hrf is one or more basis functions with dimensions t x b;
% and b_weights are linear weights with dimensions b x 1.  notice that
% this case imposes separability between the PRF and the HRF.
%
% for case 5, the predicted response is schematically
% given by conv2(divnorm(stimulus*prf),hrf) * weights, where
% stimulus is the contrast mask with dimensions time x stim;
% prf is one or more basis functions with dimensions stim x n;
% hrf is one or more basis functions with dimensions t x b;
% and weights are linear weights with dimensions n*b x 1.
% notice that this case does not impose separability between 
% the PRF and the HRF.  also, note that the conv2 is not literally
% right; really what we mean is the concatenation of the convolutions
% of divnorm(stimulus*prf) with each column of hrf.
%
% note that when <stimulus> is a cell vector, we automatically concatenate
% the different time-series together.

% HANDLE INTERNAL SPECIAL CASE (DUE TO CALL IN FITPRF.M):  [HACKY]
if iscell(stimulus) && length(stimulus)==0
  response = [];
  ts = [];
  intermediate = [];
  hrfactual = [];
  return;
end

% input
if ~iscell(stimulus)
  stimulus = {stimulus};
end
if ~iscell(flag)
  flag = {flag};
end

% calc
if iscell(prfmodel)
  numparamsprf = length(prfmodel{1});  % number of parameters in PRF
else
  if isequal(prfmodel,0)
    numparamsprf = size(stimulus{1},2);
  else
    numparamsprf = size(prfmodel,2);
  end
end
if iscell(hrfmodel)
  numparamshrf = length(hrfmodel{1});  % number of parameters in HRF
else
  numparamshrf = size(hrfmodel,2);
end

% calculate PRF (stim x 1)
if iscell(prfmodel)
  prf = feval(prfmodel{3},params(1:numparamsprf));
else
  if isequal(prfmodel,0)
    prf = 1;
  else
    prf = prfmodel;
  end
end

% calculate HRF (t x 1 or t x b).
%   we also calculate hrfactual (t x 1) which is only for internal use.
%   note that in the inseparable case, we return the first PRF component's HRF timecourse.
if iscell(hrfmodel)
  hrf = feval(hrfmodel{3},params(numparamsprf+(1:numparamshrf)));
  hrfactual = hrf;
else
  hrf = hrfmodel;
  if flag{1}==0 || flag{1}==2
    temp = params(1:numparamsprf*numparamshrf);
    hrfactual = hrf * vflatten(temp(1:numparamsprf:end));
  else
    hrfactual = hrf * vflatten(params(end-2-numparamshrf+1:end-2));
  end
end

% do it
response = []; ts = []; intermediate = {};
for p=1:length(stimulus)

  % calculate "neural" activity (time x 1 or time x n)
  if isequal(prf,1)
    neural = stimulus{p};
  else
    neural = stimulus{p}*prf;
  end
  if nargout >= 3, intermediate{end+1} = neural;, end
  
  % apply divisive-normalization nonlinearity
  neural = divnorm(neural,params(end-1),params(end));
  if nargout >= 3, intermediate{end+1} = neural;, end
  
  % finish up
  if iscell(prfmodel)

    % do HRF convolution (time x b)
    ts0 = conv2(neural,hrf);
    ts0 = ts0(1:size(neural,1),:);
    if nargout >= 3, intermediate{end+1} = ts0;, end
    
    % weighted sum
    if iscell(hrfmodel) && hrfmodel{4}==1
      response = cat(1,response,ts0);  % no weighting necessary
    else
      numweights = choose(iscell(hrfmodel),1,numparamshrf);
      response = cat(1,response,ts0 * params(end-2-numweights+1:end-2)');
    end
    if nargout >= 2, ts = cat(1,ts,ts0);, end

  else
  
    % this is the tricky case --- both PRF and HRF are basis functions and we want inseparability
    if ~iscell(hrfmodel) && (flag{1}==0 || flag{1}==2)
    
      % do convolution for each HRF basis function (time x n*b)
      ts0 = conv2(neural,upsamplematrix(hrf,size(neural,2),2,1));
        %%ts0 = blockproc(hrf,[size(hrf,1) 1],@(x) conv2(neural,x));
      ts0 = ts0(1:size(neural,1),:);
      if nargout >= 3, intermediate{end+1} = ts0;, end

      % weighted sum
      response = cat(1,response,ts0*params(1:numparamsprf*numparamshrf)');
      if nargout >= 2, ts = cat(1,ts,ts0);, end
    
    % easier cases
    else
  
      % weight the PRF basis functions and then do HRF convolution (time x b)
      ts0 = neural * params(1:numparamsprf)';
      if nargout >= 3, intermediate{end+1} = ts0;, end
      numtimepts = size(ts0,1);
      ts0 = conv2(ts0,hrf);
      ts0 = ts0(1:numtimepts,:);  % don't add to intermediate yet; check case *** below
      
      % weighted sum if necessary
      if iscell(hrfmodel)
        response = cat(1,response,ts0);  % we're done already
      else
        if nargout >= 3, intermediate{end+1} = ts0;, end  % *** ok, add to intermediate
        response = cat(1,response,ts0*params(numparamsprf+(1:numparamshrf))');  % do weighted sum
      end
      if nargout >= 2, ts = cat(1,ts,ts0);, end

    end
  
  end  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK:

% TODO: should we return the model fit?  and show scatter plot of model fit against data?
%
% TODO: does it make any sense to include regularization?  e.g. "automatic" smoothing of the FIR HRF stuff.
%       it would be painful to implement it for nonlinear optimization case (since we would have to discretize the regularization parameter and fit many time).
%       but for linear case, it's potentially feasible easily...  we regularized inseparable model, didn't seem to help.  maybe instead of L2 shrinkage, we
%       could impose smoothness prior on the FIR.  so FIR model plus penalty for roughness.  [but there are other ways to impose smoothness!]
%
% TODO: perhaps the bottleneck is the tmatrix multiplciation?  should we fit the nuisance regressors as part of the lsqnonlin??? 
%       [actually, fitting the nuisance regressors seemed to make optimization convergence slower]
%
% TODO: would back-and-forth iterative implementation (fit hrf, fit prf, etc.) be faster? [but doesn't apply to all cases (e.g. spline)]
%
% TODO: worry about local minima in prf fitting? [in concert with intelligent initial seed] convergence problem?

% ??polyparams and polymeans is somehow outside fo the bootstrap.  i need to figure this out.

% REMOVED:
% the parameter for an AR(1) noise model.  default is 0, 
% %   which in effect means to not use the noise model.

% REMOVED:
% ,meanint,driftstd,signalrms,noiserms
% %  <meanint> is the mean of the mean of each run's fitted nuisance functions.
% %  <driftstd> is the mean of the std dev of each run's fitted nuisance functions.
% %  <signalrms> is the RMS of the signal (see <signal>)
% %  <noiserms> is the RMS of the noise (defined as the raw data minus the signal and the drift)
% %    note that for the above four output variables, the following is true:
% %    - in the cross-validation case, only the training runs are considered.
% %    - when there are multiple resamplings (e.g. multiple cross-validations or multiple bootstraps),
% %      this is treated in the natural way (e.g. we just act as if the extra resamplings are just 
% %      additional cases to take the mean over).
% % calculate some signal analysis stuff
% noise = rawdata - signal - drift;
% meanint = mean(chunkfun(drift,driftchunks,@(x)nanmean(x,1)));     % mean of the mean of each run's drift
% driftstd = mean(chunkfun(drift,driftchunks,@(x)nanstd(x,[],1)));  % mean of the std of each run's drift
% signalrms = sqrt(mean(signal.^2));
% noiserms = sqrt(mean(noise.^2));

% 
% for example, suppose <wantresample> is -[1 1 -1 -2; 0 0 0 1] and suppose
% %    <derivemode> is [1 3].  then, the columns of <betas> will be obtained in the following way.  first, estimate
% %    the HRF by training on runs 1 and 2.  then derive beta weights for runs 1 and 2 as a group.  then derive 
% %    beta weights for runs 3 and 4 as a group.  then derive beta weights for run 3.  then derive beta weights 
% %    for run 4.  second, estimate HRF by training on run 4.  then derive beta weights for run 4.  so in this
% %    example, <betas> will have a total of five columns.
