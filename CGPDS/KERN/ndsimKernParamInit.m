function kern = ndsimKernParamInit(kern)

% NDSIMKERNPARAMINIT SIM kernel parameter initialisation.
%
%	Description:
%	The single input motif (SIM) kernel is specifically designed for
%	working with gene networks where there is assumed to be a single
%	transcription factor controlling several genes. If each gene is
%	related to the transcription factor through the following
%	differential equation,
%	
%	dx(t)/dt = B + S f(t-delta) - D x(t),
%	
%	where D is a decay term, S is a response term, delta is a time delay
%	and B is an initial level. Then if f(t) is assumed to come from a
%	Gaussian process with an RBF covariance function x(t) is a Gaussian
%	process with a covariance function provided by the single input
%	motif kernel.
%	
%	The kernel is designed to interoperate with the multiple output
%	block kernel so that f(t) can be inferred given several different
%	instantiations of x(t) (associated with different genes).
%	
%	By default the parameters (B, S, delta and D) are constrained positive. If
%	kern.options.isNegativeS is set true then the parameter S is allowed to go
%	negative.
%	
%
%	KERN = NDSIMKERNPARAMINIT(KERN) initialises the single input motif
%	kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%	
%
%	See also
%	KERNCREATE, KERNPARAMINIT, SIMKERNCOMPUTE


%	Copyright (c) 2006, 2009 Neil D. Lawrence
%	Copyright (c) 2011 Jaakko Peltonen
% 	ndsimKernParamInit.m SVN version 1532
% 	last update 2011-08-03T14:11:43.325503Z

if kern.inputDimension > 1
  error('NDSIM kernel is only valid for one-D input.')
end

if isfield(kern, 'options') && isfield(kern.options, 'gaussianInitial') && ...
      kern.options.gaussianInitial,
  kern.gaussianInitial = 1;
  kern.initialVariance = 1;
else
  kern.gaussianInitial = 0;
end

kern.delay = 0;
kern.initVal = 1;
kern.variance = 1;
kern.inverseWidth = 1;
sensitivityparam_index=2;

if kern.gaussianInitial,
  kern.nParams = 3;
else
  kern.nParams = 2;
end

if isfield(kern, 'options') ...
      && isfield(kern.options, 'isNegativeS') ...
      && kern.options.isNegativeS,
  kern.isNegativeS = true;  
  kern.sensitivity = 1;
  positiveparams = setdiff(1:kern.nParams, sensitivityparam_index);
else
  kern.isNegativeS = false;
  positiveparams = 1:kern.nParams;
end

for k=1:length(positiveparams),
  kern.transforms(k).type = 'sigmoidab';
  kern.transforms(k).index = positiveparams(k);
  kern.transforms(k).transformsettings = [0 1e6];
end;

kern.isStationary = false;
kern.isNormalised = false;
kern.positiveTime = true;

