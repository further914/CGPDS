function varsigma = vargplvmPosteriorVar(model, X)

% VARGPLVMPOSTERIORVAR variances of the posterior at points given by X.
%
%	Description:
%
%	[MU, SIGMA] = VARGPLVMPOSTERIORVAR(MODEL, X) returns the posterior
%	mean and variance for a given set of points.
%	 Returns:
%	  MU - the mean of the posterior distribution.
%	  SIGMA - the variances of the posterior distributions.
%	 Arguments:
%	  MODEL - the model for which the posterior will be computed.
%	  X - the input positions for which the posterior will be computed.
%	
%
%	See also
%	GPPOSTERIORMEANVAR, VARGPLVMCREATE


%	Copyright (c) 2005, 2006 Neil D. Lawrence
% 	vargplvmPosteriorVar.m SVN version 1408
% 	last update 2011-06-02T14:43:36.000000Z
 

%% ORIGINAL
model.K_uf = kernCompute(model.kern, model.X_u, X);
model.A = (1/model.beta)*model.K_uu + model.K_uf*model.K_uf';
[model.Ainv, U] = pdinv(model.A);
varsigma = gpPosteriorVar(model, X);

%%% TEMP
%[void, varsigma] = vargplvmPosteriorMeanVar(model, X);