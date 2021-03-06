function f = vargplvmSeqDynObjective(x, model, y,samd,samr)

% VARGPLVMSEQDYNOBJECTIVE Wrapper function for objective of a group of points in latent space and the output locations..
%
%	Description:
%
%	F = VARGPLVMSEQDYNOBJECTIVE(X, MODEL, Y) provides a wrapper function
%	for the negative log probability of a group of data points under the
%	posterior distribution of the Gaussian process induced by the
%	training data.
%	 Returns:
%	  F - the negative of the log probability of the given data point
%	   under the posterior distribution induced by the training data.
%	 Arguments:
%	  X - locations in input space for the point.
%	  MODEL - the model structure for which the negative log probability
%	   of the given data under the posterior is to be computed.
%	  Y - the location in data space for the points.
%	
%
%	See also
%	VARGPLVMCREATE, VARGPLVMPOINTLOGLIKELIHOOD, VARGPLVMOPTIMISEPOINT


%	Copyright (c) 2011 Michalis K. Titsias and Neil D. Lawrence
% 	vargplvmSeqDynObjective.m SVN version 1312
% 	last update 2011-04-18T11:54:44.000000Z

% % this is doing the expand 
% x = reshape(x, model.N+size(y,1), model.dynamics.q*2);
% xtrain = x(1:model.N,:);
% xtest = x(model.N+1:end,:);
% model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
% vardistx = vardistExpandParam(model.vardistx, xtest);
% % end of expand 
vardistx = model.vardistx;
vardistx = vardistExpandParam(vardistx, x);
fhandle = str2func([model.type 'SeqDynLogLikelihood']);
f = -fhandle(model, vardistx, y,samd,samr);

%f = -TEMPvargplvmPointLogLikelihoodSTATIC(model, vardistx, y);