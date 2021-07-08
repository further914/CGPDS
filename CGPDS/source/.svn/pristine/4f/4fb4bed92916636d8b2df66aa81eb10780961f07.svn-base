function [f, g] = vargplvmPointObjectiveGradient(x, model, y)

% VARGPLVMPOINTOBJECTIVEGRADIENT Wrapper function for objective and gradient of a single point in latent space and the output location..
%
%	Description:
%
%	[F, G] = VARGPLVMPOINTOBJECTIVEGRADIENT(X, MODEL, Y) provides a
%	wrapper function for the negative log probability of a given data
%	point under the posterior distribution of the Gaussian process
%	induced by the training data. Also returns the gradient of the
%	negative log probability with respect to the given latent point.
%	 Returns:
%	  F - the negative of the log probability of the given data point
%	   under the posterior distribution induced by the training data.
%	  G - the gradient of the log probability with respect to the given
%	   latent point.
%	 Arguments:
%	  X - location in input space for the point.
%	  MODEL - the model structure for which the negative log probability
%	   of the given data under the posterior is to be computed.
%	  Y - the location in data space for the point.
%	vargplvmOptimisePoint, vargplvmObjective, vargplvmGradient
%	
%
%	See also
%	VARGPLVMCREATE, VARGPLVMPOINTLOGLIKELIHOOD, 


%	Copyright (c) 2009 Michalis K. Titsias and Neil D. Lawrence
% 	vargplvmPointObjectiveGradient.m SVN version 1252
% 	last update 2011-07-04T19:08:52.274565Z

% Check how the optimiser has given the parameters
if size(xvec, 1) > size(xvec, 2)
  % As a column vector ... transpose everything.
  transpose = true;
  x = x';
else
  transpose = false;
end
f = - vargplvmPointLogLikelihood(model, x, y);

if nargout > 1
  g = - vargplvmPointLogLikeGradient(model, x, y);
end
if transpose
  g = g';
end

