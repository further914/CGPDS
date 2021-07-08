function [f, g] = vargplvmObjectiveGradient(params, model)

% VARGPLVMOBJECTIVEGRADIENT Wrapper function for VARGPLVM objective and gradient.
%
%	Description:
%
%	[F, G] = VARGPLVMOBJECTIVEGRADIENT(PARAMS, MODEL) returns the
%	negative log likelihood of a Gaussian process model given the model
%	structure and a vector of parameters. This allows the use of NETLAB
%	minimisation functions to find the model parameters.
%	 Returns:
%	  F - the negative log likelihood of the VARGPLVM model.
%	  G - the gradient of the negative log likelihood of the VARGPLVM
%	   model with respect to the parameters.
%	 Arguments:
%	  PARAMS - the parameters of the model for which the objective will
%	   be evaluated.
%	  MODEL - the model structure for which the objective will be
%	   evaluated.
%	
%	
%
%	See also
%	MINIMIZE, VARGPLVMCREATE, VARGPLVMGRADIENT, VARGPLVMLOGLIKELIHOOD, VARGPLVMOPTIMISE


%	Copyright (c) 2009 Michalis K. Titsias
%	Copyright (c) 2005, 2006 Neil D. Lawrence
% 	vargplvmObjectiveGradient.m SVN version 1473
% 	last update 2011-07-04T19:19:55.244564Z
  
% Check how the optimiser has given the parameters
if size(params, 1) > size(params, 2)
  % As a column vector ... transpose everything.
  transpose = true;
  model = vargplvmExpandParam(model, params');
else
  transpose = false;
  model = vargplvmExpandParam(model, params);
end

f = - vargplvmLogLikelihood(model);
% fprintf(1,'# F: %.13f\n',f); %%% DEBUG
if nargout > 1
  g = - vargplvmLogLikeGradients(model);
%  fprintf(1,'# G: %.13f .\n',sum(abs(g))); %%% DEBUG
end
if transpose
  g = g';
end

