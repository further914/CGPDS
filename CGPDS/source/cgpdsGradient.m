function g = cgpdsGradient(params, model ,stepSize, samd, samr)

% cgpdsGradient CGPDS gradient wrapper.
%
%	Description:
%
%	G = cgpdsGradient(PARAMS, MODEL) is a wrapper function for the
%	gradient of the negative log likelihood of the CGPDS
%	model with respect to the latent postions and parameters.
%	 Returns:
%	  G - the gradient of the negative log likelihood with respect to
%	   the latent positions and the parameters at the given point.
%	 Arguments:
%	  PARAMS - vector of parameters and latent postions where the
%	   gradient is to be evaluated.
%	  MODEL - the model structure into which the latent positions and
%	   the parameters will be placed.
%	


%	Copyright (c) 2009 - 2011 Michalis K. Titsias
%	Copyright (c) 2006, 2005, 2010-2011 Neil D. Lawrence
% 	vargplvmGradient.m SVN version 1582
% 	last update 2011-09-07T17:23:12.301181Z

model = modelExpandParam(model, params, samd, samr);
g = - modelLogLikeGradients(model ,stepSize, samd, samr);

% sum gradients of tied parameters, then assign corresponding summed gradients to each
% group of tied parameters
% if isfield( model, 'ties' )
%     g = g * model.T; % model.T == model.ties' * model.ties;
% end
% fprintf(1,'# G: %.13f\n',sum(abs(g))); %%% DEBUG
