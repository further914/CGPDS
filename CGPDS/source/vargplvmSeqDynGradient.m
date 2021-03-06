function g = vargplvmSeqDynGradient(x, model, y,samd,samr)

% VARGPLVMSEQDYNGRADIENT Wrapper function for gradient of a single point.
%
%	Description:
%
%	G = VARGPLVMSEQDYNGRADIENT(X, MODEL, Y) is a wrapper function for
%	the gradient of the log likelihood with respect to a point in the
%	latent space. The GP-LVM model is one that is assumed to have
%	already been trained.
%	 Returns:
%	  G - the gradient of the log likelihood with respect to the latent
%	   position.
%	 Arguments:
%	  X - the position in the latent space that is being optimised.
%	  MODEL - the trained GP-LVM model that is being optimised.
%	  Y - the position in data space for which the latent point is being
%	   optimised.
%	
%
%	See also
%	VARGPLVMPOINTLOGLIKEGRADIENT, VARGPLVMOPTIMISEPOINT


%	Copyright (c) 2009 % COPYRIGHT Michalis K. Titsias and Neil D. Lawrence
% 	vargplvmSeqDynGradient.m SVN version 1441
% 	last update 2011-06-02T20:53:40.000000Z

% if isfield(model, 'dynamics') && ~isempty(model.dynamics)
%    % this is doing the expand 
%    x = reshape(x, model.N+size(y,1), model.dynamics.q*2);
%    xtrain = x(1:model.N,:);
%    xtest = x(model.N+1:end,:);
%    model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
%    vardistx = vardistExpandParam(model.vardistx, xtest);
%    % end of expand 
% else 
   vardistx = model.vardistx;
   vardistx = vardistExpandParam(vardistx, x);
%end
fhandle = str2func([model.type 'SeqDynLogLikeGradient']);
g = -fhandle(model, vardistx, y,samd,samr);

