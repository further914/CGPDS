function g = vargplvmPointGradient(x, model, y)

% VARGPLVMPOINTGRADIENT Wrapper function for gradient of a single point.
%
%	Description:
%
%	G = VARGPLVMPOINTGRADIENT(X, MODEL, Y) is a wrapper function for the
%	gradient of the log likelihood with respect to a point in the latent
%	space. The GP-LVM model is one that is assumed to have already been
%	trained.
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


%	Copyright (c) 2009, 2011 % COPYRIGHT Michalis K. Titsias and Neil D. Lawrence
%	Copyright (c) 2011 % COPYRIGHT Andreas C. Damianou
% 	vargplvmPointGradient.m SVN version 1430
% 	last update 2011-06-12T16:31:32.000000Z

if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise  %%% RE-OPT-CODE-NEW
        [vardistx, model] = vargplvmPartExpand(model, x); %%% RE-OPT-CODE-NEW
    elseif isfield(model.dynamics,'onlyTest') && model.dynamics.onlyTest
        vardistx = model.vardistx;
        vardistx = vardistExpandParam(vardistx, x);
    else %%% RE-OPT-CODE-NEW
        % this is doing the expand
        x = reshape(x, model.N+size(y,1), model.dynamics.q*2);
        xtrain = x(1:model.N,:);
        xtest = x(model.N+1:end,:);
        model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
        vardistx = vardistExpandParam(model.vardistx, xtest);
        % end of expand
    end %%% RE-OPT-CODE-NEW
else
    vardistx = model.vardistx;
    vardistx = vardistExpandParam(vardistx, x);
end
fhandle = str2func([model.type 'PointLogLikeGradient']);
g = -fhandle(model, vardistx, y);
% g = - vargplvmPointLogLikeGradient(model, vardistx, y);
