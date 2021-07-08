function f = vargplvmPointObjective(x, model, y)

% VARGPLVMPOINTOBJECTIVE Wrapper function for objective of a single point in latent space and the output location..
%
%	Description:
%
%	F = VARGPLVMPOINTOBJECTIVE(X, MODEL, Y) provides a wrapper function
%	for the negative log probability of a given data point under the
%	posterior distribution of the Gaussian process induced by the
%	training data.
%	 Returns:
%	  F - the negative of the log probability of the given data point
%	   under the posterior distribution induced by the training data.
%	 Arguments:
%	  X - location in input space for the point.
%	  MODEL - the model structure for which the negative log probability
%	   of the given data under the posterior is to be computed.
%	  Y - the location in data space for the point.
%	
%
%	See also
%	VARGPLVMCREATE, VARGPLVMPOINTLOGLIKELIHOOD, VARGPLVMOPTIMISEPOINT


%	Copyright (c) 2009 Michalis K. Titsias and Neil D. Lawrence
% 	vargplvmPointObjective.m SVN version 1408
% 	last update 2011-06-02T14:43:35.000000Z

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
fhandle = str2func([model.type 'PointLogLikelihood']);
f = -fhandle(model, vardistx, y);
% f = - vargplvmPointLogLikelihood(model, vardistx, y);