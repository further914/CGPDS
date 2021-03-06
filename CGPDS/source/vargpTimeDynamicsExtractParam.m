function [params, names] = vargpTimeDynamicsExtractParam(model, samd, samr)

% VARGPTIMEDYNAMICSEXTRACTPARAM Extract parameters from the GP time dynamics model.
%
%	Description:
%	
%	Description:
%	
%	PARAMS = VARGPTIMEDYNAMICSEXTRACTPARAM(MODEL) extracts the model
%	parameters from a structure containing the information about a
%	Gaussian process dynamics model.
%	Returns:
%	PARAMS - a vector of parameters from the model.
%	Arguments:
%	MODEL - the model structure containing the information about the
%	model.
%	DESC does the same as above, but also returns parameter names.
%	ARG model : the model structure containing the information about
%	the model.
%	RETURN params : a vector of parameters from the model.
%	RETURN names : cell array of parameter names.
%	
%	
%	
%	See also
%	GPEXTRACTPARAM, GPTIMEDYNAMICSCREATE, GPTIMEDYNAMICSEXPANDPARAM, MODELEXTRACTPARAM
% 	vargpTimeDynamicsExtractParam.m SVN version 1412
% 	last update 2011-07-04T19:08:51.964568Z

%	Based on GPTIMEDYNAMICSEXTRACTPARAM

if nargout > 1
  returnNames = true; % Also return parameter names in an array
else
  returnNames = false;
end  

% Parameters are extracted in the following order: (notation: % parameter{size})
% [dynamicsVardistParams{dynamics.vardist.nParams}        % mu_bar, lambda
%  dynamics.kernParams{dynamics.kern.nParams}]            % sigmaRBF, lt,
%  sigmaWhite

% Variational parameters (reparametrized means and covariances)
if returnNames
  [varParams, varNames] = modelExtractParam(model.vardist, samd, samr);
  names = varNames;
else
  %varParams = vardistExtractParam(model.vardist);
  varParams = modelExtractParam(model.vardist, samd, samr);
end
params = varParams;


% Kernel parameters  
if returnNames
  [kernParams, kernParamNames] = kernExtractParam(model.kern); 
  for i = 1:length(kernParamNames)
    kernParamNames{i} = ['Kernel, ' kernParamNames{i}];
  end
  names = {names{:}, kernParamNames{:}};
else
  kernParams = kernExtractParam(model.kern);
end

params = [params kernParams];


