function [params, names] = vargplvmExtractParam(model, samd, samr)

% VARGPLVMEXTRACTPARAM Extract a parameter vector from a variational GP-LVM model.
%
%	Description:
%
%	PARAMS = VARGPLVMEXTRACTPARAM(MODEL) extracts a parameter vector
%	from a given VARGPLVM structure.
%	 Returns:
%	  PARAMS - the parameter vector extracted from the model.
%	 Arguments:
%	  MODEL - the model from which parameters are to be extracted.
%	DESC does the same as above, but also returns parameter names.
%	ARG model : the model structure containing the information about
%	the model.
%	RETURN params : a vector of parameters from the model.
%	RETURN names : cell array of parameter names.
%	
%	
%	Modifications: Andreas C. Damianou, 2010-2011
%	
%
%	See also
%	VARGPLVMCREATE, VARGPLVMEXPANDPARAM, MODELEXTRACTPARAM


%	Copyright (c) 2009-2011 Michalis K. Titsias and Neil D. Lawrence
% 	vargplvmExtractParam.m SVN version 1572
% 	last update 2011-08-30T14:57:48.000000Z

%%% Parameters must be returned as a vector in the following order (left to right) 
% - parameter{size} -
% vardistParams{model.vardist.nParams} % mu, S
%       OR
% [dynamicsVardistParams{dynamics.vardist.nParams} dynamics.kernParams{dynamics.kern.nParams}] % mu_bar, lambda
% inducingInputs{model.q*model.k}
% kernelParams{model.kern.nParams}
% beta{prod(size(model.beta))}


if nargout > 1
  returnNames = true;
else
  returnNames = false;
end 

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    % [VariationalParameters(reparam)   dynKernelParameters]
    if returnNames
        [dynParams, dynParamNames] = modelExtractParam(model.dynamics,samd,samr);%vardistExtractParam.m
        names = dynParamNames;
    else
        %%%
        dynParams = modelExtractParam(model.dynamics,samd,samr);%vardistExtractParam.m
    end
    params = dynParams;% add variational parameters :mu and lambda
else
    % Variational parameters 
    if returnNames
        %[varParams, varNames] = vardistExtractParam(model.vardist);
        [varParams, varNames] = modelExtractParam(model.vardist,samd,samr);
        %names = varNames{:}; %%% ORIGINAL
        names = varNames; %%% NEW
    else
        %varParams = vardistExtractParam(model.vardist);
        varParams = modelExtractParam(model.vardist,samd,samr);
    end
    params = varParams;
end


%if does not fix Inducing inputs 
if ~model.fixInducing 
    params =  [params model.X_v(:)'];
    params =  [params model.X_u(:)'];
    
%     for i = 1:size(model.X_v, 1)
%           for j = 1:size(model.X_v, 2)
%               X_vNames{i, j} = ['X_v(' num2str(i) ', ' num2str(j) ')'];
%           end
%     end
%       
%     if returnNames 
%       for i = 1:size(model.X_u, 1)
%           for j = 1:size(model.X_u, 2)
%               X_uNames{i, j} = ['X_u(' num2str(i) ', ' num2str(j) ')'];
%           end
%       end
%       
%       names = {names{:}, X_uNames{:},X_vNames{i, j}};
%     end
end





% beta in the likelihood 
if model.optimiseBeta
   
   if ~isstruct(model.betaTransform)
       fhandle = str2func([model.betaTransform 'Transform']);
       betaParam = fhandle(model.beta, 'xtoa');
   else
      if isfield(model.betaTransform,'transformsettings') && ~isempty(model.betaTransform.transformsettings)
          fhandle = str2func([model.betaTransform.type 'Transform']);
          betaParam = fhandle(model.beta, 'xtoa', model.betaTransform.transformsettings);
      else
          error('vargplvmExtractParam: Invalid transform specified for beta.'); 
      end
   end   

   params = [params betaParam(:)'];% add beta into param
   
   
   params = [params reshape(model.W(samd,:),1,[])];% add beta into param
   
%    if returnNames
%      for i = 1:length(betaParam)
%        betaParamNames{i} = ['Beta ' num2str(i)];
%      end
%      names = {names{:}, betaParamNames{:}};
%    end
end
% Kernel parameters  
% if returnNames
%     
%   [kernParams_v, kernParamNames_v] = kernExtractParam(model.kern_v); 
%   [kernParams_u, kernParamNames_u] = kernExtractParam(model.kern_u); 
%   
%   for i = 1:length(kernParamNames_v)
%     kernParamNames_v{i} = ['Kernel, ' kernParamNames_v{i}];
%   end
%   
%   for i = 1:length(kernParamNames_u)
%     kernParamNames_u{i} = ['Kernel, ' kernParamNames_u{i}];
%   end
%   
%   names = {names{:}, kernParamNames_u{:} , kernParamNames_v{:}};
% else
  kernParams_v = kernExtractParam(model.kern_v);
  kernParams_u = kernExtractParam(model.kern_u);
  
% end
params = [params kernParams_v kernParams_u];
