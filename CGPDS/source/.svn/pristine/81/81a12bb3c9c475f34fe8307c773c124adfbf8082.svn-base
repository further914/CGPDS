function [muqOrig SqOrig] = modelPriorKernGrad(dynModel)

% MODELPRIORKERNGRAD
%
%	Description:
%	see vargpTimeDynamicsPriorKernGrad. This function here is just a wrapper.
% 	modelPriorKernGrad.m SVN version 1430
% 	last update 2011-06-12T16:31:31.000000Z
fhandle = str2func([dynModel.type 'PriorKernGrad']);
[muqOrig SqOrig] = fhandle(dynModel);

% if isfield(model, 'paramGroups')
%  g = g*model.paramGroups;
% end
