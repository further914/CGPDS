function model = vargpTimeDynamicsExpandParam(model, params, samd, samr)

% VARGPTIMEDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP time dynamics.
%
%	Description:
%	
%	Description:
%	
%	MODEL = VARGPTIMEDYNAMICSEXPANDPARAM(MODEL, PARAMS) takes the given
%	vector of parameters and places them in the model structure, it then
%	updates any stored representations that are dependent on those
%	parameters, for example kernel matrices etc..
%	Returns:
%	MODEL - a returned model structure containing the updated
%	parameters.
%	Arguments:
%	MODEL - the model structure for which parameters are to be
%	updated.
%	PARAMS - a vector of parameters for placing in the model
%	structure.
%	
%	
%	See also
%	vargpTimeDynamicsExtractParam, vargpTimeDynamicsCreate,
%	modelExpandParam


%	Copyright (c) 2010-2011 Andreas C. Damianou
%	Copyright (c) 2010-2011 Michalis K. Titsias
%	Copyright (c) 2010-2011 Neil D. Lawrence
% 	vargpTimeDynamicsExpandParam.m SVN version 1412
% 	last update 2011-07-04T19:08:51.784563Z
% Parameters are received in the following order: (notation: % parameter{size})
% [dynamicsVardistParams{dynamics.vardist.nParams}        % mu_bar, lambda
%  dynamics.kernParams{dynamics.kern.nParams}]            % sigmaRBF, lt,  sigmaWhite


% variational parameters (means and covariances) reparametrized
startVal = 1;

% if model.Rtr == 1
    endVal = model.vardist.nParams; 
% elseif model.Rtr > 1
%     partN = 0;
%     for rr = 1 : length(samr)
%         r = samr(rr);
%         partN = partN + model.seqTrain(r+1)-model.seqTrain(r);
%     end
%     endVal = partN * model.q*2;
% end





model.vardist = modelExpandParam(model.vardist, params(startVal:endVal),samd,samr); 
% dynamic model's kernel's parameters
startVal = endVal+1;
endVal = endVal + model.kern.nParams;

if isempty(find(isnan(params(startVal:endVal))))&&isempty(find(isinf(params(startVal:endVal))))
    model.kern = kernExpandParam(model.kern, params(startVal:endVal));
end


    %%%% DEBUG_
    %fprintf(1,'In vargpTimeDynamicsExpandParam, params=%d %d %d\n', model.kern.comp{1}.inverseWidth, model.kern.comp{1}.variance,model.kern.comp{2}.variance );
    %%% _DEBUG

model.nParams = endVal;


