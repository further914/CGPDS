function vardist = vardistExpandParam(vardist, params,samd, samr)

% VARDISTEXPANDPARAM Expand a parameter vector into a vardist structure.
%
%	Description:
%
%	VARDISTEXPANDPARAM(MODEL, PARAMS) takes an VARDIST structure and a
%	vector of parameters, and fills the structure with the given
%	parameters. Also performs any necessary precomputation for
%	likelihood and gradient computations, so can be computationally
%	intensive to call.
%	 Arguments:
%	  MODEL - the VARDIST structure to put the parameters in.
%	  PARAMS - parameter vector containing the parameters to put in the
%	   VARDIST structure.
%	
%	
%	
%
%	See also
%	VARDISTCREATE, VARDISTEXTRACTPARAM, MODELEXPANDPARAM


%	Copyright (c) 2009 Michalis K. Titsias
%	Copyright (c) 2009 Neil D. Lawrence
% 	vardistExpandParam.m SVN version 583
% 	last update 2009-11-08T13:07:34.000000Z

if ~isempty(vardist.transforms)
  for i = 1:length(vardist.transforms)
%     index = vardist.transforms(i).index;
    index = length(params)/2+1:length(params);
    fhandle = str2func([vardist.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'atox');
  end
end

% if model.Rtr == 1
    means = params(1:(vardist.numData*vardist.latentDimension));
    st = vardist.numData*vardist.latentDimension + 1;
    covs = params(st:end);
    vardist.means = reshape(means, vardist.numData, vardist.latentDimension);
    vardist.covars = reshape(covs, vardist.numData, vardist.latentDimension);
    
% elseif model.Rtr > 1
%     means = params(1:length(params)/2);
%     covs = params(1:length(params)/2+1:end);
%     vardist.means = reshape(means, [], vardist.latentDimension);
%     vardist.covars = reshape(covs, [], vardist.latentDimension);
%     
% end





