function [params, names] = vardistExtractParam(vardist,samd,samr)

% VARDISTEXTRACTPARAM Extract a parameter vector from a vardist structure.
%
%	Description:
%
%	PARAMS = VARDISTEXTRACTPARAM(MODEL) extracts a parameter vector from
%	a given VARDIST structure.
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
%
%	See also
%	VARDISTCREATE, VARDISTEXPANDPARAM, MODELEXTRACTPARAM


%	Copyright (c) 2005, 2006 Neil D. Lawrence
% 	vardistExtractParam.m SVN version 583
% 	last update 2009-11-08T13:07:34.000000Z



%means = vardist.means'; 
%covs = vardist.covars';

% the variational means and diagonal covariances obtained COLUMN-WISE 

% if model.Rtr == 1
    params = [vardist.means(:)' vardist.covars(:)'];
% elseif model.Rtr>1
%     params = [];
%     for rr = 1:length(samr)
%         r = samr(rr);
%         subseq = model.seqTrain(r)+1:model.seqTrain(r+1);
%         params = [params reshape(vardist.means(subseq,:),1,[])];
%     end
%     for rr = 1:length(samr)
%         r = samr(rr);
%         subseq = model.seqTrain(r)+1:model.seqTrain(r+1);
%         params = [params reshape(vardist.covars(subseq,:),1,[])];
%     end
% end
    





% names
% if nargout > 1  
%     for i=1:size(vardist.means,1)
%     for j=1:size(vardist.means,2)
%         varmeanNames{i,j} = ['varmean(' num2str(i) ', ' num2str(j) ')'];
%     end
%     end
%     for i=1:size(vardist.means,1)
%     for j=1:size(vardist.means,2)
%         varcovNames{i,j} = ['varcov(' num2str(i) ', ' num2str(j) ')'];
%     end
%     end
%     names = {varmeanNames{:}, varcovNames{:}}; 
% end

% Check if parameters are being optimised in a transformed space.
if ~isempty(vardist.transforms)
  for i = 1:length(vardist.transforms)
    %     index = vardist.transforms(i).index;
    %expTransform????, expTransform to ensure covar is positive, the second
    %half need to be transformed
    index = length(params)/2+1:length(params);
    fhandle = str2func([vardist.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'xtoa');
  end
end
