function model = vargplvmRestorePrunedModel(model, Ytr, onlyData, options)

% VARGPLVMRESTOREPRUNEDMODEL Restore a pruned var-GPLVM model.
%
%	Description:
%
%	MODEL = VARGPLVMRESTOREPRUNEDMODEL(MODEL, YTR, ONLYDATA) restores a
%	vargplvm model which has been pruned and it brings it in the same
%	state that it was before pruning.
%	 Returns:
%	  MODEL - the variational GP-LVM model after being restored
%	 Arguments:
%	  MODEL - the model to be restored
%	  YTR - the training data
%	  ONLYDATA - only pruned the data parts. Useful when saving a model
%	   which is updated after predictions.
%	
%
%	See also
%	VARGPLVMREDUCEMODEL, VARGPLVMPRUNEMODEL


%	Copyright (c) Michalis Titsias, Neil Lawrence, 2011 Andreas Damianou
% 	vargplvmRestorePrunedModel.m SVN version 1687
% 	last update 2011-11-15T01:51:25.836341Z


% TODO: We could also prune bias and scale (big matrices) but then the user
% should also give the rule for making them again.

if exist('onlyData') && onlyData
    %  model.mOrig = model.m;
    model.bias = mean(Ytr); % Default, has to be changed later if it was different
    
    if (nargin > 3) && ~isempty(options) && isfield(options,'scale2var1')
        if(options.scale2var1)
            model.scale = std(Ytr);
            model.scale(find(model.scale==0)) = 1;
            if(model.learnScales)
                warning('Both learn scales and scale2var1 set for GP');
            end
            if(isfield(options, 'scaleVal'))
                warning('Both scale2var1 and scaleVal set for GP');
            end
        end
    elseif  (nargin > 3) && ~isempty(options) && isfield(options, 'scaleVal')
        model.scale = repmat(options.scaleVal, 1, size(Ytr,2));
    else
        model.scale = ones(1,size(Ytr,2));
    end
end

model.y = Ytr;
model.m= gpComputeM(model);
model.y=reshape(model.y,model.d*model.N,1);
model.m=reshape(model.m,model.d*model.N,1);
if model.DgtN
    model.mOrig = model.m;
    YYT = model.m * model.m'; % NxN
    % Replace data with the cholesky of Y*Y'.Same effect, since Y only appears as Y*Y'.
    %%% model.m = chol(YYT, 'lower');  %%% Put a switch here!!!!
    [U S V]=svd(YYT);
    model.m=U*sqrt(abs(S));
    model.TrYY = sum(diag(YYT)); % scalar
else
    model.TrYY = sum(sum(model.m .* model.m));
end

if exist('onlyData') && onlyData
    model.P = model.P1 * (model.Psi1' * model.m);
    model.B = model.P1' * model.P;
    model = orderfields(model);
    return
end

params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);


model = orderfields(model);

