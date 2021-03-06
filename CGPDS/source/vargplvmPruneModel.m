function mm = vargplvmPruneModel(model, onlyData)

% VARGPLVMPRUNEMODEL Prune a var-GPLVM model.
%
%	Description:
%
%	MM = VARGPLVMPRUNEMODEL(MODEL, ONLYDATA) prunes a VAR-GPLVM model by
%	removing some fields which can later be reconstructed based on what
%	is being kept. Used when storing a model.
%	 Returns:
%	  MM - the variational GP-LVM model after being pruned
%	 Arguments:
%	  MODEL - the model to be pruned
%	  ONLYDATA - only prune the data parts. Useful when saving a model
%	   which is updated after predictions.
%	
%
%	See also
%	VARGPLVMREDUCEMODEL, VARGPLVMRESTOREPRUNEDMODEL


%	Copyright (c) Michalis Titsias, Neil Lawrence, 2011 Andreas Damianou
% 	vargplvmPruneModel.m SVN version 1687
% 	last update 2011-11-15T01:51:25.827041Z

if exist('onlyData') && onlyData
    model.m = [];
    model.y = [];
   if isfield(model, 'mOrig')
        model.mOrig = [];
   end
   if isfield(model, 'scale')
       model.scale = [];
   end
   if isfield(model, 'bias')
       model.bias = [];
   end
   model.P = [];
   model.B = [];
   mm = model;
   return
end


fieldsToKeep = ...
    {'type','approx','learnScales', 'optimiseBeta','betaTransform','q','d','N','optimiser','nPrivateParams','learnSigmaf',...
    'bias','scale','X','DgtN', 'kern_u','kern_v','k','fixInducing','inducingIndices','X_u','X_v','beta','prior','vardist','numParams', ...
    'nParams','K_uu','K_vv','learnBeta', 'dynamics', 'iters', 'date', 'fixedBetaIters', 'dataSetInfo','id','initVardist'}';

mm=[];
for i=1:length(fieldsToKeep)
    if isfield(model, fieldsToKeep{i})
        f = getfield(model,fieldsToKeep{i});
        mm=setfield(mm,fieldsToKeep{i},f);
    else
        %if ~strcmp(fieldsToKeep{i}, 'inducingIndices') && ~strcmp(fieldsToKeep{i}, 'iters') && ~strcmp(fieldsToKeep{i}, 'fixedBetaIters')
        %    fprintf(['??? Field ' fieldsToKeep{i} ' was missing from model! This field was skipped...\n']);
        %end
    end
end

