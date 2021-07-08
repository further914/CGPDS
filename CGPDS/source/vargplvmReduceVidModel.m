
% VARGPLVMREDUCEVIDMODEL Take a model computed on a video dataset and return a model which is computed on
%
%	Description:
%	the same video but with lower resolution
%	DESC Take a model computed on a video dataset and return a model which is computed on
%	the same video but with lower resolution.
%
%	See also
%	VARGPLVMREDUCEVIDEO.M
% 	vargplvmReduceVidModel.m SVN version 1430
% 	last update 2011-06-12T16:31:32.000000Z
% varargin may include "factor1" and "factor2" optionally
function [model, newHeight, newWidth] = vargplvmReduceVidModel(model, height, width, varargin)

[model.m, newHeight, newWidth] = vargplvmReduceVideo(model.m, height, width, varargin{:});
[model.scale, newHeight, newWidth] = vargplvmReduceVideo(model.scale, height, width, varargin{:});
[model.bias, newHeight, newWidth] = vargplvmReduceVideo(model.bias, height, width, varargin{:});
model.d = newHeight * newWidth;

if ~isempty(model.y)
    [model.y, newHeight, newWidth] = vargplvmReduceVideo(model.y, height, width, varargin{:});
end

model.TrYY = sum(sum(model.m .* model.m));
model.P = model.P1 * (model.Psi1' * model.m);
model.B = model.P1' * model.P;

params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);
