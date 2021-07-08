function model = cgpdsAddDynamics(model, type, varargin)

% cgpdsAddDynamics Add a dynamics structure to the model.

%	Copyright (c) 2010-2011 Andreas C. Damianou
%	Copyright (c) 2010-2011 Michalis K. Titsias
%	Copyright (c) 2010-2011 Neil D. Lawrence
% 	vargplvmAddDynamics.m SVN version 1748
% 	last update 2011-11-17T19:40:02.276622Z
type = [type 'Dynamics'];
samd = varargin{6};
samr = varargin{7};
model.dynamics = modelCreate(type, model.q, model.q, model.X, varargin{:}); % e.g. vargpTimeDynamicsCreate
params = cgpdsExtractParam(model,samd,samr);
model.dynamics.nParams = length(modelExtractParam(model.dynamics,samd,samr));
model.nParams = model.nParams + model.dynamics.nParams;
model = cgpdsExpandParam(model, params,samd,samr);

med = median(model.vardist.covars(:));
if med < 0.3 || med > 0.7
   fprintf('WARNING: Median value of variational covariances is %1.2f.\n', med);
end
