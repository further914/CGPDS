function f = cgpdsSeqDynObjective(x, model, y,samd,samr)

% cgpdsSeqDynObjective Wrapper function for objective of a group of points in latent space and the output locations..

%	Copyright (c) 2011 Michalis K. Titsias and Neil D. Lawrence
% 	vargplvmSeqDynObjective.m SVN version 1312
% 	last update 2011-04-18T11:54:44.000000Z

% % this is doing the expand 
% x = reshape(x, model.N+size(y,1), model.dynamics.q*2);
% xtrain = x(1:model.N,:);
% xtest = x(model.N+1:end,:);
% model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
% vardistx = vardistExpandParam(model.vardistx, xtest);
% % end of expand 
vardistx = model.vardistx;
vardistx = vardistExpandParam(vardistx, x);
fhandle = str2func([model.type 'SeqDynLogLikelihood']);
f = -fhandle(model, vardistx, y,samd,samr);

%f = -TEMPvargplvmPointLogLikelihoodSTATIC(model, vardistx, y);