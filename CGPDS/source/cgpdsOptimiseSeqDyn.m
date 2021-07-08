function [X varX] = cgpdsOptimiseSeqDyn(model, vardistx, y, display, iters,samd,samr)

% cgpdsOptimiseSeqDyn Optimise vraitaional means and covar for testing data


if nargin < 5
  iters = 2000;
  %if nargin < 5
    display = true;
  %end
end

options = optOptions;
if display
  options(1) = 1;
  %options(9) = 1;
end
options(14) = iters;


if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('scg');
end


% % augment the training and testing variational distributions
% vardist = vardistCreate(zeros(model.N+size(y,1), model.q), model.q, 'gaussian');
% vardist.means = [model.dynamics.vardist.means; vardistx.means];
% vardist.covars = [model.dynamics.vardist.covars; vardistx.covars]; 
% vardist.numData = size(vardist.means,1);
% vardist.nParams = 2*prod(size(vardist.means));
x = modelExtractParam(vardistx,samd,samr);

% NETLAB style optimization.
x = optim('cgpdsSeqDynObjective', x,  options, ...
    'cgpdsSeqDynGradient', model, y,samd,samr);

% now separate the variational disribution into the training part and the
% testing part and update the original training model (only with the new training 
% variational distribution) and the test variational distribution
% this is doing the expand 

% x = reshape(x, vardist.numData, model.dynamics.q*2);
% xtrain = x(1:model.N,:);
% xtest = x(model.N+1:end,:);
% model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
% vardistx = vardistExpandParam(model.vardistx, xtest);
% % end of expand 

vardistx = vardistExpandParam(vardistx,x);
X = vardistx.means;
varX = vardistx.covars;
