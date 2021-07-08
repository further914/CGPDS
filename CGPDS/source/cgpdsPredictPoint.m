function [X, varX] = cgpdsPredictPoint(dynModel, t_star)

% cgpdsPredictPoint Predict the postions of a number of latent points.
%
%	Description:
%	[X, varX] = cgpdsPredictPoint(dynModel, t_star)
%% 	vargplvmPredictPoint.m SVN version 1312
% 	last update 2011-04-18T11:54:45.000000Z
% COPYRIGHT: Michalis Titsias, Neil Lawrence, Andreas Damianou 2011
%
% SEEALSO : vargplvmOptimisePoint

N_star = size(t_star,1); % number of test points
K_ss = kernDiagCompute(dynModel.kern, t_star);
K_star = kernCompute(dynModel.kern, dynModel.t, t_star);

X = K_star' * dynModel.vardist.means; % mean

varX = zeros(N_star, dynModel.q); % initialize variances
for q=1:dynModel.q
    invLambda = 1./dynModel.vardist.covars(:,q); 
    Lq = chol(dynModel.Kt + diag(invLambda))'; 
    vq = Lq \ K_star;
    varX(:,q) = K_ss - sum(vq .* vq)';
end
