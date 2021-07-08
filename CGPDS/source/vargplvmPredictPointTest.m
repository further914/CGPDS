function [X, varX] = vargplvmPredictPointTest(dynModel, t_star)

% VARGPLVMPREDICTPOINT Predict the postions of a number of latent points.
%
%	Description:
%	[X, varX] = vargplvmPredictPoint(dynModel, t_star)
%% 	vargplvmPredictPoint.m SVN version 1312
% 	last update 2011-04-18T11:54:45.000000Z
% COPYRIGHT: Michalis Titsias, Neil Lawrence, Andreas Damianou 2011
%
% SEEALSO : vargplvmOptimisePoint
N = size(dynModel.t_star,1);
N_star = size(t_star,1); % number of test points
K_ss = kernDiagCompute(dynModel.kern, t_star);
K_star = kernCompute(dynModel.kern, dynModel.t_star, t_star);
K_t = kernCompute(dynModel.kern, dynModel.t_star,  dynModel.t_star);
X = K_star' * dynModel.vardist.means(end-N+1:end,:); % mean

varX = zeros(N_star, dynModel.q); % initialize variances
for q=1:dynModel.q
    invLambda = 1./dynModel.vardist.covars(end-N+1:end,q); 
    Lq = chol(K_t + diag(invLambda))'; 
    vq = Lq \ K_star;
    varX(:,q) = K_ss - sum(vq .* vq)';
end
