function [gKern, gVarmeans, gVarcovars, gInd] = conbiasVardistPsi1Gradient(conbiaskern, vardist, Z, covGrad)

% CONBIASVARDISTPSI1GRADIENT Compute gradient of conbias variational PSI1.
%
%	Description:



%	Copyright (c) 2013 ZhaoJing
% 	last update 2013-12-24
D = conbiaskern.outputDimension;
M = size(Z,1);
N  = size(vardist.means,1);
%% 加入缺失数据情况
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end

Ktmp=repmat(reshape(vardist.map,[],1),1,M);% 0-1 matrix with size of ND*M
gKern = sum(sum(Ktmp.*covGrad)); 

gVarmeans = zeros(1,prod(size(vardist.means))); 

gInd = zeros(1,prod(size(Z))); 

gVarcovars = zeros(1,prod(size(vardist.covars))); 
