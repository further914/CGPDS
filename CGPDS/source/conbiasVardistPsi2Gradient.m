function [gKern, gVarmeans, gVarcovars, gInd] = conbiasVardistPsi2Gradient(conbiaskern, vardist, Z, covGrad)

% CONBIASVARDISTPSI2GRADIENT Compute gradient of bias variational PSI2.

%	
%	
%
%	See also
%	


%	Copyright (c) 2013 ZhaoJing
% 	last update 2013-12-24
D = conbiaskern.outputDimension;

N  = size(vardist.means,1);
%% 加入缺失数据情况
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
Nexist=sum(sum(vardist.map,1));% 存在观测数据的个数

gKern = (2*Nexist*conbiaskern.variance)*sum(sum(ones(size(Z,1),size(Z,1)).*covGrad)); 

gVarmeans = zeros(1,prod(size(vardist.means))); 

gInd = zeros(1,prod(size(Z))); 

gVarcovars = zeros(1,prod(size(vardist.covars))); 

