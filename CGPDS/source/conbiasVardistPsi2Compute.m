function [Psi2 P] = conbiasVardistPsi2Compute(conbiaskern, vardist, Z)

% CONBIASVARDISTPSI2COMPUTE one line description



%	Copyright (c) 2013 ZhaoJing
% 	last update 2013-12-24
%% 为了计算每一维的独立方差
if isfield(vardist,'chooseD')&& vardist.chooseD 
    conbiaskern.outputDimension=1;
end
%%
D = conbiaskern.outputDimension;

N  = size(vardist.means,1);
%% 加入缺失数据情况
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
Nexist=sum(sum(vardist.map,1));% 存在观测数据的个数
Psi2 = repmat(Nexist*(conbiaskern.variance^2),size(Z,1),size(Z,1)); 
P = [];




