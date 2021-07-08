function [K, P] = conbiasVardistPsi1Compute(conbiaskern, vardist, Z)

% CONBIASVARDISTPSI1COMPUTE one line description
%	Description:



%	Copyright (c) 2013 ZhaoJing
% 	last update 2013-12-24


P = [];
%% 为了计算每一维的独立方差
if isfield(vardist,'chooseD')&& vardist.chooseD 
    conbiaskern.outputDimension=1;
end
%%

M = size(Z,1); 
D = conbiaskern.outputDimension;

N  = size(vardist.means,1);
%% 加入缺失数据情况
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
% K = repmat(conbiaskern.variance,N*D,size(Z,1));
K=repmat(reshape(vardist.map,[],1),1,M)*conbiaskern.variance;
