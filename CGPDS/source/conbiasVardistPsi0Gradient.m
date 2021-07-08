function [gKern, gVarmeans, gVarcovars] = conbiasVardistPsi0Gradient(conbiaskern, vardist, covGrad)

% CONBIASVARDISTPSI0GRADIENT one line description
% 
%
%	Description:
%
%	[GKERN, GVARMEANS, GVARCOVARS] = CONBIASVARDISTPSI0GRADIENT(CONBIASKERN, VARDIST, COVGRAD) description
%	 Returns:
%	  [GKERN, GVARMEANS, GVARCOVARS] - 计算关于conbiaskern的Psi0的梯度 
%	 Arguments:
%	  CONBIASKERN - the kernel structure associated with the conbias kernel.可与convolved kern 搭配使用的
%	  VARDIST - description
%	
%	
%	
%	
%
%	See also
%	OTHERS

%	Copyright (c) 2013 ZhaoJing
% 	last update 2013-12-24

%%
N = size(vardist.means,1);
D = conbiaskern.outputDimension;


%% 加入缺失数据情况
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
Nexist=sum(sum(vardist.map,1));% 存在观测数据的个数
gKern = covGrad*Nexist;
 
gVarmeans = zeros(1,prod(size(vardist.means))); 
gVarcovars = zeros(1,prod(size(vardist.means))); 