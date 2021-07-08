function Psi0 = conbiasVardistPsi0Compute(conbiaskern, vardist)

% CONBIASVARDISTPSI0COMPUTE one line description
%
%	Description:
%
%	PSI0 = CONBIASVARDISTPSI0COMPUTE(CONBIASKERN, VARDIST) description
%	 Returns:
%	  PSI0 - conbias kern 计算得出的Psi0
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


%% 为了计算每一维的独立方差
if isfield(vardist,'chooseD')&& vardist.chooseD 
    conbiaskern.outputDimension=1;
end
%%
N = size(vardist.means,1);
D = conbiaskern.outputDimension;


%% 加入缺失数据情况
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
Nexist=sum(sum(vardist.map,1));% 存在观测数据的个数
Psi0 =Nexist*conbiaskern.variance; 

    
