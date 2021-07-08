function [K, P] = conwhiteVardistPsi1Compute(conwhitekern, vardist, Z)

% CONWHITEVARDISTPSI1COMPUTE one line description
%
%	Description:
%
%	[K, P] = CONWHITEVARDISTPSI1COMPUTE(WHITEKERN, VARDIST, Z) description
%	 Returns:
%	  K - description
%	  P - description
%	 Arguments:
%	  WHITEKERN - the kernel structure associated with the white kernel.
%	  VARDIST - description
%	  Z - description
%	
%	
%	
%	
%
%	See also
%	OTHERS


%	Copyright (c) 2009 Michalis K. Titsias
%	Copyright (c) 2013 ZhaoJing
% 	whiteVardistPsi1Compute.m SVN version 583
% 	last update 2009-11-08T13:07:35.000000Z
% 	last update 2013-12-24


%% 为了计算每一维的独立方差
if isfield(vardist,'chooseD')&& vardist.chooseD 
    conwhitekern.outputDimension=1;
end
%%

M = size(Z,1); 
D = conwhitekern.outputDimension;

N  = size(vardist.means,1);
  K = zeros(size(vardist.means,1)*conwhitekern.outputDimension, size(Z,1));
  
  P = [];