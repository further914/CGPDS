function [K, P] = biasVardistPsi1Compute(biaskern, vardist, Z)

% BIASVARDISTPSI1COMPUTE one line description
%
%	Description:
%
%	[K, P] = BIASVARDISTPSI1COMPUTE(BIASKERN, VARDIST, Z) description
%	 Returns:
%	  K - description
%	  P - description
%	 Arguments:
%	  BIASKERN - the kernel structure associated with the white kernel.
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
% 	biasVardistPsi1Compute.m SVN version 583
% 	last update 2009-11-08T13:07:32.000000Z

K = repmat(biaskern.variance,size(vardist.means,1),size(Z,1));

P = [];