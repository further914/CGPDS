function [K, P] = whiteVardistPsi1Compute(whitekern, vardist, Z)

% WHITEVARDISTPSI1COMPUTE one line description
%
%	Description:
%
%	[K, P] = WHITEVARDISTPSI1COMPUTE(WHITEKERN, VARDIST, Z) description
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
% 	whiteVardistPsi1Compute.m SVN version 583
% 	last update 2009-11-08T13:07:35.000000Z



  K = zeros(size(vardist.means,1), size(Z,1));
  
  P = [];