function [K, outKern, sumKern, Kgvar] = rbfardjitVardistPsi2Compute(rbfardKern, vardist, Z)

% RBFARDJITVARDISTPSI2COMPUTE one line description
%
%	Description:
%
%	[K, OUTKERN, SUMKERN, KGVAR] =
%	RBFARDJITVARDISTPSI2COMPUTE(RBFARDKERN, VARDIST, Z) description
%	 Returns:
%	  K - description
%	  OUTKERN - description
%	  SUMKERN - description
%	  KGVAR - description
%	 Arguments:
%	  RBFARDKERN - the kernel structure associated with the rbfard2
%	   kernel.
%	  VARDIST - description
%	  Z - description
%	
%	
%	
%	
%	
%
%	See also
%	OTHERS


%	Copyright (c) 2009 Michalis K. Titsias
%	Copyright (c) 2009 Neil D. Lawrence
% 	rbfardjitVardistPsi2Compute.m SVN version 1460
% 	last update 2011-07-04T19:08:51.625186Z

[K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2Compute(rbfardKern, vardist, Z);