function [gKern, gVarmeans, gVarcovars] = biasVardistPsi0Gradient(biaskern, vardist, covGrad)

% BIASVARDISTPSI0GRADIENT one line description
%
%	Description:
%
%	[GKERN, GVARMEANS, GVARCOVARS, GIND] =
%	BIASVARDISTPSI0GRADIENT(BIASKERN, VARDIST, Z, COVGRAD) description
%	here.
%	 Returns:
%	  GKERN - 
%	  GVARMEANS - 
%	  GVARCOVARS - 
%	  GIND - 
%	 Arguments:
%	  BIASKERN - the kernel structure associated with the bias kernel.
%	  VARDIST - 
%	  Z - 
%	  COVGRAD - 
%	
%	
%	
%	
%
%	See also
%	OTHERS


%	Copyright (c) 2009 Michalis K. Titsias
% 	biasVardistPsi0Gradient.m SVN version 583
% 	last update 2009-11-08T13:07:32.000000Z

gKern = covGrad*vardist.numData;
 
gVarmeans = zeros(1,prod(size(vardist.means))); 
gVarcovars = zeros(1,prod(size(vardist.means))); 