function [gKern, gVarmeans, gVarcovars, gInd] = whiteVardistPsi2Gradient(whitekern, vardist, Z, covGrad)

% WHITEVARDISTPSI2GRADIENT Compute gradient of white variational PSI2.
%
%	Description:
%
%	[GKERN, GVARMEANS, GVARCOVARS, GIND] =
%	WHITEVARDISTPSI2GRADIENT(WHITEKERN, VARDIST, Z, COVGRAD) description
%	here.
%	 Returns:
%	  GKERN - 
%	  GVARMEANS - 
%	  GVARCOVARS - 
%	  GIND - 
%	 Arguments:
%	  WHITEKERN - the kernel structure associated with the white kernel.
%	  VARDIST - 
%	  Z - 
%	  COVGRAD - 
%	
%	
%
%	See also
%	


%	Copyright (c) 2009 Michalis K. Titsias
% 	whiteVardistPsi2Gradient.m SVN version 583
% 	last update 2009-11-08T13:07:35.000000Z
  
  gKern = 0;
  gVarmeans = zeros(1,prod(size(vardist.means))); 
  gInd = zeros(1,prod(size(Z))); 
  gVarcovars = zeros(1,prod(size(vardist.covars))); 

