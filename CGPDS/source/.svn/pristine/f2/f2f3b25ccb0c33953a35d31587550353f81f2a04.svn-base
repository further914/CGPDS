function [gKern, gVarmeans, gVarcovars, gInd] = conwhiteVardistPsi2Gradient(conwhitekern, vardist, Z, covGrad)

% CONWHITEVARDISTPSI2GRADIENT Compute gradient of white variational PSI2.
%
%	Description:
%
%	[GKERN, GVARMEANS, GVARCOVARS, GIND] =
%	CONWHITEVARDISTPSI2GRADIENT(CONWHITEKERN, VARDIST, Z, COVGRAD) description
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

%	Copyright (c) 2013 ZhaoJing
% 	conwhiteVardistPsi2Gradient.m SVN version 583
% 	last update 2009-11-08T13:07:35.000000Z
%   last update 2013-12-24
  
  gKern = 0;
  gVarmeans = zeros(1,prod(size(vardist.means))); 
  gInd = zeros(1,prod(size(Z))); 
  gVarcovars = zeros(1,prod(size(vardist.covars))); 

