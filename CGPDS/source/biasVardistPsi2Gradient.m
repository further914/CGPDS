function [gKern, gVarmeans, gVarcovars, gInd] = biasVardistPsi2Gradient(biaskern, vardist, Z, covGrad)

% BIASVARDISTPSI2GRADIENT Compute gradient of bias variational PSI2.
%
%	Description:
%
%	[GKERN, GVARMEANS, GVARCOVARS, GIND] =
%	BIASVARDISTPSI2GRADIENT(BIASKERN, VARDIST, Z, COVGRAD) description
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
%	See also
%	


%	Copyright (c) 2009 Michalis K. Titsias
% 	biasVardistPsi2Gradient.m SVN version 583
% 	last update 2009-11-08T13:07:32.000000Z

gKern = (2*vardist.numData*biaskern.variance)*sum(sum(ones(size(Z,1),size(Z,1)).*covGrad)); 

gVarmeans = zeros(1,prod(size(vardist.means))); 

gInd = zeros(1,prod(size(Z))); 

gVarcovars = zeros(1,prod(size(vardist.covars))); 

  