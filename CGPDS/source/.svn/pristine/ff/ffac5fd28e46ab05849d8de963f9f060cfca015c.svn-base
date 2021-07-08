function [gKern, gVarmeans, gVarcovars, gInd] = biasVardistPsi1Gradient(biaskern, vardist, Z, covGrad)

% BIASVARDISTPSI1GRADIENT Compute gradient of bias variational PSI1.
%
%	Description:
%
%	[GKERN, GVARMEANS, GVARCOVARS, GIND] =
%	BIASVARDISTPSI1GRADIENT(BIASKERN, VARDIST, Z, COVGRAD) description
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
% 	biasVardistPsi1Gradient.m SVN version 583
% 	last update 2009-11-08T13:07:32.000000Z
  
gKern = sum(sum(ones(vardist.numData,size(Z,1)).*covGrad)); 

gVarmeans = zeros(1,prod(size(vardist.means))); 

gInd = zeros(1,prod(size(Z))); 

gVarcovars = zeros(1,prod(size(vardist.covars))); 
