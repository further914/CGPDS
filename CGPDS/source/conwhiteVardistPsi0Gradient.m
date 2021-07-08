function [gKern, gVarmeans, gVarcovars] = conwhiteVardistPsi0Gradient(conwhitekern, vardist, covGrad)

% CONWHITEVARDISTPSI0GRADIENT one line description
%
%	Description:
%
%	[GKERN, GVARMEANS, GVARCOVARS, GIND] =
%	CONWHITEVARDISTPSI0GRADIENT(WHITEKERN, VARDIST, Z, COVGRAD) description
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
%	
%
%	See also
%	OTHERS


%	Copyright (c) 2009 Michalis K. Titsias
%	Copyright (c) 2013 ZhaoJing
% 	whiteVardistPsi0Gradient.m SVN version 583
% 	last update 2009-11-08T13:07:35.000000Z
% 	last update 2013-12-24
% the "white" kernel only affects the K_uu matrix (jitter inducing variables)"
% that's why the gKern = 0 
gKern = 0; % covGrad*vardist.numData;
%gKern = covGrad*vardist.numData;
gVarmeans = zeros(1,prod(size(vardist.means))); 
gVarcovars = zeros(1,prod(size(vardist.means))); 