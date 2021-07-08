function [gKern, gVarmeans, gVarcovars] = linard2VardistPsi0Gradient(linard2Kern, vardist, covGrad)

% LINARD2VARDISTPSI0GRADIENT description.
%
%	Description:
%	[gKern, gVarmeans, gVarcovars] = linard2VardistPsi0Gradient(linard2Kern, vardist, covGrad)
%% 	linard2VardistPsi0Gradient.m SVN version 583
% 	last update 2011-07-04T19:08:50.534563Z
  
A = linard2Kern.inputScales;
gKern = covGrad*sum((vardist.means.*vardist.means) + vardist.covars,1); 
 
gVarmeans = 2*(vardist.means*sparse(diag(A))); 
%gVarmeans1 = 2*(repmat(A,size(vardist.means,1),1).*vardist.means); 

gVarcovars = ones(size(vardist.means,1),1)*A; 

gVarmeans = covGrad*gVarmeans(:)'; 
gVarcovars = covGrad*gVarcovars(:)';


