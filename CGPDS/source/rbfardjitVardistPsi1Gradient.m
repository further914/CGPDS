function [gKern, gVarmeans, gVarcovars, gInd] = rbfardjitVardistPsi1Gradient(rbfard2Kern, vardist, Z, covGrad)

% RBFARDJITVARDISTPSI1GRADIENT description.
%
%	Description:
%	[gKern, gVarmeans, gVarcovars, gInd] = rbfardjitVardistPsi1Gradient(rbfard2Kern, vardist, Z, covGrad)
%% 	rbfardjitVardistPsi1Gradient.m SVN version 1460
% 	last update 2011-07-04T19:08:51.585185Z
  
 [gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi1Gradient(rbfard2Kern, vardist, Z, covGrad);
