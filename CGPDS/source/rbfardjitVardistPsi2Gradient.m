function [gKern, gVarmeans, gVarcovars, gInd] = rbfardjitVardistPsi2Gradient(rbfardKern, vardist, Z, covGrad)

% RBFARDJITVARDISTPSI2GRADIENT description.
%
%	Description:
%	[gKern, gVarmeans, gVarcovars, gInd] = rbfardjitVardistPsi2Gradient(rbfardKern, vardist, Z, covGrad)
%% 	rbfardjitVardistPsi2Gradient.m SVN version 1473
% 	last update 2011-07-04T19:08:51.724564Z
  
[gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi2Gradient(rbfardKern, vardist, Z, covGrad);
