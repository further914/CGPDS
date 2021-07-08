function [gKern, gVarmeans, gVarcovars] = rbfardjitVardistPsi0Gradient(rbfardjitKern, vardist, covGrad)

% RBFARDJITVARDISTPSI0GRADIENT Description
%
%	Description:
%	[gKern, gVarmeans, gVarcovars] = rbfardjitVardistPsi0Gradient(rbfardjitKern, vardist, covGrad)
%% 	rbfardjitVardistPsi0Gradient.m SVN version 1460
% 	last update 2011-07-04T19:08:51.494563Z

[gKern, gVarmeans, gVarcovars] = rbfard2VardistPsi0Gradient(rbfardjitKern, vardist, covGrad);

