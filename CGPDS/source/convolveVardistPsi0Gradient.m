function [gKern, gVarmeans, gVarcovars] = convolveVardistPsi0Gradient(convolveKern, vardist, covGrad)

% RBFARD2VARDISTPSI0GRADIENT Description
%
%	Description:
% 	last update 2013-09-17
%   Copyright (c) ZhaoJing
[N,Q]  = size(vardist.means);
D = convolveKern.outputDimension;
sigma=repmat(convolveKern.Lambda_k,[D,1])+2*convolveKern.P_d;


% S_d的梯度
gS_d=2*N*convolveKern.S.*prod(sqrt(convolveKern.Lambda_k),2)./prod(sqrt(sigma),2)*covGrad;

%
[k0,k0part] =convolveVardistPsi0Compute(convolveKern, vardist);

% P_d的梯度 
gP_d=-repmat(k0part,[1,Q])./sigma*covGrad;

% Lambda的梯度 
gLambda=0.5*sum(repmat(k0part,[1,Q]).*(repmat(convolveKern.Lambda_k.^-1,[D,1])-sigma.^-1),1)*covGrad;


gKern=[gS_d' gP_d(:)' gLambda];
gVarmeans = zeros(1,numel(vardist.means)); 
gVarcovars = zeros(1,numel(vardist.means)); 