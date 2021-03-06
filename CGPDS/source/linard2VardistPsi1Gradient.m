function [gKern, gVarmeans, gVarcovars, gInd] = linard2VardistPsi1Gradient(linard2Kern, vardist, Z, covGrad)

% LINARD2VARDISTPSI1GRADIENT description.
%
%	Description:
%	[gKern, gVarmeans, gVarcovars, gInd] = linard2VardistPsi1Gradient(linard2Kern, vardist, Z, covGrad)
%% 	linard2VardistPsi1Gradient.m SVN version 583
% 	last update 2011-07-04T19:08:50.584563Z

A = linard2Kern.inputScales;
%
% TYPICAL WAY
%for q=1:vardist.latentDimension
%   % 
%   gKern(q) = sum(sum((vardist.means(:,q)*(Z(:,q)')).*covGrad));
%   
%   gVarmeans(:,q) = A(q)*sum((ones(vardist.numData,1)*Z(:,q)').*covGrad,2);
%   gInd(:,q) = A(q)*sum((vardist.means(:,q)*ones(1,size(Z,1))).*covGrad,1)';
%   %
%end
%%% end of typical way 


% FAST WAY
AA = ones(size(vardist.means,1),1)*A; 
covVarm = covGrad'*vardist.means;
gKern = sum(covVarm.*Z,1); 
gVarmeans = AA.*(covGrad*Z);
AA = ones(size(Z,1),1)*A;
gInd = AA.*covVarm;

%sum(sum(abs(gKern1-gKern)))
%sum(sum(abs(gVarmeans1 - gVarmeans)))
%sum(sum(abs(gInd1 - gInd)))
%pause

gKern = gKern(:)';  
% gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
gVarmeans = gVarmeans(:)'; 

% gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
% this will unfold this matrix column-wise 
gInd = gInd(:)'; 

gVarcovars = zeros(1,prod(size(vardist.covars))); 
