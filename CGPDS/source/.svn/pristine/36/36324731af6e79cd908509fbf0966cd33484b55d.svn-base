function [K, Knovar, argExp] = rbfardVardistPsi1Compute(rbfardKern, vardist, Z)

% RBFARDVARDISTPSI1COMPUTE one line description
%
%	Description:
%
%	[K, KNOVAR, ARGEXP] = RBFARDVARDISTPSI1COMPUTE(RBFARDKERN, VARDIST,
%	Z) description
%	 Returns:
%	  K - description
%	  KNOVAR - description
%	  ARGEXP - description
%	 Arguments:
%	  RBFARDKERN - the kernel structure associated with the rbfard
%	   kernel.
%	  VARDIST - description
%	  Z - description
%	
%	
%	
%	
%
%	See also
%	OTHERS


%	Copyright (c) 2009 Michalis K. Titsias
% 	rbfardVardistPsi1Compute.m SVN version 583
% 	last update 2009-11-08T13:07:34.000000Z

% variational means
N  = size(vardist.means,1);
%  inducing variables 
M = size(Z,1); 

A = rbfardKern.inverseWidth;
         
argExp = zeros(N,M); 
normfactor = ones(N,1);
for q=1:vardist.latentDimension
%
    S_q = vardist.covars(:,q);  
    normfactor = normfactor.*(A(q)*S_q + 1);
    Mu_q = vardist.means(:,q); 
    Z_q = Z(:,q)'; 
    distan = (repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1])).^2;
    size(distan)
    size(repmat(A(q)/(A(q)*S_q + 1), [1 M]))
    argExp = argExp + repmat(A(q)/(A(q)*S_q + 1), [1 M]).*distan;
%
end
Knovar = repmat(normfactor,[1 M]).*exp(-0.5*argExp); 
K = rbfardKern.variance.*Knovar; 


