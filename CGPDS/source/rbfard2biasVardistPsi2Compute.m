function [Psi2, Pnobias, Psi1] = rbfard2biasVardistPsi2Compute(rbfardKern, biasKern, vardist, Z)

% RBFARD2BIASVARDISTPSI2COMPUTE description.
%
%	Description:
%	[Psi2, Pnobias, Psi1] = rbfard2biasVardistPsi2Compute(rbfardKern, biasKern, vardist, Z)
%% 	rbfard2biasVardistPsi2Compute.m SVN version 583
% 	last update 2011-07-04T19:08:51.124563Z

Psi1 = rbfard2VardistPsi1Compute(rbfardKern, vardist, Z); 

sumPsi1 = sum(Psi1,1); 

Psi2 = ones(size(Z,1),1)*sumPsi1; 

Pnobias = Psi2 + Psi2'; 
Psi2 = biasKern.variance*Pnobias; 

% second naive way
%Psi22 = zeros(size(Z,1),size(Z,1));
%for j=1:size(Z,1)
%    for i=1:size(Z,1)
%        Psi22(j,i) = biasKern.variance*(sum(Psi1(:,j)) + sum(Psi1(:,i)));
%    end
%end
%sum(sum(abs(Psi2 - Psi22))) 

