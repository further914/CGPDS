function [Psi2, Pnobias, Psi1] = linard2biasVardistPsi2Compute(linardKern, biasKern, vardist, Z)

% LINARD2BIASVARDISTPSI2COMPUTE one line description
%
%	Description:
%
%	[PSI2, PNOBIAS, PSI1] = LINARD2BIASVARDISTPSI2COMPUTE(LINARD2KERN,
%	BIASKERN, VARDIST, Z) description
%	 Returns:
%	  PSI2 - description
%	  PNOBIAS - description
%	  PSI1 - description
%	 Arguments:
%	  LINARD2KERN - the kernel structure associated with the linard2
%	   kernel.
%	  BIASKERN - description
%	  VARDIST - description
%	  Z - description
%	
%	
%	
%	
%	
%
%	See also
%	OTHERS


%	Copyright (c) 2009 Michalis K. Titsias
%	Copyright (c) 2009 Neil D. Lawrence
% 	linard2biasVardistPsi2Compute.m SVN version 583
% 	last update 2011-07-04T19:08:50.744563Z


Psi1 = linard2VardistPsi1Compute(linardKern, vardist, Z); 

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

