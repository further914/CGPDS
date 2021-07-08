function Psi0 = linard2VardistPsi0Compute(linard2kern, vardist)

% LINARD2VARDISTPSI0COMPUTE Description
%
%	Description:
%	Psi0 = linard2VardistPsi0Compute(linard2kern, vardist)
%% 	linard2VardistPsi0Compute.m SVN version 583
% 	last update 2011-07-04T19:08:50.324563Z
  
A = linard2kern.inputScales;
Psi0 = sum(A.*sum((vardist.means.*vardist.means) + vardist.covars,1));

%Psi0covs = sum(A.*sum(vardist.covars,1)); 

%psi00 = kernCompute(linard2kern,vardist.means);
%ok = 0;
%for n=1:size(vardist.means,1)
%    ok = ok + trace(diag(A)*(vardist.means(n,:)'*vardist.means(n,:))) + trace(diag(A)*diag(vardist.covars(n,:))); 
%end
%ok
%pause
%Psi0 = Psi0means + Psi0covs;



