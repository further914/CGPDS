function [Psi1, P] = kernVardistPsi1Compute(kern, vardist, Z)

% KERNVARDISTPSI1COMPUTE description.
%
%	Description:
%	[Psi1, P] = kernVardistPsi1Compute(kern, vardist, Z)
%% 	kernVardistPsi1Compute.m SVN version 583
% 	last update 2009-11-08T13:07:33.000000Z


if ~strcmp(kern.type,'cmpnd')
  %  
  fhandle = str2func([kern.type 'VardistPsi1Compute']);
  [Psi1 P] = fhandle(kern, vardist, Z);
  %  
else % the kernel is cmpnd
  %     
  fhandle = str2func([kern.comp{1}.type 'VardistPsi1Compute']);
  [Psi1 P] = fhandle(kern.comp{1}, vardist, Z);
  %
  for i = 2:length(kern.comp)
      %
      fhandle = str2func([kern.comp{i}.type 'VardistPsi1Compute']);
      [Ptmp P] = fhandle(kern.comp{i}, vardist, Z);
      Psi1 = Psi1 + Ptmp;
      %
  end
  %
end
