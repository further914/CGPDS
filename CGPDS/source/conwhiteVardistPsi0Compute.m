function Psi0 = conwhiteVardistPsi0Compute(conwhitekern, vardist)

% CONWHITEVARDISTPSI0COMPUTE one line description
%
%	Description:


%	Copyright (c) 2009 Michalis K. Titsias
%	Copyright (c) 2013 ZhaoJing
% 	last update 2013-12-24

% the "white" kernel only affects the K_uu matrix (jitter inducing variables)"
% that's why the psi0_white = 0 
Psi0 = 0; % vardist.numData*whitekern.variance; 
%Psi0 =  vardist.numData*whitekern.variance;
