function  [h, compUpXP] =  lfmComputeH4JV(gamma1_p, gamma1_m, sigma2, t1, ...
    preFactor, preExp, mode)

% LFMCOMPUTEH4JV Helper function for computing part of the LFMJV kernel.
%
%	Description:
%
%	H = LFMCOMPUTEH4JV(GAMMA1, GAMMA2, SIGMA2, T1, T2, MODE) computes a
%	portion of the LFMJV kernel.
%	 Returns:
%	  H - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA1 - Gamma value for first system.
%	  GAMMA2 - Gamma value for second system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  MODE - indicates in which way the vectors t1 and t2 must be
%	   transposed
%	
%
%	See also
%	LFMCOMPUTEH4.M, LFMCOMPUTEH4AV.M


%	Copyright (c) 2010 Mauricio Alvarez
% 	lfmComputeH4JV.m SVN version 1519
% 	last update 2011-07-22T13:09:00.000000Z

if mode==0
    if nargout > 1
        compUpXP{1} = lfmjpComputeUpsilonVector(gamma1_p,sigma2, t1);
        compUpXP{2} = lfmjpComputeUpsilonVector(gamma1_m,sigma2, t1);
        h =  compUpXP{1}*( preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
            + compUpXP{2}*( preExp(:,1)/preFactor(4)- preExp(:,2)/preFactor(3)).';
    else
        h =  lfmjpComputeUpsilonVector(gamma1_p,sigma2, t1)*(preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
            + lfmjpComputeUpsilonVector(gamma1_m,sigma2, t1)*(preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
    end
else
    if nargout > 1
        compUpXP{1} = lfmvpComputeUpsilonVector(gamma1_p,sigma2, t1);
        compUpXP{2} = lfmvpComputeUpsilonVector(gamma1_m,sigma2, t1);
        h =  compUpXP{1}*(preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
            + compUpXP{2}*(preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
    else
        h =  lfmvpComputeUpsilonVector(gamma1_p,sigma2, t1)*(preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
            + lfmvpComputeUpsilonVector(gamma1_m,sigma2, t1)*(preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
    end
end