function k = lfmaKernDiagCompute(lfmKern, t)

% LFMAKERNDIAGCOMPUTE Compute diagonal of LFMAXLFMA kernel.
%
%	Description:
%
%	K = LFMAKERNDIAGCOMPUTE(KERN, T) computes the diagonal of the kernel
%	matrix for acceleration - acceleration in the switching dynamical
%	latent force kernel.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T - input data matrix in the form of a design matrix.


%	Copyright (c) 2010 Mauricio Alvarez
% 	lfmaKernDiagCompute.m SVN version 809
% 	last update 2010-05-28T06:01:33.000000Z



% Get length scale out.
sigma2 = 2/lfmKern.inverseWidth;
sigma = sqrt(sigma2);

% Parameters of the kernel
alpha(1) = lfmKern.damper./(2*lfmKern.mass);
omega(1) = sqrt(lfmKern.spring./lfmKern.mass - alpha(1)*alpha(1));

% Precomputations to increase the speed
preExp1 = zeros(length(t),2);
gamma1_p = alpha(1) + j*omega(1);
gamma1_m = alpha(1) - j*omega(1);
preGamma(1) = gamma1_p + gamma1_p;
preGamma(2) = gamma1_p + gamma1_m;
preGamma(3) = gamma1_m + gamma1_p;
preGamma(4) = gamma1_m + gamma1_m;
preConst = 1./preGamma;
preFactors(1) = preConst(2) - preConst(1);
preFactors(2) = preConst(3) - preConst(4);
preFactors(3) = preConst(3) - preConst(1);
preFactors(4) = preConst(2) - preConst(4);
preExp1(:,1) = (gamma1_p^2)*exp(-gamma1_p*t);
preExp1(:,2) = (gamma1_m^2)*exp(-gamma1_m*t);
% Actual computation of the kernel
sk =lfmDiagComputeH3AA(gamma1_p, gamma1_m, sigma2, t, preFactors([1 2])) + ...
    lfmDiagComputeH3AA(gamma1_p, gamma1_m, sigma2, t, preFactors([3 4])) + ...
    lfmDiagComputeH4AA(gamma1_p, gamma1_m, sigma2, t, preGamma([1 2 4 3]), preExp1) + ...
    lfmDiagComputeH4AA(gamma1_p, gamma1_m, sigma2, t, preGamma([1 3 4 2]), preExp1);

if lfmKern1.isNormalised
    k0 = kern.sensitivity^2/(8*sqrt(2)*kern.mass^2*omega^2);
else
    k0 = sqrt(pi)*sigma*kern.sensitivity^2/(8*kern.mass^2*omega^2);
end
k = k0*sk;


function h = lfmDiagComputeH3AA(gamma1_p, gamma1_m, sigma2, t,preFactor)

h = preFactor(1)*lfmaaComputeUpsilonDiagVector(gamma1_p,sigma2, t,0) ...
    + preFactor(2)*lfmaaComputeUpsilonDiagVector(gamma1_m,sigma2, t,0);


function h = lfmDiagComputeH4AA(gamma1_p, gamma1_m, sigma2, t, ...
   preFactor, preExp)


    h =  lfmapComputeUpsilonVector(gamma1_p,sigma2, t).*(preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)) ...
        + lfmapComputeUpsilonVector(gamma1_m,sigma2, t).*(preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4));



