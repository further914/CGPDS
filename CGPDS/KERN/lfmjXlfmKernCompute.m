function K = lfmjXlfmKernCompute(lfmKern1, lfmKern2, t1, t2)

% LFMJXLFMKERNCOMPUTE Jolt and position LFM kernel
%
%	Description:
%
%	K = LFMJXLFMKERNCOMPUTE(LFMKERN1, LFMKERN2, T) computes cross kernel
%	terms between jolt and position LFM kernels.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  LFMKERN1 - the kernel structure associated with the first LFM
%	   kernel.
%	  LFMKERN2 - the kernel structure associated with the second LFM
%	   kernel.
%	  T - inputs for which kernel is to be computed.
%
%	K = LFMJXLFMKERNCOMPUTE(LFMKERN1, LFMKERN2, T1, T2) computes cross
%	kernel terms between jolt and position LFM kernels.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  LFMKERN1 - the kernel structure associated with the first LFM
%	   kernel.
%	  LFMKERN2 - the kernel structure associated with the second LFM
%	   kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.


%	Copyright (c) 2010 Mauricio Alvarez
% 	lfmjXlfmKernCompute.m SVN version 1519
% 	last update 2011-07-22T13:08:16.000000Z

if nargin < 4
    t2 = t1;
end
if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end

if lfmKern1.inverseWidth ~= lfmKern2.inverseWidth
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

% Get length scale out.
sigma2 = 2/lfmKern1.inverseWidth;
sigma = sqrt(sigma2);

% Parameters of the kernel
alpha(1) = lfmKern1.damper./(2*lfmKern1.mass);
alpha(2) = lfmKern2.damper./(2*lfmKern2.mass);
omega(1) = sqrt(lfmKern1.spring./lfmKern1.mass - alpha(1)*alpha(1));
omega(2) = sqrt(lfmKern2.spring./lfmKern2.mass - alpha(2)*alpha(2));

% Precomputations to increase the speed
preExp1 = zeros(length(t1),2);
preExp2 = zeros(length(t2),2);
gamma1_p = alpha(1) + j*omega(1);
gamma1_m = alpha(1) - j*omega(1);
gamma2_p = alpha(2) + j*omega(2);
gamma2_m = alpha(2) - j*omega(2);
preGamma(1) = gamma1_p + gamma2_p;
preGamma(2) = gamma1_p + gamma2_m;
preGamma(3) = gamma1_m + gamma2_p;
preGamma(4) = gamma1_m + gamma2_m;
preConst = 1./preGamma;
preFactors(1) = preConst(2) - preConst(1);
preFactors(2) = preConst(3) - preConst(4);
preFactors(3) = preConst(3) - preConst(1);
preFactors(4) = preConst(2) - preConst(4);
preExp1(:,1) = (gamma1_p^3)*exp(-gamma1_p*t1);
preExp1(:,2) = (gamma1_m^3)*exp(-gamma1_m*t1);
preExp2(:,1) = exp(-gamma2_p*t2);
preExp2(:,2) = exp(-gamma2_m*t2);
% Actual computation of the kernel
sK = lfmComputeH3JP(gamma1_p, gamma1_m, sigma2, t1,t2,preFactors([1 2]), 0) + ...
    lfmComputeH3JP(gamma2_p, gamma2_m, sigma2, t2,t1,preFactors([3 4]), 1).' + ...
    lfmComputeH4JP(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExp2, 0 ) + ...
    lfmComputeH4JP(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExp1, 1 ).';

if lfmKern1.isNormalised
    K0 =  lfmKern1.sensitivity*lfmKern2.sensitivity/(8*sqrt(2)*lfmKern1.mass*lfmKern2.mass*prod(omega));
else
    K0 =  sigma*sqrt(pi)*lfmKern1.sensitivity*lfmKern2.sensitivity/(8*lfmKern1.mass*lfmKern2.mass*prod(omega));    
end

K = K0*sK;




