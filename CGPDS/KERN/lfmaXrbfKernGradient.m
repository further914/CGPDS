function [g1, g2] = lfmaXrbfKernGradient(lfmKern, rbfKern, t1, t2, covGrad, meanVector)

% LFMAXRBFKERNGRADIENT Compute gradient between the LFMA and RBF kernels.
%
%	Description:
%
%	[G1, G2] = LFMAXRBFKERNGRADIENT(LFMKERN, RBFKERN, T) computes the
%	gradient of an objective function with respect to cross kernel terms
%	between LFMA and RBF kernels for the multiple output kernel.
%	 Returns:
%	  G1 - gradient of objective function with respect to kernel
%	   parameters of LFMA kernel.
%	  G2 - gradient of objective function with respect to kernel
%	   parameters of RBF kernel.
%	 Arguments:
%	  LFMKERN - the kernel structure associated with the LFMA kernel
%	   (acceleration).
%	  RBFKERN - the kernel structure associated with the RBF kernel.
%	  T - inputs for which kernel is to be computed.
%
%	[G1, G2] = LFMAXRBFKERNGRADIENT(LFMKERN, RBFKERN, T1, T2) computes
%	the gradient of an objective function with respect to cross kernel
%	terms between LFMA and RBF kernels for the multiple output kernel.
%	 Returns:
%	  G1 - gradient of objective function with respect to kernel
%	   parameters of LFMA kernel.
%	  G2 - gradient of objective function with respect to kernel
%	   parameters of RBF kernel.
%	 Arguments:
%	  LFMKERN - the kernel structure associated with the LFMV kernel
%	   (acceleration).
%	  RBFKERN - the kernel structure associated with the RBF kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%
%	[G1, G2] = LFMAXRBFKERNGRADIENT(LFMKERN, RBFKERN, T1, T2, MEANVEC)
%	computes the gradient of an objective function with respect to cross
%	kernel terms between LFMA and RBF kernels for the multiple output
%	kernel.
%	 Returns:
%	  G1 - gradient of objective function with respect to kernel
%	   parameters of LFMA kernel.
%	  G2 - gradient of objective function with respect to kernel
%	   parameters of RBF kernel.
%	 Arguments:
%	  LFMKERN - the kernel structure associated with the LFMV kernel
%	   (acceleration).
%	  RBFKERN - the kernel structure associated with the RBF kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	  MEANVEC - precomputed factor that is used for the switching
%	   dynamical latent force model.


%	Copyright (c) 2010 Mauricio A. Alvarez
% 	lfmaXrbfKernGradient.m SVN version 809
% 	last update 2010-05-28T06:01:33.000000Z

subComponent = false; % This is just a flag that indicates if this kernel is part of a bigger kernel (SDLFM)

if nargin == 4
    covGrad = t2;
    t2 = t1;
elseif nargin == 6
    subComponent = true;
    if numel(meanVector)>1
        if size(meanVector,1) == 1,
            if size(meanVector, 2)~=size(covGrad, 2)
                error('The dimensions of meanVector don''t correspond to the dimensions of covGrad')
            end
        else
            if size((meanVector'), 2)~=size(covGrad,2)
                error('The dimensions of meanVector don''t correspond to the dimensions of covGrad')
            end
        end
    else
        if numel(t1)==1 && numel(t2)>1
            % matGrad will be row vector and so should be covGrad
            dimcovGrad = length(covGrad);
            covGrad = reshape(covGrad, [1 dimcovGrad]);
        elseif numel(t1)>1 && numel(t2)==1
            % matGrad will be column vector and sp should be covGrad
            dimcovGrad = length(covGrad);
            covGrad = reshape(covGrad, [dimcovGrad 1]);
        end
    end
end

if size(t1, 2) > 1 || size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern.inverseWidth ~= rbfKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end

m = lfmKern.mass;
D = lfmKern.spring;
C = lfmKern.damper;
S = lfmKern.sensitivity;

alpha = C/(2*m);
omega = sqrt(D/m-alpha^2);

sigma2 = 2/lfmKern.inverseWidth;% Tamporarly changed by MA
sigma = sqrt(sigma2);

gamma1 = alpha + j*omega;
gamma2 = alpha - j*omega;
[ComputeUpsilon1, ComputeUpsilon1Local] = lfmapComputeUpsilonMatrix(gamma2,sigma2,t1, t2, 0);
[ComputeUpsilon2, ComputeUpsilon2Local] = lfmapComputeUpsilonMatrix(gamma1,sigma2,t1, t2, 0);

GradientUpsilon1 = lfmapGradientUpsilonMatrix(gamma2,sigma2, t1, t2, 0, ComputeUpsilon1Local);
GradientUpsilon2 = lfmapGradientUpsilonMatrix(gamma1,sigma2, t1, t2, 0, ComputeUpsilon2Local);

if lfmKern.isNormalised
   K0 = lfmKern.sensitivity/(j*4*sqrt(2)*lfmKern.mass*omega);
   K02 = 1/(j*4*sqrt(2)*lfmKern.mass*omega);
else
   K0 = sqrt(pi)*sigma*lfmKern.sensitivity/(j*4*lfmKern.mass*omega);
   K02 = sqrt(pi)*sigma/(j*4*lfmKern.mass*omega);
end

g1 = zeros(1,5);

% Gradient with respect to m, C and D

for ind = 1:3 % Parameter (m, D or C)
    switch ind
        case 1  % Gradient wrt m
            gradThetaM = 1;
            gradThetaAlpha = -C/(2*(m^2));
            gradThetaOmega = (C^2-2*m*D)/(2*(m^2)*sqrt(4*m*D-C^2));
        case 2  % Gradient wrt D
            gradThetaM = 0;
            gradThetaAlpha = 0;
            gradThetaOmega = 1/sqrt(4*m*D-C^2);
        case 3  % Gradient wrt C
            gradThetaM = 0;
            gradThetaAlpha = 1/(2*m);
            gradThetaOmega = -C/(2*m*sqrt(4*m*D-C^2));
    end

    % Gradient evaluation
    gradThetaGamma1 = gradThetaAlpha + j*gradThetaOmega;
    gradThetaGamma2 = gradThetaAlpha - j*gradThetaOmega;

    matGrad = K0*(GradientUpsilon1*gradThetaGamma2 -  GradientUpsilon2*gradThetaGamma1 ...
        - (gradThetaM/lfmKern.mass + gradThetaOmega/omega) ...
        * (ComputeUpsilon1 - ComputeUpsilon2));

    if subComponent
        if size(meanVector,1) ==1,
            matGrad = matGrad*meanVector;
        else
            matGrad = (meanVector*matGrad).';
        end
    end
    
    g1(ind) = sum(sum(matGrad.*covGrad));
end

% Gradient with respect to sigma
gamma1 = alpha + j*omega;
gamma2 = alpha - j*omega;

if lfmKern.isNormalised
    matGrad = K0*(lfmapGradientSigmaUpsilonMatrix(gamma2,sigma2,t1,t2,0) ...
        - lfmapGradientSigmaUpsilonMatrix(gamma1,sigma2,t1,t2,0));
else
    matGrad = (sqrt(pi)*S/(j*4*m*omega)) ...
        *(ComputeUpsilon1 - ComputeUpsilon2 ...
        + sigma*(lfmapGradientSigmaUpsilonMatrix(gamma2,sigma2,t1,t2,0) ...
        - lfmapGradientSigmaUpsilonMatrix(gamma1,sigma2,t1,t2,0)));
end

if subComponent
    if size(meanVector,1) ==1,
        matGrad = matGrad*meanVector;
    else
        matGrad = (meanVector*matGrad).';
    end
end

g1(4) = sum(sum(matGrad.*covGrad))*(-(sigma^3)/4); % temporarly introduced by MA
g2(1) = g1(4);

% Gradient with respect to S

matGrad = K02*(ComputeUpsilon1 - ComputeUpsilon2);

if subComponent
    if size(meanVector,1) ==1,
        matGrad = matGrad*meanVector;
    else
        matGrad = (meanVector*matGrad).';
    end
end

g1(5) = sum(sum(matGrad.*covGrad));
g1 = real(g1);

% Gradient with respect to the "variance" of the RBF
g2(1) = 0; % Otherwise is counted twice, temporarly changed by MA
g2(2) = 0;
