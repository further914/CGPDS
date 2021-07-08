function [g1,g2, covGradLocal] = sdlfmaXsdlfmvKernGradient(sdlfmaKern1, sdlfmvKern2, ...
    t1, t2, covGrad, covIC)

% SDLFMAXSDLFMVKERNGRADIENT Gradients of cross kernel between 2 SDLFM kernels.
%
%	Description:
%
%	[G1, G2, COVGRADLOCAL] = SDLFMAXSDLFMVKERNGRADIENT(SDLFMVKERN1,
%	SDLFMKERN2, T1, COVGRAD, COVIC) computes a cross gradient for a
%	cross kernel between two switching dynamical LFM kernels for the
%	multiple output kernel. The first SDLFM corresponds to the accel.
%	and the second to a velocity.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see lfmKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see lfmKernExtractParam.
%	  COVGRADLOCAL - partial derivative wrt the initial conditions of
%	   the firts interval.
%	 Arguments:
%	  SDLFMVKERN1 - the kernel structure associated with the first SDLFM
%	   kernel (acceleration).
%	  SDLFMKERN2 - the kernel structure associated with the second SDLFM
%	   kernel (velocity).
%	  T1 - inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%	  COVIC - covariance for the initial conditions
%
%	[G1, G2, COVGRADLOCAL] = SDLFMAXSDLFMVKERNGRADIENT(SDLFMVKERN1,
%	SDLFMKERN2, T1, T2, COVGRAD, COVIC) computes a cross gradient for a
%	cross kernel between two switching dynamical LFM kernels for the
%	multiple output kernel. The first SDLFM corresponds to the accel.
%	and the second to a velocity.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see lfmKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see lfmKernExtractParam.
%	  COVGRADLOCAL - partial derivative wrt the initial conditions of
%	   the first interval.
%	 Arguments:
%	  SDLFMVKERN1 - the kernel structure associated with the first SDLFM
%	   kernel (acceleration).
%	  SDLFMKERN2 - the kernel structure associated with the second SDLFM
%	   kernel (velocity).
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%	  COVIC - covariance for the initial conditions


%	Copyright (c) 2010 Mauricio A. Alvarez
% 	sdlfmaXsdlfmvKernGradient.m SVN version 807
% 	last update 2010-05-28T06:01:33.000000Z

if nargin == 4
    covGrad = t2;
    t2 = t1;
end

[g1, g2, covGradLocal] = sdlfmXsdlfmKernGradient(sdlfmaKern1, sdlfmvKern2, t1, ...
    t2, covGrad, covIC, 'AccelVel');