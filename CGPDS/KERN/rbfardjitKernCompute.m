function [k, k2] = rbfardjitKernCompute(kern, x, x2)

% RBFARDJITKERNCOMPUTE Compute the RBFARD kernel given the parameters and X.
%
%	Description:
%
%	RBFARDJITKERNCOMPUTE Implementation of the rbfard2 kernel that allow
%	for numerically robust computations (especially wrt to the variance
%	of the kernel, ie. sigmaf^2)
%	sigmaf^2 * [ exp( - 0.5*sum_k a_k (x_ik - x_jk)^2)  + jit*delta_ij ]
%	
%	The trick here is to have exp(..)  kernel and jitter white
%	kernel have a shared variance sigma_f^2
%	
%
%	K = RBFARDJITKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	automatic relevance determination radial basis function kernel given
%	a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	RBFARDJITKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, RBFARDJITKERNDIAGCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias
% 	rbfardjitKernCompute.m SVN version 1460
% 	last update 2011-07-04T18:55:17.584561Z


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
    
if nargin < 3
  n2 = dist2(x, x);
  k2 = kern.variance*exp(-n2*0.5);
  k = k2 +  (kern.variance*kern.jitter)*eye(size(x,1));
else
  x2 = x2*scales;
  n2 = dist2(x, x2);
  k = kern.variance*exp(-n2*0.5);
  k2 = k;
end
