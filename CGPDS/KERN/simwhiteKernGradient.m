function g = simwhiteKernGradient(kern, t1, varargin)

% SIMWHITEKERNGRADIENT Gradient of SIM-WHITE kernel's parameters.
%
%	Description:
%
%	G = SIMWHITEKERNGRADIENT(KERN, T1, PARTIAL) computes the gradient of
%	functions with respect to the SIM-White (Single Input Motif - White)
%	kernel's parameters. As well as the kernel structure and the input
%	positions, the user provides a matrix PARTIAL which gives the
%	partial derivatives of the function with respect to the relevant
%	elements of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  T1 - the input locations for which the gradients are being
%	   computed.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in t.
%
%	G = SIMWHITEKERNGRADIENT(KERN, T1, T2, PARTIAL) computes the
%	derivatives as above, but input locations are now provided in two
%	matrices associated with rows and columns of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  T1 - the input locations associated with the rows of the kernel
%	   matrix.
%	  T2 - the input locations associated with the columns of the kernel
%	   matrix.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as t1 and the same number of columns as t2
%	   has rows.
%	
%
%	See also
%	% SEEALSO SIMWHITEKERNPARAMINIT, KERNGRADIENT, SIMWHITEKERNDIAGGRADIENT, KERNGRADX


%	Copyright (c) 2009 David Luengo
% 	simwhiteKernGradient.m SVN version 198
% 	last update 2011-06-16T07:23:44.000000Z


if length(varargin)<2
  t2 = t1;
else
  t2 = varargin{1};
end

[g1, g2] = simwhiteXsimwhiteKernGradient(kern, kern, t1, t2, varargin{end});

g = real(g1 + g2);
