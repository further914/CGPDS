function g = whiteKernGradient(kern, x, varargin)

% WHITEKERNGRADIENT Gradient of white-noise (WHITE) kernel's
%
%	Description:
%	parameters.
%	
%
%	WHITEKERNGRADIENT
%	DESC computes the gradient of functions with respect to the
%	white noise kernel's parameters. As well as the kernel structure
%	and the input positions, the user provides a matrix PARTIAL which
%	gives the partial derivatives of the function with respect to the
%	relevant elements of the kernel matrix.
%	
%	ARG kern : the kernel structure for which the gradients are being
%	computed.
%	
%	ARG x : the input locations for which the gradients are being
%	computed.
%	
%	ARG partial (only element of varargin): matrix of partial
%	derivatives of the function of interest with respect to the
%	kernel matrix. The argument takes the form of a square matrix of
%	dimension  numData, where numData is the number of rows in X.
%	
%	RETURN g : gradients of the function of interest with respect to
%	the kernel parameters. The ordering of the vector should match
%	that provided by the function kernExtractParam.
%	
%	ALTERNATIVE FORMAT:
%	
%	DESC computes the derivatives as above, but input locations are
%	now provided in two matrices associated with rows and columns of
%	the kernel matrix.
%	
%	ARG kern : the kernel structure for which the gradients are being
%	computed.
%	
%	ARG x1 : the input locations associated with the rows of the
%	kernel matrix.
%	
%	ARG x2 (first element of varargin): the input locations
%	associated with the columns of the kernel matrix.
%	
%	ARG partial (second element of varargin): matrix of partial
%	derivatives of the function of interest with respect to the
%	kernel matrix. The matrix should have the same number of rows as
%	X1 and the same number of columns as X2 has rows.
%	
%	RETURN g : gradients of the function of interest with respect to
%	the kernel parameters.
%	
%	
%	
%
%	See also
%	% SEEALSO WHITEKERNPARAMINIT, KERNGRADIENT, WHITEKERNDIAGGRADIENT, KERNGRADX


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2011 Jaakko Peltonen
% 	whiteKernGradient.m CVS version 1.5
% 	whiteKernGradient.m SVN version 1566
% 	last update 2011-08-03T14:11:43.275502Z


  if nargin < 4
    g(1, 1) = trace(varargin{end});
  else
    g = 0;
  end
end