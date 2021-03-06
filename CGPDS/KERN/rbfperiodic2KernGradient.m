function g = rbfperiodic2KernGradient(kern, x, varargin)

% RBFPERIODIC2KERNGRADIENT Gradient of RBFPERIODIC2 kernel's parameters.
%
%	Description:
%
%	G = RBFPERIODIC2KERNGRADIENT(KERN, X, PARTIAL) computes the gradient
%	of functions with respect to the RBF periodic covariance with
%	variying period kernel's parameters. As well as the kernel structure
%	and the input positions, the user provides a matrix PARTIAL which
%	gives the partial derivatives of the function with respect to the
%	relevant elements of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  X - the input locations for which the gradients are being
%	   computed.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in X.
%
%	G = RBFPERIODIC2KERNGRADIENT(KERN, X1, X2, PARTIAL) computes the
%	derivatives as above, but input locations are now provided in two
%	matrices associated with rows and columns of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  X1 - the input locations associated with the rows of the kernel
%	   matrix.
%	  X2 - the input locations associated with the columns of the kernel
%	   matrix.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as X1 and the same number of columns as X2
%	   has rows.
%	
%	
%	
%
%	See also
%	% SEEALSO RBFPERIODIC2KERNPARAMINIT, KERNGRADIENT, RBFPERIODIC2KERNDIAGGRADIENT, KERNGRADX


%	Copyright (c) 2007, 2009 Neil D. Lawrence


%	With modifications by Andreas C. Damianou 2011


%	With modifications by Michalis K. Titsias 2011
% 	rbfperiodic2KernGradient.m SVN version 1519
% 	last update 2011-07-22T12:55:52.000000Z

% The last argument is covGrad

if nargin < 4
  [k, sk, expArg] = rbfperiodic2KernCompute(kern, x);
else
  [k, sk, expArg] = rbfperiodic2KernCompute(kern, x, varargin{1});
end
g(1) = - 2*sum(sum(varargin{end}.*k.*expArg)); %d(kern)/d(inv.width)
g(2) =  sum(sum(varargin{end}.*sk)); %d(kern)/d(variance)

% d(kern)/d(factor):
% Note: The folowing can also be written easier using the trigonometric
% formulas for the half of an angle, and get rid of the cosine
factor = kern.factor;
if nargin<4
	xx=repmat(x, 1, size(x, 1)) - repmat(x', size(x, 1), 1);
else
    x2 = varargin{1};
    xx = repmat(x, 1, size(x2, 1)) - repmat(x2', size(x, 1), 1);
end

arg = (factor*xx) / 2;
sinarg = sin(arg);
n2 = sinarg.*sinarg;
wi2 = (2 .* kern.inverseWidth);
rbfPart = kern.variance*exp(-n2*wi2);


g(3)= -2*sum(sum(varargin{end}.*xx.*kern.inverseWidth.*cos(arg).*sinarg.*rbfPart));

% g(3) = -xx*kern.inverseWidth*cos(arg).*sinarg.*2.*rbfPart;


%gX = zeros(size(X2));
%arg = factor*(X2 - x)/2;
%sinarg = sin(arg);
%n2 = sinarg.*sinarg;
%wi2 = (2 .* kern.inverseWidth);
%rbfPart = kern.variance*exp(-n2*wi2);

%gX = factor*kern.inverseWidth*cos(arg).*sinarg.*2.*rbfPart;


