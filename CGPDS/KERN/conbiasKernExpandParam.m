function kern = conbiasKernExpandParam(kern, params)

% BIASKERNEXPANDPARAM Create kernel structure from BIAS kernel's parameters.
%
%	Description:
%
%	KERN = BIASKERNEXPANDPARAM(KERN, PARAM) returns a bias kernel
%	structure filled with the parameters in the given vector. This is
%	used as a helper function to enable parameters to be optimised in,
%	for example, the NETLAB optimisation functions.
%	 Returns:
%	  KERN - kernel structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  KERN - the kernel structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the kernel
%	   structure.
%	
%
%	See also
%	BIASKERNPARAMINIT, BIASKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
% 	biasKernExpandParam.m CVS version 1.3
% 	biasKernExpandParam.m SVN version 1
% 	last update 2011-06-16T07:23:44.000000Z


kern.variance = params(1);