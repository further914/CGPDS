function kern = indexKernParamInit(kern)

% INDEXKERNPARAMINIT INDEX kernel parameter initialisation.
%
%	Description:
%	The index covariance function returns a similarity if two indices from the input match. In other words if x_i is one index and x_j is another index then
%	
%	k(x_i, x_j) = sigma2 * \delta_{x_i, x_j}
%	
%	where \delta_{k, l} is the Kronecker delta which is one if k=l and
%	zero otherwise.
%	
%
%	KERN = INDEXKERNPARAMINIT(KERN) initialises the index based
%	covariance function kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2011 Neil D. Lawrence
% 	indexKernParamInit.m SVN version 1566
% 	last update 2011-08-07T06:05:18.519297Z

% These parameters are restricted to lie between 0 and 1.
  kern.variance = 1;
  kern.nParams = 1;
  
  kern.transforms(1).index = [1];
  kern.transforms(1).type = optimiDefaultConstraint('positive');

  kern.isStationary = true;
end