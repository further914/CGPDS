function kern = indexardKernParamInit(kern, indices)

% INDEXARDKERNPARAMINIT INDEXARD kernel parameter initialisation.
%
%	Description:
%	The index covariance function returns a similarity if two indices from the input match. In other words if x_i is one index and x_j is another index then
%	
%	k(x_i, x_j) = sigma2_{x_i} * \delta_{x_i, x_j}
%	
%	where \delta_{k, l} is the Kronecker delta which is one if k=l and
%	zero otherwise. In the ARD covaraince the value of the variance
%	associated with each index can also different.
%	
%
%	KERN = INDEXARDKERNPARAMINIT(KERN) initialises the index based
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
% 	indexardKernParamInit.m SVN version 1566
% 	last update 2011-08-07T06:30:20.608661Z
  
  
  kern.indexValues = true;
  if ~isfield(kern, 'indices')
    kern.indices = [0 1 2 3];
  end
  
  kern.indexScales = ones(1, length(kern.indices));
  kern.nParams = length(kern.indices);
  
  kern.transforms(1).index = [1:kern.nParams];
  kern.transforms(1).type = optimiDefaultConstraint('positive');

  kern.isStationary = false;
end