function f = vargplvmObjective(params, model ,stepSize, samd, samr)

% VARGPLVMOBJECTIVE Wrapper function for variational GP-LVM objective.
%
%	Description:
%
%	F = VARGPLVMOBJECTIVE(PARAMS, MODEL) provides a wrapper function for
%	the variational GP-LVM, it takes the negative of the log likelihood,
%	feeding the parameters correctly to the model.
%	 Returns:
%	  F - the negative of the log likelihood of the model.
%	 Arguments:
%	  PARAMS - the parameters of the variational GP-LVM model.
%	  MODEL - the model structure in which the parameters are to be
%	   placed.
%	
%
%	See also
%	VARGPLVMCREATE, VARGPLVMLOGLIKELIHOOD, VARGPLVMEXPANDPARAM


%	Copyright (c) 2009-2011 Michalis K. Titsias
%	Copyright (c) 2009-2011 Neil D. Lawrence
% 	vargplvmObjective.m SVN version 1582
% 	last update 2011-09-07T17:23:12.321178Z


model = modelExpandParam(model, params,samd,samr);
fhandle = str2func([model.type 'LogLikelihood']);
f = - fhandle(model ,stepSize, samd, samr);
