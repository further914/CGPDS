function f = cgpdsObjective(params, model ,stepSize, samd, samr)

% cgpdsObjective Wrapper function for CGPDS objective.
%
%	Description:
%
%	F = cgpdsObjective(PARAMS, MODEL) provides a wrapper function for
%	the CGPDS, it takes the negative of the log likelihood,
%	feeding the parameters correctly to the model.
%	 Returns:
%	  F - the negative of the log likelihood of the model.
%	 Arguments:
%	  PARAMS - the parameters of the CGPDS.
%	  MODEL - the model structure in which the parameters are to be
%	   placed.
%	
%	Copyright (c) 2009-2011 Michalis K. Titsias
%	Copyright (c) 2009-2011 Neil D. Lawrence


model = modelExpandParam(model, params,samd,samr);
fhandle = str2func([model.type 'LogLikelihood']);

f = - fhandle(model ,stepSize, samd, samr);
