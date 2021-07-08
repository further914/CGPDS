function ll = modelVarPriorBound(model)

% MODELVARPRIORBOUND Wrapper function for the various types of the
%
%	Description:
%	variational GPLVM bound when dynamics is used.
%
%	MODELVARPRIORBOUND provides a wrapper function for the variational
%	GP-LVM, when there are dynamics. It takes a model structure and
%	according to its dynamics' type it calls the appropriate function to
%	calculate the bound that corresponds only to the dynamics part (a KL
%	term). The bound is a lower Jensens bound for variational
%	approximation.
%	ARG model : the model structure which contains the dynamics one, for
%	which the corresponding part of the bound is to be computed.
%	RETURN ll : the term of the variational bound which corresponds to the
%	dynamics structure found in the model.
%	
%	
%
%	See also
%	VARGPLVMLOGLIKELIHOOD, VARGPTIMEDYNAMICSVARPRIORBOUND


%	Copyright (c) 2010-2011 Michalis K. Titsias
%	Copyright (c) 2010-2011 Neil D. Lawrence
%	Copyright (c) 2010-2011 Andreas C. Damianou
% 	modelVarPriorBound.m SVN version 1312
% 	last update 2011-04-18T11:54:45.000000Z


fhandle = str2func([model.dynamics.type 'VarPriorBound']);
ll = fhandle(model.dynamics);