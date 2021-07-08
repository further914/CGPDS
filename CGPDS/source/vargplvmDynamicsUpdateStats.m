function model = vargplvmDynamicsUpdateStats(model)

% VARGPLVMDYNAMICSUPDATESTATS wrapper function which according to the type
%
%	Description:
%	of the model dynamics calls the appropriate function to perform the
%	precomputations needed for the dynamics model
%	DESC
%	ARG model: the model that contains the dynamics for which the statistics
%	must be updated
%	RETURN model : the updated model
%	
%	


%	Copyright (c) Michalis Titsias, Neil Lawrence, 2010-2011 % COPYRIGHT Andreas C. Damianou
% 	vargplvmDynamicsUpdateStats.m SVN version 1412
% 	last update 2011-06-02T19:53:06.000000Z
fhandle = str2func([model.dynamics.type 'UpdateStats']);
model = fhandle(model);
