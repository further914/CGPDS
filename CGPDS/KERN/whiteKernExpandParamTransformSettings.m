function kern = whiteKernExpandParamTransformSettings(kern, paramtransformsettings)

% WHITEKERNEXPANDPARAMTRANSFORMSETTINGS Create kernel structure from WHITE kernel's
%
%	Description:
%	parameter transformation settings.
%	
%
%	WHITEKERNEXPANDPARAMTRANSFORMSETTINGS
%	DESC returns a white noise kernel structure filled with the
%	parameter transformation settings in the given cell array.
%	This is used as a helper function to enable parameters to
%	be optimised in, for example, the NETLAB optimisation functions.
%	
%	ARG kern : the kernel structure in which the parameter
%	transformation settings are to be placed.
%	
%	ARG paramtransformsettings : cell array of parameter
%	transformation settings which are to be placed in the
%	kernel structure. Each setting needs to correspond to
%	the transformations used by the kernel. For example, some
%	transformation might need knowledge of the desired output range,
%	so the transformation setting could be e.g. [0 100], whereas some
%	other transformations might need more complicated setting
%	information.
%	
%	RETURN kern : kernel structure with the given parameter
%	transformation settings in the relevant locations.
%	
%	
%	
%
%	See also
%	WHITEKERNPARAMINIT, WHITEKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2011 Jaakko Peltonen
% 	whiteKernExpandParamTransformSettings.m SVN version 1567
% 	last update 2011-08-09T18:31:29.899803Z


% The "white" kernel uses just one transformation, for the variance 
% parameter of the kernel.
if length(paramtransformsettings~=1),
  error(sprintf('Problem in whiteKernExpandParamTransformSettings: expected 1 transformation setting, received %d\n',...
        length(paramtransformsettings)));
end;

kern.transforms(1).transformsettings=paramtransformsettings{1};
