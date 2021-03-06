function kern = kernExpandParam(kern, params)

% KERNEXPANDPARAM Expand parameters to form a kernel structure.
%
%	Description:
%
%	KERN = KERNEXPANDPARAM(KERN, PARAM) returns a kernel structure
%	filled with the parameters in the given vector. This is used as a
%	helper function to enable parameters to be optimised in, for
%	example, the NETLAB optimisation functions.
%	 Returns:
%	  KERN - kernel structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  KERN - the kernel structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the kernel
%	   structure.
%	As well as extracting the parameters, some transformation of
%	parameters is also undertaken in this file. If the field
%	transforms is not empty, it dictactes how the kernel parameters
%	are to be transformed (for example by a exponential to keep them
%	positive). The field transforms may also include custom settings
%	for each transformation, for example a custom output range.
%	
%	
%	
%
%	See also
%	KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2003, 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2011 Jaakko Peltonen
% 	kernExpandParam.m CVS version 1.5
% 	kernExpandParam.m SVN version 1537
% 	last update 2011-08-03T14:11:43.315503Z


% Check if any parameters are being optimised in a transformed space,
% and if they are, apply the transformations to the parameters
% before giving them to the kernel.
if ~isempty(kern.transforms)
  % Process each transformation. Each transformation may affect
  % several parameters.
  for i = 1:length(kern.transforms)
    % Get the parameter indices affected by the i:th transformation.
    index = kern.transforms(i).index;
    
    % Get the function handle of the transformation
    fhandle = str2func([kern.transforms(i).type 'Transform']);
    
    % If custom settings have been provided for the transformation,
    % use them, otherwise call the transformation with no settings
    if isfield(kern.transforms(i),'transformsettings'),    
      params(index) = fhandle(params(index), 'atox', kern.transforms(i).transformsettings);    
    else
      params(index) = fhandle(params(index), 'atox');    
    end
    
    % fprintf(1, 'applying transformation, index %d: %f, A %f, B %f, result %f\n',index,origparams(index),kern.transforms(i).transformsettings(1),kern.transforms(i).transformsettings(2),params(index));
  end
end


% Now that any transformation have been applied, give the
% transformed parameters to the kernel using the function
% specific to the current kernel type.
fhandle = str2func([kern.type 'KernExpandParam']);
kern = fhandle(kern, params);
