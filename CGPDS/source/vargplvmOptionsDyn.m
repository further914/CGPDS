function optionsDyn = vargplvmOptionsDyn(optionsDyn)

% VARGPLVMOPTIONSDYN Fill in an options structure with default parameters
%
%	Description:
%
%	VARGPLVMOPTIONSDYN(OPTIONSDYN) takes an VARGPLVM options structure
%	amends it with default values wherever the field is empty.
%	 Arguments:
%	  OPTIONSDYN - the VARGPLVM structure to fill in with default
%	   values.
%	
%	
%
%	See also
%	VARGPLVMINITDYNAMICS, VARGPLVMADDDYNAMICS


%	Copyright (c) 2011 Michalis K. Titsias
%	Copyright (c) 2011 Neil D. Lawrence
%	Copyright (c) 2011 Andreas C. Damianou
% 	vargplvmOptionsDyn.m SVN version 1570
% 	last update 2011-08-30T14:57:48.000000Z


if isfield(optionsDyn, 'type') && ~isempty(optionsDyn.type)
    optionsDyn.type = 'vargpTime';
end


if strcmp(optionsDyn.type, 'vargpTime')
    if ~isfield(optionsDyn, 'initX')
        optionsDyn.initX = 'ppca';
    end
    
    % Create time vector for each dimension; if t is not given, assume an
    % equally spaced time vector.
    if  ~isfield(optionsDyn, 't') || isempty(optionsDyn.t)
        fprintf(1, '# Time vector unknown; creating random, equally spaced time vector\n');
        t = linspace(0, 2*pi, size(model.X, 1)+1)';
        t = t(1:end-1, 1);
        optionsDyn.t = t;
    else
        t=optionsDyn.t;
    end
    
    
    % A better way to initialize the  kernel hyperparameter,
    % especially lengthscale, should be to learn it by running few iterations of
    % GP regression marginal likelihood maximization given as data the PCA output
    % (the default value is jsut reasonable and it sets the inverse lenthscale to quite
    % small value so the GP dynamic prior is weaker (less smoothed) at
    % initializations
    if  ~isfield(optionsDyn, 'kern') || isempty(optionsDyn.kern)
       % kern = kernCreate(t, {'rbf','white'});
        kern = kernCreate(t, {'rbf','white', 'bias'});
        kern.comp{2}.variance = 1e-3; % 1e-1
        
        if  isfield(optionsDyn, 'inverseWidth')
            invWidth=optionsDyn.inverseWidth;
        else
            invWidth=5;
        end
        
        % The following is related to the expected number of zero-crossings.
        kern.comp{1}.inverseWidth = invWidth./(((max(t)-min(t))).^2);
        kern.comp{1}.variance = 1;
        
        optionsDyn.kern = kern;
    end
    
    
    if ~isfield(optionsDyn,'seq')
        optionsDyn.seq=[];
    end
    
    % Set to 1 to reoptimise all inducing points during test time
     if ~isfield(optionsDyn, 'testReoptimise')
     	optionsDyn.testReoptimise = 1;  % !!! SET THE DEFAULT TO 1
     end
    
     % Set to 0 so that if you use a kernel from the exp. family for the
     % dynamics, then its variance is not learned and the lenghtscale
     % compensates this (we have one less free parameter then).
    if ~isfield(optionsDyn, 'learnVariance')
        optionsDyn.learnVariance = 0; % !!! SET THE DEFAULT TO 0
    end
    
    % If set to 1 (not recommended) then the initial means have some
    % smoothing in the first and last values (i.e. closer to zero). Useful
    % mostly for periodic kernels.
    if ~isfield(optionsDyn, 'regularizeMeans')
        optionsDyn.regularizeMeans = 0; 
     end
    
end

