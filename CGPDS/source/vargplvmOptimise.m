function model = vargplvmOptimise(model, display, iters, stepSize, samd, samr, varargin)

% VARGPLVMOPTIMISE Optimise the VARGPLVM.
%
%	Description:
%
%	MODEL = VARGPLVMOPTIMISE(MODEL, DISPLAY, ITERS) takes a given GP-LVM
%	model structure and optimises with respect to parameters and latent
%	positions.
%	 Returns:
%	  MODEL - the optimised model.
%	 Arguments:
%	  MODEL - the model to be optimised.
%	  DISPLAY - flag dictating whether or not to display optimisation
%	   progress (set to greater than zero) (default value 1).
%	  ITERS - number of iterations to run the optimiser for (default
%	   value 2000).
%	vargplvmLogLikeGradients, vargplvmObjective, vargplvmGradient
%	
%	
%
%	See also
%	VARGPLVMCREATE, VARGPLVMLOGLIKELIHOOD, 


%	Copyright (c) 2009 Michalis K. Titsias
%	Copyright (c) 2005, 2006 Neil D. Lawrence
% 	vargplvmOptimise.m SVN version 1720
% 	last update 2011-11-15T01:51:25.827041Z


if nargin < 3
  iters = 2000;
  if nargin < 2
    display = 1;
  end
end

options = optOptions;
if length(varargin) == 2
    if strcmp(varargin{1}, 'gradcheck')
        assert(islogical(varargin{2}));
        options(9) = varargin{2};
        if options(9)
            [params, names] = modelExtractParam(model, samd, samr);
            for i=1:length(names)
                fprintf('%d\t%s\n', i, names{i});
            end
            feval('gradchek2', params, @vargplvmObjective, @vargplvmGradient, model, names);
        else
            params = modelExtractParam(model, samd, samr);
        end
    end
else
    params = modelExtractParam(model, samd, samr);
end

options(2) = 0.1*options(2); 
options(3) = 0.1*options(3);

if display
  options(1) = 1;
  if length(params) <= 100
    options(9) = 1;
  end
end
options(14) = iters;

if isfield(model, 'optimiser') && ~isa(model.optimiser, 'function_handle')
    
    if isfield(model, 'optimiser')
      optim = str2func(model.optimiser);
    else
      optim = str2func('scg');
    end
    

    if strcmp(func2str(optim), 'optimiMinimize')
        % Carl Rasmussen's minimize function 
        params = optim('vargplvmObjectiveGradient', params, options, model);
    elseif strcmp(func2str(optim), 'scg2')
        % NETLAB style optimization with a slight modification so that an
        % objectiveGradient can be used where applicable, in order to re-use
        % precomputed quantities.
        params = optim('vargplvmObjectiveGradient', params,  options,  'vargplvmGradient', model);
    elseif strcmp(func2str(optim), 'scg')
        % NETLAB style optimization.
        params = optim('vargplvmObjective', params,  options,  'vargplvmGradient', model ,stepSize, samd, samr);
    elseif strcmp(func2str(optim), 'Adam')
        params = optim('vargplvmObjective', params,options, 'vargplvmGradient', model ,stepSize, samd, samr);
    end
    
elseif isfield(model, 'optimiser') && isa(model.optimiser, 'function_handle')
    f = fcnchk(model.optimiser);
    params = f(model);
else
    error('vargplvmOptimise: Invalid optimiser setting.');
end

%model = vargplvmExpandParam(model, params);
model = modelExpandParam(model, params,samd,samr);

