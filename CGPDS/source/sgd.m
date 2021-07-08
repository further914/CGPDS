function [x, options, flog, pointlog, scalelog] = sgd(f, x, options, gradf, varargin)
%SGD Stochastic Gradient Desent

% global experimentNo
%  Set up the options.
if length(options) < 18
  error('Options vector too short')
end

if(options(14))
  niters = options(14);
else
  niters = 100;
end

display = options(1);
gradcheck = options(9);

% Set up strings for evaluating function and gradient
f = fcnchk(f, length(varargin));
gradf = fcnchk(gradf, length(varargin));

%  Check gradients
if (gradcheck)
  feval('gradchek', x, f, gradf, varargin{:});
end

sigma0 = 1.0e-4;
% global OuterIter;

% ELBOName = strcat('TrainELBO',num2str(OuterIter));

% eval([ELBOName '= zeros(' num2str(niters) '+1,1);']);

% global ELBO;
% global lastIter;

j = 1;

stepSize = varargin{3};
while j < niters
    fnow = feval(f, x, varargin{:});
    grad = feval(gradf, x, varargin{:});
    
    x = x - stepSize*grad;
    fprintf(1, 'Cycle %4d  Error %11.6f  Scale %e  Gradient %f\n', j, fnow, stepSize,grad*grad');
    
    j = j + 1;
end




