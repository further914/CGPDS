function model = cgpdsExpandParam(model, params, samd, samr)

% cgpdsExpandParam Expand a parameter vector into a CGPDS.
%
%	Description:
%
%	cgpdsExpandParam(MODEL, PARAMS) takes an CGPDS structure and a
%	vector of parameters, and fills the structure with the given
%	parameters. Also performs any necessary precomputation for
%	likelihood and gradient computations, so can be computationally
%	intensive to call.
%	 Arguments:
%	  MODEL - the VARGPLVM structure to put the parameters in.
%	  PARAMS - parameter vector containing the parameters to put in the
%	   VARGPLVM structure.
%	Copyright (c) 2009-2011 Michalis K. Titsias
%	Copyright (c) 2009-2011 Neil D. Lawrence
% 	vargplvmExpandParam.m SVN version 1572
% 	last update 2011-08-30T14:57:48.000000Z

startVal = 1;
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    % variational parameters (reparametrized) AND dyn.kernel's parameters
    endVal = model.dynamics.nParams;
    if isempty(find(isnan(params(startVal:endVal))))&&isempty(find(isinf(params(startVal:endVal))))
        model.dynamics = modelExpandParam(model.dynamics, params(startVal:endVal),samd,samr);
    end
%     if max(max(model.dynamics.vardist.means))>20
%         model.dynamics.vardist.means = model.dynamics.vardist.means./repmat(max(model.dynamics.vardist.means),size(model.dynamics.vardist.means,1),1);
%     end
%     if max(max(model.dynamics.vardist.covars))>1
%         model.dynamics.vardist.covars = model.dynamics.vardist.covars./repmat(max(model.dynamics.vardist.covars),size(model.dynamics.vardist.covars,1),1);
%     end
    
%     if ~isempty(find(isnan(model.dynamics.vardist.means)|isinf(model.dynamics.vardist.means)))
%         model.dynamics.vardist.means(find(isnan(model.dynamics.vardist.means)|isinf(model.dynamics.vardist.means)))...
%             = model.X(find(isnan(model.dynamics.vardist.means)|isinf(model.dynamics.vardist.means)));
%     end
%     
%     if ~isempty(find(isnan(model.dynamics.vardist.covars)|isinf(model.dynamics.vardist.covars)))
%         naninfLen = length(find(isinf(model.dynamics.vardist.covars)|isnan(model.dynamics.vardist.covars)));
%         model.dynamics.vardist.covars(find(isnan(model.dynamics.vardist.covars)|isinf(model.dynamics.vardist.covars)))...
%             = 0.5*ones(naninfLen,1) + 0.001*randn(naninfLen,1);
%     end
else
    % variational parameters (means and covariances), original ones
    endVal = model.vardist.nParams;
    model.vardist = modelExpandParam(model.vardist, params(startVal:endVal)); 
end

% inducing inputs
startVal = endVal+1;
% Parameters include inducing variables.
endVal = endVal + model.q*model.k;
if isempty(find(isnan(params(startVal:endVal))))&&isempty(find(isinf(params(startVal:endVal))))
    model.X_v = reshape(params(startVal:endVal),model.k,model.q);
end
%     model.X_v = reshape(params(startVal:endVal),model.k,model.q);

startVal = endVal+1;
endVal = endVal + model.q*model.k;

if isempty(find(isnan(params(startVal:endVal))))&&isempty(find(isinf(params(startVal:endVal))))
    model.X_u = reshape(params(startVal:endVal),model.k,model.q);
end

% likelihood beta parameters
if model.optimiseBeta
  startVal = endVal + 1;
  endVal = endVal + prod(size(model.beta));
  if ~isstruct(model.betaTransform)
    fhandle = str2func([model.betaTransform 'Transform']);
    
    if isempty(find(isnan(params(startVal:endVal))))&&isempty(find(isinf(params(startVal:endVal))))
        model.beta = fhandle(params(startVal:endVal), 'atox');
    end
  else
      if isfield(model.betaTransform,'transformsettings') && ~isempty(model.betaTransform.transformsettings)
          fhandle = str2func([model.betaTransform.type 'Transform']);
          model.beta = fhandle(params(startVal:endVal), 'atox', model.betaTransform.transformsettings);
      else
          error('vargplvmExtractParam: Invalid transform specified for beta.'); 
      end
  end
end

startVal = endVal + 1;
endVal = endVal + prod(size(model.W(samd,:)));
lend = length(samd);
if isempty(find(isnan(params(startVal:endVal))))&&isempty(find(isinf(params(startVal:endVal))))
        model.W(samd,:) = reshape(params(startVal:endVal),lend,model.J);
end
  
% kernel hyperparameters 
startVal = endVal+1; 
endVal = endVal + model.kern_v.nParams;
if isempty(find(isnan(params(startVal:endVal))))&&isempty(find(isinf(params(startVal:endVal))))
    model.kern_v = kernExpandParam(model.kern_v, params(startVal:endVal));
end

startVal = endVal+1; 
endVal = endVal + model.kern_u.nParams;

if isempty(find(isnan(params(startVal:endVal))))&&isempty(find(isinf(params(startVal:endVal))))
    model.kern_u = kernExpandParam(model.kern_u, params(startVal:endVal));
end
model.nParams = endVal;

% Update statistics
model = vargplvmUpdateStats(model, model.X_v ,model.X_u, samd, samr);
% %%%TEMP: This is not needed, probably. If yes, it should be merged with
% %%%the above code for fixInducing.
% if model.fixInducing
%     model.X_u=model.X(model.inducingIndices, :);
% end
% %%%
