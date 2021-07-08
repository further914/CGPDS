function ll = cgpdsLogLikelihood(model ,stepSize, samd, samr)

% cgpdsLogLikelihood Log-likelihood for a CGPDS
%
%	Description:
%
%	LL = cgpdsLogLikelihood(MODEL) returns the log likelihood for a
%	given CGPDS.
%	 Returns:
%	  LL - the log likelihood of the data given the model.
%	 Arguments:
%	  MODEL - the model for which the log likelihood is to be computed.
%	   The model contains the data for which the likelihood is being
%	   computed in the 'y' component of the structure.


%	Copyright (c) 2009-2011 Michalis K. Titsias
%	Copyright (c) 2009-2011 Neil D. Lawrence
%	Copyright (c) 2010-2011 Andreas Damianou
% 	vargplvmLogLikelihood.m SVN version 1570
% 	last update 2011-08-30T14:57:49.000000Z


% Note: The 'onlyKL' and 'onlyLikelihood' fields can be set by external
% wrappers and cause the function to only calculate the likelihood or the
% KL part of the variational bound.

% Likelihood term

lend = length(samd);
lenr = length(samr);

tempMatrix = calculateMatrix(model,samd,samr);

if ~(isfield(model, 'onlyKL') && model.onlyKL)
        if length(model.beta)==1
            ll = -0.5*(lend*(-(model.N-model.k)*log(model.beta) ...
                + tempMatrix.logDetAt_v ) ...
                - (tempMatrix.TrPP_v- tempMatrix.TrYYdr)*model.beta);
            %if strcmp(model.approx, 'dtcvar')
            ll = ll - 0.5*model.beta*lend*model.Psi2 + 0.5*lend*model.beta*tempMatrix.TrC_v;

%             for d = 1:model.d
%                 for j = 1:model.J
%                     ll = ll - 0.5*(model.logDetAt_u(d,j) - (-model.k)*log(model.beta));
%                     %if strcmp(model.approx, 'dtcvar')
%                     ll = ll - 0.5*model.beta*model.W(d,j)^2*model.Psi3 + 0.5*model.beta*model.TrC_u(d,j);
%                 end
%             end

            ll = ll - 0.5*(sum(sum(tempMatrix.logDetAt_u)) - lend*model.J*(-model.k)*log(model.beta));
            %if strcmp(model.approx, 'dtcvar')
            ll = ll - 0.5*model.beta*sum(sum((model.W(samd,:)).^2))*model.Psi3 + 0.5*model.beta*sum(sum(tempMatrix.TrC_u));
            ll = ll + 0.5*model.beta*tempMatrix.TrPP_u;
        %end
        else
            error('Not implemented variable length beta yet.');
        end
        ll = ll-lend*model.N/2*log(2*pi);
else
    ll=0;
end

ll = model.d/lend*ll;
clear tempMatrix;

% KL divergence term
if ~(isfield(model, 'onlyLikelihood') && model.onlyLikelihood)
    if isfield(model, 'dynamics') && ~isempty(model.dynamics)
        % A dynamics model is being used.
        KLdiv = modelVarPriorBound(model);
    else
        varmeans = sum(sum(model.vardist.means.*model.vardist.means));
        varcovs = sum(sum(model.vardist.covars - log(model.vardist.covars)));
        KLdiv = -0.5*(varmeans + varcovs) + 0.5*model.q*model.N;
    end
else
    KLdiv=0;
end

% Obtain the final value of the bound by adding the likelihood
% and the KL term.
ll = ll + KLdiv;
