function ll = varmgplvmLogLikelihood(model)

% VARGPLVMLOGLIKELIHOOD Log-likelihood for a variational GP-LVM.
%
%	Description:
%
%	LL = VARGPLVMLOGLIKELIHOOD(MODEL) returns the log likelihood for a
%	given GP-LVM model.
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
if ~(isfield(model, 'onlyKL') && model.onlyKL)
    if length(model.beta)==1
        ll = -0.5*(-(model.N*model.d-model.k)*log(model.beta) ...
            + model.logDetAt_u + model.logDetAt_v ...
            - (model.TrPP_u + model.TrPP_v ...
            - model.TrYY)*model.beta);
        %if strcmp(model.approx, 'dtcvar')
        ll = ll - 0.5*model.beta*model.Psi2 - 0.5*W^2*model.beta*model.Psi3  ...
            + 0.5*model.beta*model.TrC_v + 0.5*W^2*model.beta*model.TrC_u;
        %end
    else
        error('Not implemented variable length beta yet.');
    end
    ll = ll-model.d*model.N/2*log(2*pi);
else
    ll=0;
end

%{
% %%%%%%%%TEMP_ (for analysing the "difficult" trace term)
%try
%    load TEMPLikelihoodTrace
%catch exception
%    TEMPLikelihoodTrace=[];
%end
%TEMPLikelihoodTrace=[TEMPLikelihoodTrace ;-0.5*model.beta*model.d*model.Psi0 + 0.5*model.d*model.beta*model.TrC];
%save 'TEMPLikelihoodTrace.mat' 'TEMPLikelihoodTrace';
%%% _TEMP
% load TEMPbetaLikTrC;
% TEMPbetaLikTrC = [TEMPbetaLikTrC 0.5*model.d*model.beta*model.TrC];
% save 'TEMPbetaLikTrC.mat' 'TEMPbetaLikTrC';
%
% load TEMPbetaLikNksigm;
% TEMPbetaLikNksigm=[TEMPbetaLikNksigm (model.d*(-(model.N-model.k)*log(model.beta)))];
% save 'TEMPbetaLikNksigm.mat' 'TEMPbetaLikNksigm';
%
% load TEMPbetaLikPsi0;
% TEMPbetaLikPsi0=[TEMPbetaLikPsi0 (- 0.5*model.beta*model.d*model.Psi0)];
% save 'TEMPbetaLikPsi0.mat' 'TEMPbetaLikPsi0';
%
% load TEMPbetaLikTrPP;
% TEMPbetaLikTrPP=[TEMPbetaLikTrPP 0.5*model.TrPP];
% save 'TEMPbetaLikTrPP.mat' 'TEMPbetaLikTrPP';
%
% load TEMPbetaLikLat;
% TEMPbetaLikLat=[TEMPbetaLikLat -0.5*model.d*model.logDetAt];
% save 'TEMPbetaLikLat.mat' 'TEMPbetaLikLat';
%%%%%%%%%%%%
%}

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
