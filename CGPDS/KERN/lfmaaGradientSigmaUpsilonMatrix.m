function dUpsilonS = lfmaaGradientSigmaUpsilonMatrix(gamma, sigma2, ...
    t1, t2, mode)

% LFMAAGRADIENTSIGMAUPSILONMATRIX Gradient of upsilon matrix aa wrt sigma
%
%	Description:
%
%	UPSILON = LFMAAGRADIENTSIGMAUPSILONMATRIX(GAMMA, SIGMA2, T1, T2,
%	MODE) computes the gradient wrt sigma of a portion of the LFMAA
%	kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  MODE - operation mode, according to the derivative (mode 0,
%	   derivative wrt t1, mode 1 derivative wrt t2)
%	
%
%	See also
%	LFMAACOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio Alvarez
% 	lfmaaGradientSigmaUpsilonMatrix.m SVN version 807
% 	last update 2010-05-28T06:01:33.000000Z

sigma = sqrt(sigma2);
gridt1 = repmat(t1, 1, length(t2));
gridt2 = repmat(t2', length(t1), 1);
timeGrid = gridt1 - gridt2;

dUpsilon = lfmapGradientSigmaUpsilonMatrix(gamma, sigma2, t1, t2, 1-mode);

if mode == 0
    dUpsilonS = (gamma^2)*dUpsilon - (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2).* ...
        ((2/sigma2- (2*timeGrid/sigma2).^2).*((gamma + 2*timeGrid/sigma2).* ...
        (1/sigma - 2*(timeGrid.^2)/sigma^3) + 4*timeGrid/sigma^3) ...
        - (gamma + 2*timeGrid/sigma2).*(-4/sigma^3 + 16*timeGrid.^2/sigma^5)) ...
        - (16/(sqrt(pi)*sigma^6))*timeGrid.*exp(-(timeGrid.^2)./sigma2).* ...
        (5 - 2*timeGrid.^2/sigma2);
else
    dUpsilonS = (gamma^2)*dUpsilon - (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2).* ...
        ((2/sigma2- (2*timeGrid/sigma2).^2).*((gamma + 2*timeGrid/sigma2).* ...
        (1/sigma - 2*(timeGrid.^2)/sigma^3) + 4*timeGrid/sigma^3) ...
        - (gamma + 2*timeGrid/sigma2).*(-4/sigma^3 + 16*timeGrid.^2/sigma^5)) ...
        - (16/(sqrt(pi)*sigma^6))*timeGrid.*exp(-(timeGrid.^2)./sigma2).* ...
        (5 - 2*timeGrid.^2/sigma2) - (2*gamma^2/(sqrt(pi)*sigma))*exp(-gamma*t1)* ...
        (((gamma-2*t2/sigma2).*(1/sigma - 2*t2.^2/sigma^3) - 4*t2/sigma^3).*exp(-t2.^2/sigma2)).';  
end
