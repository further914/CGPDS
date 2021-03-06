function [g1, g2] = lfmaXlfmaKernGradient(lfmKern1, lfmKern2, t1, t2, covGrad, meanVector)

% LFMAXLFMAKERNGRADIENT Compute a cross gradient between a LFMA and a LFMA.
%
%	Description:
%
%	[G1, G2] = LFMAXLFMAKERNGRADIENT(LFMKERN1, LFMKERN2, T1, COVGRAD)
%	computes cross gradient of parameters of a cross kernel between two
%	lfm kernels, both lfm kernels correspond to the accelerations. with
%	the multiple output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see lfmKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see lfmKernExtractParam.
%	 Arguments:
%	  LFMKERN1 - the kernel structure associated with the first LFM
%	   kernel (acceleration).
%	  LFMKERN2 - the kernel structure associated with the second LFM
%	   kernel (acceleration).
%	  T1 - inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
%	[G1, G2] = LFMAXLFMAKERNGRADIENT(LFMKERN1, LFMKERN2, T1, T2,
%	COVGRAD) computes cross gradient of parameters of a cross kernel
%	between two lfm kernels, both lfm kernels correspond to the
%	accelerations. with the multiple output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see lfmKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see lfmKernExtractParam.
%	 Arguments:
%	  LFMKERN1 - the kernel structure associated with the first LFM
%	   kernel (acceleration).
%	  LFMKERN2 - the kernel structure associated with the second LFM
%	   kernel (acceleration).
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
%	[G1, G2] = LFMAXLFMAKERNGRADIENT(LFMKERN1, LFMKERN2, T1, T2,
%	COVGRAD, MEANVEC) computes cross gradient of parameters of a cross
%	kernel between two lfm kernels. Both kernels correspond to
%	accelerations. It is supposed to be used together with the multiple
%	output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see lfmKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see lfmKernExtractParam.
%	 Arguments:
%	  LFMKERN1 - the kernel structure associated with the first LFM
%	   kernel (acceleration).
%	  LFMKERN2 - the kernel structure associated with the second LFM
%	   kernel (acceleration).
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%	  MEANVEC - precomputed factor that is used for the switching
%	   dynamical latent force model.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, LFMVKERNPARAMINIT


%	Copyright (c) 2010 Mauricio Alvarez
% 	lfmaXlfmaKernGradient.m SVN version 1519
% 	last update 2011-07-22T13:08:31.000000Z

subComponent = false; % This is just a flag that indicates if this kernel is part of a bigger kernel (SDLFM)

if nargin == 4
    covGrad = t2;
    t2 = t1;
elseif nargin == 6
    subComponent = true;
    if numel(meanVector)>1
        if size(meanVector,1) == 1,
            if size(meanVector, 2)~=size(covGrad, 2)
                error('The dimensions of meanVector don''t correspond to the dimensions of covGrad')
            end
        else
            if size((meanVector'), 2)~=size(covGrad,2)
                error('The dimensions of meanVector don''t correspond to the dimensions of covGrad')
            end
        end
    else
        if numel(t1)==1 && numel(t2)>1
            % matGrad will be row vector and so should be covGrad
            dimcovGrad = length(covGrad);
            covGrad = reshape(covGrad, [1 dimcovGrad]);
        elseif numel(t1)>1 && numel(t2)==1
            % matGrad will be column vector and sp should be covGrad
            dimcovGrad = length(covGrad);
            covGrad = reshape(covGrad, [dimcovGrad 1]);
        end
    end  
end

if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end
if lfmKern1.inverseWidth ~= lfmKern2.inverseWidth
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

% Parameters of the simulation (in the order providen by kernExtractParam)

m = [lfmKern1.mass lfmKern2.mass];                  % Par. 1
D = [lfmKern1.spring lfmKern2.spring];              % Par. 2
C = [lfmKern1.damper lfmKern2.damper];              % Par. 3
sigma2 = 2/lfmKern1.inverseWidth;                   % Par. 4
sigma = sqrt(sigma2);
S = [lfmKern1.sensitivity lfmKern2.sensitivity];    % Par. 5

alpha = C./(2*m);
omega = sqrt(D./m-alpha.^2);

% Initialization of vectors and matrices

g1 = zeros(1,5);
g2 = zeros(1,5);

% Precomputations
computeH  = cell(4,1);
computeUpsilonMatrix = cell(2,1);
computeUpsilonVector = cell(2,1);
gradientUpsilonMatrix = cell(4,1);
gradientUpsilonVector = cell(4,1);
gamma1_p = alpha(1) + j*omega(1);
gamma1_m = alpha(1) - j*omega(1);
gamma2_p = alpha(2) + j*omega(2);
gamma2_m = alpha(2) - j*omega(2);
preExp1 = zeros(length(t1),2);
preExp2 = zeros(length(t2),2);
preExpg1 = zeros(length(t1),2);
preExpgg1 = zeros(length(t1),2);
preExpg2 = zeros(length(t2),2);
preExpgg2 = zeros(length(t2),2);
preExpt1 = zeros(length(t1),2);
preExpt2 = zeros(length(t2),2);
preGamma(1) = gamma1_p + gamma2_p;
preGamma(2) = gamma1_p + gamma2_m;
preGamma(3) = gamma1_m + gamma2_p;
preGamma(4) = gamma1_m + gamma2_m;
preGamma2 = preGamma.^2;
preConst = 1./preGamma;
preConst2 = 1./(preGamma2);
preFactors(1) = preConst(2) - preConst(1);
preFactors(2) = preConst(3) - preConst(4);
preFactors(3) = preConst(3) - preConst(1);
preFactors(4) = preConst(2) - preConst(4);
preFactors2(1) = -preConst2(2) + preConst2(1);
preFactors2(2) = -preConst2(3) + preConst2(4);
preFactors2(3) = -preConst2(3) + preConst2(1);
preFactors2(4) = -preConst2(2) + preConst2(4);
preExp1(:,1) = exp(-gamma1_p*t1);
preExp1(:,2) = exp(-gamma1_m*t1);
preExp2(:,1) = exp(-gamma2_p*t2);
preExp2(:,2) = exp(-gamma2_m*t2);
preExpg1(:,1) = (2*gamma1_p)*preExp1(:,1);
preExpg1(:,2) = (2*gamma1_m)*preExp1(:,2);
preExpg2(:,1) = (2*gamma2_p)*preExp2(:,1);
preExpg2(:,2) = (2*gamma2_m)*preExp2(:,2);
preExpgg1(:,1) = (gamma1_p^2)*preExp1(:,1);
preExpgg1(:,2) = (gamma1_m^2)*preExp1(:,2);
preExpgg2(:,1) = (gamma2_p^2)*preExp2(:,1);
preExpgg2(:,2) = (gamma2_m^2)*preExp2(:,2);
preExpt1(:,1) = t1.*preExpgg1(:,1);
preExpt1(:,2) = t1.*preExpgg1(:,2);
preExpt2(:,1) = t2.*preExpgg2(:,1);
preExpt2(:,2) = t2.*preExpgg2(:,2);
[computeH{1}, computeUpsilonMatrix{1}, computeUpsilonMatrixLocal{1}] =  lfmComputeH3AA(gamma1_p, gamma1_m, sigma2, t1,t2,preFactors([1 2]), 0);
[computeH{2}, computeUpsilonMatrix{2}, computeUpsilonMatrixLocal{2}] =  lfmComputeH3AA(gamma2_p, gamma2_m, sigma2, t2,t1,preFactors([3 4]), 1);
[computeH{3}, computeUpsilonVector{1}, computeUpsilonVectorLocal{1}] =  lfmComputeH4AA(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExpgg2, 0 );
[computeH{4}, computeUpsilonVector{2}, computeUpsilonVectorLocal{2}] =  lfmComputeH4AA(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExpgg1, 1 );
preKernel = ( computeH{1} + computeH{2}.' + computeH{3} + computeH{4}.');

% Precompute derivatives
gradientUpsilonMatrix{1} = lfmaaGradientUpsilonMatrix(gamma1_p,sigma2,t1, t2, 0, computeUpsilonMatrixLocal{1}{1});
gradientUpsilonMatrix{2} = lfmaaGradientUpsilonMatrix(gamma1_m,sigma2,t1, t2, 0, computeUpsilonMatrixLocal{1}{2});
gradientUpsilonMatrix{3} = lfmaaGradientUpsilonMatrix(gamma2_p,sigma2,t2, t1, 1, computeUpsilonMatrixLocal{2}{1});
gradientUpsilonMatrix{4} = lfmaaGradientUpsilonMatrix(gamma2_m,sigma2,t2, t1, 1, computeUpsilonMatrixLocal{2}{2});
gradientUpsilonVector{1} = lfmapGradientUpsilonVector(gamma1_p,sigma2,t1, computeUpsilonVectorLocal{1}{1});
gradientUpsilonVector{2} = lfmapGradientUpsilonVector(gamma1_m,sigma2,t1, computeUpsilonVectorLocal{1}{2});
gradientUpsilonVector{3} = lfmapGradientUpsilonVector(gamma2_p,sigma2,t2, computeUpsilonVectorLocal{2}{1});
gradientUpsilonVector{4} = lfmapGradientUpsilonVector(gamma2_m,sigma2,t2, computeUpsilonVectorLocal{2}{2});

if lfmKern1.isNormalised
    K0 =  lfmKern1.sensitivity*lfmKern2.sensitivity/(8*sqrt(2)*lfmKern1.mass*lfmKern2.mass*prod(omega));
    K02 = 1/(8*sqrt(2)*prod(m)*prod(omega));
else
    K0 =  sigma*sqrt(pi)*lfmKern1.sensitivity*lfmKern2.sensitivity/(8*lfmKern1.mass*lfmKern2.mass*prod(omega)); 
    K02 = sigma*sqrt(pi)/(8*prod(m)*prod(omega));
end

% Gradient with respect to m, D and C
for ind_theta = 1:3 % Parameter (m, D or C)
    for ind_par = 0:1 % System (1 or 2)
        % Choosing the right gradients for m, omega, gamma1 and gamma2
        switch ind_theta
            case 1  % Gradient wrt m
                gradThetaM = [1-ind_par ind_par];
                gradThetaAlpha = -C./(2*(m.^2));
                gradThetaOmega = (C.^2-2*m.*D)./(2*(m.^2).*sqrt(4*m.*D-C.^2));
            case 2  % Gradient wrt D
                gradThetaM = zeros(1,2);
                gradThetaAlpha = zeros(1,2);
                gradThetaOmega = 1./sqrt(4*m.*D-C.^2);
            case 3  % Gradient wrt C
                gradThetaM = zeros(1,2);
                gradThetaAlpha = 1./(2*m);
                gradThetaOmega = -C./(2*m.*sqrt(4*m.*D-C.^2));
        end
        gradThetaGamma1 = gradThetaAlpha + j*gradThetaOmega;
        gradThetaGamma2 = gradThetaAlpha - j*gradThetaOmega;
        % Gradient evaluation
        gradThetaGamma11 = [gradThetaGamma1(1) gradThetaGamma2(1)];
        gradThetaGamma2 = [gradThetaGamma1(2) gradThetaGamma2(2)];
        gradThetaGamma1 = gradThetaGamma11;

        if ~ind_par % ind_par = k or d
            matGrad = K0 ...
                * ( lfmGradientH31( preFactors([1 2]), preFactors2([1 2]), gradThetaGamma1, ...
                gradientUpsilonMatrix{1}, gradientUpsilonMatrix{2}, computeUpsilonMatrix{1}{1}, ...
                computeUpsilonMatrix{1}{2}, 1) + ...
                lfmGradientH32( preGamma2,  gradThetaGamma1, computeUpsilonMatrix{2}{1}, ...
                computeUpsilonMatrix{2}{2}, 1).' + ...
                lfmGradientH41( preGamma, preGamma2, gradThetaGamma1, preExpgg2, ...
                gradientUpsilonVector{1}, gradientUpsilonVector{2}, computeUpsilonVector{1}{1},...
                computeUpsilonVector{1}{2}, 1) + ...
                lfmGradientH42AP(preGamma, preGamma2, gradThetaGamma1, preExpg1, preExpgg1, preExpt1, ...
                computeUpsilonVector{2}{1}, computeUpsilonVector{2}{2}).'...
                - (gradThetaM(1+ind_par)/m(1+ind_par) ...
                + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                *preKernel);
        else % ind_par = r or d'
            matGrad = K0 ...
                * ( lfmGradientH31( preFactors([3 4]), preFactors2([3 4]), gradThetaGamma2, ...
                gradientUpsilonMatrix{3}, gradientUpsilonMatrix{4}, computeUpsilonMatrix{2}{1}, ...
                computeUpsilonMatrix{2}{2}, 1).' + ...
                lfmGradientH32( preGamma2([1 3 2 4]),  gradThetaGamma2, computeUpsilonMatrix{1}{1}, ...
                computeUpsilonMatrix{1}{2}, 1) + ...
                lfmGradientH41( preGamma([1 3 2 4]), preGamma2([1 3 2 4]), gradThetaGamma2, preExpgg1, ...
                    gradientUpsilonVector{3}, gradientUpsilonVector{4}, computeUpsilonVector{2}{1},...
                    computeUpsilonVector{2}{2}, 1).' + ...
                lfmGradientH42AP(preGamma([1 3 2 4]), preGamma2([1 3 2 4]), gradThetaGamma2, preExpg2, preExpgg2, preExpt2, ...
                computeUpsilonVector{1}{1}, computeUpsilonVector{1}{2})...
                - (gradThetaM(1+ind_par)/m(1+ind_par) ...
                + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                *preKernel);            
        end

        if subComponent
            if size(meanVector,1) ==1,
                matGrad = matGrad*meanVector;
            else
                matGrad = (meanVector*matGrad).';
            end
        end
         
        % Check the parameter to assign the derivative
        if ind_par == 0
            g1(ind_theta) = sum(sum(matGrad.*covGrad));
        else
            g2(ind_theta) = sum(sum(matGrad.*covGrad));
        end       
    end
end

% Gradients with respect to sigma
if lfmKern1.isNormalised
    matGrad = K0*(lfmGradientSigmaH3AA(gamma1_p, gamma1_m, sigma2, t1, t2, preFactors([1 2]), 0)...
        +  lfmGradientSigmaH3AA(gamma2_p, gamma2_m, sigma2, t2, t1, preFactors([3 4]), 1).'...
        +  lfmGradientSigmaH4AA(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExpgg2)...
        +  lfmGradientSigmaH4AA(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExpgg1).');
else
    matGrad = (prod(S)*sqrt(pi)/(8*prod(m)*prod(omega))) ...
        * (preKernel ...
        + sigma*(lfmGradientSigmaH3AA(gamma1_p, gamma1_m, sigma2, t1, t2, preFactors([1 2]), 0)...
        +  lfmGradientSigmaH3AA(gamma2_p, gamma2_m, sigma2, t2, t1, preFactors([3 4]), 1).'...
        +  lfmGradientSigmaH4AA(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExpgg2)...
        +  lfmGradientSigmaH4AA(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExpgg1).' ));
end

if subComponent
    if size(meanVector,1) ==1,
        matGrad = matGrad*meanVector;
    else
        matGrad = (meanVector*matGrad).';
    end
end

g1(4) = sum(sum(matGrad.*covGrad))*(-(sigma^3)/4);
g2(4) = g1(4);

% Gradients with respect to S


matGrad = K02*preKernel;

if subComponent
    if size(meanVector,1) ==1,
        matGrad = matGrad*meanVector;
    else
        matGrad = (meanVector*matGrad).';
    end
end

g1(5) = sum(sum(S(2)*matGrad.*covGrad));
g2(5) = sum(sum(S(1)*matGrad.*covGrad));


g2(4) = 0; % Otherwise is counted twice, temporarly changed by Mauricio Alvarez

g1 = real(g1);
g2 = real(g2);
