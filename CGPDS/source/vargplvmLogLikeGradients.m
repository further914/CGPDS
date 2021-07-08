function g = vargplvmLogLikeGradients(model ,stepSize, samd, samr)

% VARGPLVMLOGLIKEGRADIENTS Compute the gradients for the variational GPLVM.
%
%	Description:
%
%	G = VARGPLVMLOGLIKEGRADIENTS(MODEL) returns the gradients of the log
%	likelihood with respect to the parameters of the GP-LVM model and
%	with respect to the latent positions of the GP-LVM model.
%	 Returns:
%	  G - the gradients of the latent positions (or the back
%	   constraint's parameters) and the parameters of the GP-LVM model.
%	 Arguments:
%	  MODEL - the FGPLVM structure containing the parameters and the
%	   latent positions.
%
%	[GX, GPARAM] = VARGPLVMLOGLIKEGRADIENTS(MODEL) returns the gradients
%	of the log likelihood with respect to the parameters of the GP-LVM
%	model and with respect to the latent positions of the GP-LVM model
%	in seperate matrices.
%	 Returns:
%	  GX - the gradients of the latent positions (or the back
%	   constraint's parameters).
%	  GPARAM - gradients of the parameters of the GP-LVM model.
%	 Arguments:
%	  MODEL - the FGPLVM structure containing the parameters and the
%	   latent positions.
%	
%	
%	


%	Copyright (c) 2009, 2010, 2011 Michalis K. Titsias
%	Copyright (c) 2009, 2010, 2011 Mauricio Alvarez
%	Copyright (c) 2009, 2010, 2011 Neil D. Lawrence
%	Copyright (c) 2010, 2011 Andreas Damianou
% 	vargplvmLogLikeGradients.m SVN version 1767
% 	last update 2011-11-21T23:01:18.331219Z
% SEEALSO : vargplvmLogLikelihood, vargplvmCreate, modelLogLikeGradients



% The gradient of the kernel of the dynamics (e.g. temporal prior)
gDynKern = [];

% KL divergence terms: If there are no dynamics, the only params w.r.t to
% which we need derivatives are the var. means and covars. With dynamics
% (see below) we also need derivs. w.r.t theta_t (dyn. kernel's params).
if ~isfield(model, 'dynamics') || isempty(model.dynamics)
    if ~(isfield(model, 'onlyLikelihood') && model.onlyLikelihood)
        gVarmeansKL = - model.vardist.means(:)';
        % !!! the covars are optimized in the log space
        gVarcovsKL = 0.5 - 0.5*model.vardist.covars(:)';
    else
        gVarmeansKL = 0;
        gVarcovsKL = 0;
    end
end

%fprintf(' %d \n',sum([gVarmeansKL gVarcovsKL]));%%%%%TEMP

% Likelihood terms (coefficients)
[gK_vv , gK_uu, gPsi0, gPsi1, gPsi2,gPsi3, gPsi4, gPsi5, g_Lambda, gBeta, gW , tmpV_v ,tmpV_u] = vargpCovGrads(model ,stepSize, samd, samr);


% Get (in three steps because the formula has three terms) the gradients of
% the likelihood part w.r.t the data kernel parameters, variational means
% and covariances (original ones). From the field model.vardist, only
% vardist.means and vardist.covars and vardist.lantentDimension are used.
[gKern0, gVarmeans0, gVarcovs0, gInd0] = kernVardistPsi1Gradient(model.kern_v, model.vardist, model.X_v, gPsi0');
[gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern_u, model.vardist, model.X_u, gPsi1');

[gKern2, gVarmeans2, gVarcovs2] = kernVardistPsi0Gradient(model.kern_v, model.vardist, gPsi2);
[gKern3, gVarmeans3, gVarcovs3] = kernVardistPsi0Gradient(model.kern_u, model.vardist, gPsi3);

[gKern4, gVarmeans4, gVarcovs4, gInd4] = kernVardistPsi2Gradient(model.kern_v, model.vardist, model.X_v, gPsi4);
[gKern5, gVarmeans5, gVarcovs5, gInd5] = kernVardistPsi2Gradient(model.kern_u, model.vardist, model.X_u, gPsi5);

gKern6 = kernGradient(model.kern_v, model.X_v, gK_vv);
gKern7 = kernGradient(model.kern_u, model.X_u, gK_uu);

% At this point, gKern gVarmeansLik and gVarcovsLik have the derivatives for the
% likelihood part. Sum all of them to obtain the final result.
gKern_v = gKern0 + gKern2 + gKern4 + gKern6;
gKern_u = gKern1 + gKern3 + gKern5 + gKern7;
gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2 + gVarmeans3 + gVarmeans4 + gVarmeans5;


% if strcmp(model.kern.type, 'rbfardjit')
%     % different derivatives for the variance, which is super-numerically stable for
%     % this particular kernel
%     gKern(1) = 0.5*model.d*( - model.k+ sum(sum(model.invLat.*model.invLat))/model.beta - model.beta*(model.Psi0-model.TrC)  )...
%                     + 0.5*tmpV;
% end

% if strcmp(model.kern.type, 'rbfardjit')
%     % different derivatives for the variance, which is super-numerically stable for
%     % this particular kernel
%     if model.learnSigmaf == 1
%         gKern(1) = 0.5*model.d*( - model.k+ sum(sum(model.invLat.*model.invLat))/model.beta - model.beta*(model.Psi0-model.TrC)  )...
%             + 0.5*tmpV;
%         
%         if ~isstruct(model.kern.transforms(1))
%             fhandle = str2func([model.kern.transform(1) 'Transform']);
%             gKern(1) = gKern(1).*fhandle(model.kern.variance, 'gradfact');
%         else
%             fhandle = str2func([model.kern.transforms(1).type 'Transform']);
%             if ~isfield(model.kern.transforms(1), 'transformsettings')
%                 gKern(1) = gKern(1).*fhandle(model.kern.variance, 'gradfact');
%             else
%                 gKern(1) = gKern(1).*fhandle(model.kern.variance, 'gradfact', model.kern.transforms(1).transformsettings);
%             end
%         end
%     else
%         gKern(1) = 0;
%     end
% end

%%% Compute Kvv and Kuu Gradients with respect to X_v and X_u %%%
gKX_v = kernGradX(model.kern_v, model.X_v, model.X_v);
gKX_u = kernGradX(model.kern_u, model.X_u, model.X_u);

% The 2 accounts for the fact that covGrad is symmetric
gKX_v = gKX_v*2;
dgKX_v = kernDiagGradX(model.kern_v, model.X_v);
for i = 1:model.k
    gKX_v(i, :, i) = dgKX_v(i, :);
end

gKX_u = gKX_u*2;
dgKX_u = kernDiagGradX(model.kern_u, model.X_u);
for i = 1:model.k
    gKX_u(i, :, i) = dgKX_u(i, :);
end


% Allocate space for gX_u
gX_u = zeros(model.k, model.q);
% Compute portion associated with gK_u
for i = 1:model.k
    for j = 1:model.q
        gX_u(i, j) = gKX_u(:, j, i)'*gK_uu(:, i);%gKuu Q*M
    end
end

% Allocate space for gX_u
gX_v = zeros(model.k, model.q);
% Compute portion associated with gK_u
for i = 1:model.k
    for j = 1:model.q
        gX_v(i, j) = gKX_v(:, j, i)'*gK_vv(:, i);
    end
end

% This should work much faster
%gX_u2 = kernKuuXuGradient(model.kern, model.X_u, gK_uu);

%sum(abs(gX_u2(:)-gX_u(:)))
%pause

gInd_v = gInd0 + gInd4 + gX_v(:)';
gInd_u = gInd1 + gInd5 + gX_u(:)';

% If the inducing points are fixed (tied to the latent points) then
% X_u=K_t*dynamics.vardist.means and the derivatives w.r.t theta_t must be
% amended with the appropriate partial derivatives. gInd must be passed,
% in that case, as an argument to the function which calculates the
% derivatives for the reparametrized quantities.
if isfield(model, 'fixInducing') & model.fixInducing
    gIndRep_v = gInd_v;
    gIndRep_u = gInd_u;
else
    gIndRep_v=[];
    gIndRep_u=[];
end


% If we only want to exclude the derivatives for the variational
% distribution, the following big block will be skipped.
if ~(isfield(model, 'onlyKernel') && model.onlyKernel)
    if isfield(model, 'dynamics') && ~isempty(model.dynamics)
        % Calculate the derivatives for the reparametrized variational and Kt parameters.
        % The formulae for these include in a mixed way the derivatives of the KL
        % term w.r.t these., so gVarmeansKL and gVarcovsKL are not needed now. Also
        % the derivatives w.r.t kernel parameters also require the derivatives of
        % the likelihood term w.r.t the var. parameters, so this call must be put
        % in this part.
        
        % For the dynamical GPLVM further the original covs. must be fed,
        % before amending with the partial derivative due to exponing to enforce
        % positiveness.
        gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2 + gVarcovs3 + gVarcovs4 + gVarcovs5;
        [gVarmeans gVarcovs gDynKern] = modelPriorReparamGrads(model.dynamics, gVarmeansLik, gVarcovsLik, gIndRep_v , gIndRep_u);
        % Variational variances are positive: Now that the final covariances
        % are obtained we amend with the partial derivatives due to the
        % exponential transformation to ensure positiveness.
         if ~isfield(model, 'notransform') || (isfield(model,'notransform') && model.notransform == false)
       		 gVarcovs = (gVarcovs(:).*model.dynamics.vardist.covars(:))';
    	end
    else
        % For the non-dynamical GPLVM these cov. derivatives are the final, so
        % it is time to amend with the partial derivative due to exponing them
        % to force posigiveness.
        gVarcovs0 = (gVarcovs0(:).*model.vardist.covars(:))';
        gVarcovs1 = (gVarcovs1(:).*model.vardist.covars(:))';
        gVarcovs2 = (gVarcovs2(:).*model.vardist.covars(:))';
        gVarcovs3 = (gVarcovs3(:).*model.vardist.covars(:))';
        gVarcovs4 = (gVarcovs4(:).*model.vardist.covars(:))';
        gVarcovs5 = (gVarcovs5(:).*model.vardist.covars(:))';
        
        gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2 + gVarcovs3 + gVarcovs4 + gVarcovs5;
        gVarmeans = gVarmeansLik + gVarmeansKL;
        %gVarcovsLik = (gVarcovsLik(:).*model.vardist.covars(:))';
        gVarcovs = gVarcovsLik + gVarcovsKL;
    end
else
    gVarmeans = [];
    gVarcovs = [];
    gDynKern = [];
end



%{
    % %%% Compute Gradients with respect to X_u %%%
    % gKX = kernGradX(model.kern, model.X_u, model.X_u);
    %
    % % The 2 accounts for the fact that covGrad is symmetric
    % gKX = gKX*2;
    % dgKX = kernDiagGradX(model.kern, model.X_u);
    % for i = 1:model.k
    %   gKX(i, :, i) = dgKX(i, :);
    % end
    %
    % % Allocate space for gX_u
    % gX_u = zeros(model.k, model.q);
    % % Compute portion associated with gK_u
    % for i = 1:model.k
    %   for j = 1:model.q
    %     gX_u(i, j) = gKX(:, j, i)'*gK_uu(:, i);
    %   end
    % end
    %
    % % This should work much faster
    % %gX_u2 = kernKuuXuGradient(model.kern, model.X_u, gK_uu);
    %
    % %sum(abs(gX_u2(:)-gX_u(:)))
    % %pause
    %
    % gInd = gInd1 + gInd2 + gX_u(:)';
%}



%%% TEMP:  the following needs some more testing...
% If fixInducing is true then the inducing points are not optimised but are
% rather reparametrized as X_u=X=mu (static case). Then, there is no derivative
% for X_u any more but the one for mu must be amended by adding the partial
% derivatives of the inducing points (which is just gInd because X_u is
% mapped to mu with the identity function whose deriv. is 1). In the
% dynamics case X_u=X=mu=Kt*mu_bar so we further need to amend with
% d mu/ d mu_bar = K_t because we need
% d F/ d mu_bar instead of d F/ d mu.
if isfield(model, 'fixInducing') & model.fixInducing
    % If there are dynamics the derivative must further be amended with the
    % partial deriv. due to the mean reparametrization.
    if isfield(model, 'dynamics') && ~isempty(model.dynamics)
        gInd_v = reshape(gInd_v,model.k,model.q);
        %gInd = gInd' * model.dynamics.Kt;
        gInd_v =  model.dynamics.Kt * gInd_v;
        gInd_v = gInd_v(:)';
        
        gInd_u = reshape(gInd_u,model.k,model.q);
        %gInd = gInd' * model.dynamics.Kt;
        gInd_u =  model.dynamics.Kt * gInd_u;
        gInd_u = gInd_u(:)';
    end
    %gVarmeans(model.inducingIndices, :) = gVarmeans(model.inducingIndices,
    %:) + gInd; % This should work AFTER reshaping the matrices...but here
    %we use all the indices anyway.
    gVarmeans = gVarmeans + gInd_v + gInd_u;
    gInd_v = []; % Inducing points are not free variables anymore, they dont have derivatives on their own.
    gInd_u = [];
end



gVar = [gVarmeans gVarcovs];

% gVarmeans = gVarmeans0 + gVarmeans1 + gVarmeans2 + gVarmeansKL;
% gVarcovs = gVarcovs0 + gVarcovs1 + gVarcovs2 + gVarcovsKL;



if isfield(model.vardist,'paramGroups')
    gVar = gVar*model.vardist.paramGroups;
end


% If we only want to exclude the derivatives for the variational
% distribution, the following big block will be skipped.
if ~(isfield(model, 'onlyKernel') && model.onlyKernel)
    % It may better to start optimize beta a bit later so that
    % so that the rest parameters can be initialized
    % (this could help to escape from the trivial local
    % minima where the noise beta explains all the data)
    % The learnBeta option deals with the above.
    
    % This constrains the variance of the dynamics kernel to one
    % (This piece of code needs to be done in better way with unit variance dynamic
    %  kernels. The code below also will only work for rbf dynamic kernel)
    % Assume that the base rbf/matern/etc kernel is first in the compound
    % structure
    if isfield(model, 'dynamics') && ~isempty(model.dynamics)
        if strcmp(model.dynamics.kern.comp{1}.type,'rbf') || strcmp(model.dynamics.kern.comp{1}.type,'matern32') || strcmp(model.dynamics.kern.comp{1}.type,'rbfperiodic') || strcmp(model.dynamics.kern.comp{1}.type,'rbfperiodic2')
            if ~isfield(model.dynamics, 'learnVariance') || ~model.dynamics.learnVariance
                gDynKern(2) = 0;
            end
        end
        
        
        %___NEW: assume that the second rbf/matern etc kernel is last in the
        %compound kernel
        %if numel(model.dynamics.kern.comp) > 3
        if isfield(model.dynamics, 'learnSecondVariance') && ~model.dynamics.learnSecondVariance   %%%%% NEW
            gDynKern(end) = 0;
        end
        %end
        %___
    end
end
%% added by Jingjing Fei
% if isfield(model, 'onlyVariationalParam') && model.onlyVariationalParam
%     gBeta=0*gBeta;
%     gKern_u=0*gKern_u;
%     gKern_v=0*gKern_v;
%     gDynKern = 0*gDynKern;
%     gW = 0*gW;
% end
% 
% if isfield(model, 'onlyModelParam') && model.onlyModelParam
%     gVar = 0*gVar;
%     gInd_v = 0*gInd_v;
%     gInd_u = 0*gInd_u;
% end

% if isfield(model, 'optiLocalParam') && model.optiLocalParam
%     gBeta=0*gBeta;
%     gKern_u=0*gKern_u;
%     gKern_v=0*gKern_v;
%     gDynKern = 0*gDynKern;
%     gInd_v = 0*gInd_v;
%     gInd_u = 0*gInd_u;
% end
% 
% if isfield(model, 'optiGlobalParam') && model.optiGlobalParam
%     gVar = 0*gVar;
%     gW = 0*gW;
% end




if isempty(gVar) 
    gVar = zeros(1,model.q*model.N*2); 
end
if isempty(gInd_v) 
    gInd_v = 0*model.X_v(:); 
end
if isempty(gInd_u) 
    gInd_u = 0*model.X_u(:); 
end
if isempty(gDynKern)
    gDynKern= zeros(1,model.dynamics.kern.nParams); 
end
%%
% In case we are in the phase where the vardistr. is initialised (see above
% for the variance of the kernel), beta is kept fixed. For backwards
% compatibility this can be controlled either with the learnBeta field or
% with the initVardist field. The later overrides the first.
if isfield(model, 'learnBeta') && model.learnBeta
    gBetaFinal = gBeta;
else
    gBetaFinal = 0*gBeta;
end
if isfield(model, 'initVardist')
    if model.initVardist == 1
        gBetaFinal = 0*gBeta;
    else
        gBetaFinal = gBeta;
    end
end

% At this point, gDynKern will be [] if there are no dynamics.
% gVar=0*gVar;
% gDynKern=0*gDynKern;
% gInd=0*gInd;
% gBetaFinal=0*gBetaFinal;
% gKern=0*gKern;
jg = 1;
% indexPresent = setdiff([1:model.d],samd);
% gWmean = mean(gW(samd,:),1); %1*model.J
% gW(indexPresent,:) = repmat(gWmean,length(indexPresent),1);

g = [jg*gVar jg*gDynKern jg*gInd_v jg*gInd_u 0*gBetaFinal 0*gW(:)' jg*gKern_v gKern_u];   
% gNozero = find(g);
% fprintf('optimised parameters from %d to %d\n',gNozero(1),gNozero(length(gNozero)));
% fprintf('gInd_u is %f\n',gInd_u*gInd_u');
% fprintf('gKern_u is %f\n',gKern_u*gKern_u');
% fprintf('gBetaFinal is %f\n',gBetaFinal*gBetaFinal');


% fprintf(' gVar is %f\n',gVar*gVar');
% fprintf(' gBetaFinal is %f\n',gBetaFinal*gBetaFinal');
% gpart=zeros(1,numel([gDynKern gInd gKern gBetaFinal]));
% g = [gVar gpart];
% if model.learnBeta == 1
%     % optimize all parameters including beta
%     g = [gVar gDynKern gInd gKern gBeta];
% else
%     % keep beta fixed
%     g = [gVar gDynKern gInd gKern 0*gBeta];
% end


%if model.kern.comp{1}.variance > 100 | model.kern.comp{1}.variance < 0.001
%     model.kern.comp{1}
%end

% delete afterwards
%g = [gVarmeans2 gVarcovs2 (gX_u(:)'+ gInd2) (gKern3+gKern2) 0*gBeta];

end


function [gK_vv , gK_uu, gPsi0, gPsi1, gPsi2,gPsi3, gPsi4, gPsi5, g_Lambda, gBeta, gW , tmpV_v ,tmpV_u] = vargpCovGrads(model ,stepSize, samd, samr)

sigm = 1/model.beta; % beta^-1
lend = length(samd);
lenr = length(samr);

tempMatrix = calculateMatrix( model,samd,samr );

%with svi
gBeta_v = 0;

    gPsi0 = model.beta * model.m * model.B_v'; %Npart*M
    gPsi0 = gPsi0'; % because it is passed to "kernVardistPsi1Gradient" as gPsi1'...M*N
    gPsi2 = -0.5 * model.beta * model.d; %1*1
    gPsi4 = (model.beta/2) * model.T1_v; %30*30s
    gK_vv = 0.5 * (model.T1_v - (model.beta * model.d) * model.invLmT_v * model.C_v * model.invLm_v);%30*30
    PLm_v = model.invLatT_v*model.P_v; %30*59
    tmpV_v = sum(sum(PLm_v.*PLm_v)); %1*1
    gBeta_v = 0.5*(model.d*(model.TrC_v + (model.N-model.k)*sigm -model.Psi2 ) ...
        - model.TrYY + model.TrPP_v + model.d/lend*tempMatrix.TrPP_u  ...
        + (1/model.beta^2 * model.d * sum(sum(model.invLat_v.*model.invLat_v)))...
        + sigm*tmpV_v);

%     gPsi0 = model.beta * model.m(:,samd) * tempMatrix.B_v'; %Npart*M
%     gPsi0 = gPsi0'; % because it is passed to "kernVardistPsi1Gradient" as gPsi1'...M*N
%     gPsi2 = -0.5 * model.beta * lend; %1*1
%     gPsi4 = (model.beta/2) * tempMatrix.T1_v; %30*30s
%     gK_vv = 0.5 * (tempMatrix.T1_v - (model.beta * lend) * model.invLmT_v * tempMatrix.C_v * model.invLm_v);%30*30
%     PLm_v = tempMatrix.invLatT_v*tempMatrix.P_v; %30*59
%     tmpV_v = sum(sum(PLm_v.*PLm_v)); %1*1
%     gBeta_v = 0.5*(lend*(tempMatrix.TrC_v + (model.N-model.k)*sigm -model.Psi2 ) ...
%         - tempMatrix.TrYYdr + tempMatrix.TrPP_v + tempMatrix.TrPP_u  ...
%         + (1/model.beta^2 * lend * sum(sum(tempMatrix.invLat_v.*tempMatrix.invLat_v)))...
%         + sigm*tmpV_v);
%     gPsi0 = model.d/lend * gPsi0;
%     gPsi2 = model.d/lend * gPsi2;
%     gPsi4 = model.d/lend * gPsi4;
%     gK_vv = model.d/lend * gK_vv;
%     gBeta_v = model.d/lend * gBeta_v;


gPsi1_temp = zeros(model.N,model.k);
gPsi3 = 0;
gK_uu = zeros(model.k,model.k);
gW = zeros(model.d,model.J);
gBeta_u = 0;

for d = 1:lend
    dd = samd(d);
    for j = 1:model.J
            
            gPsi1_temp = gPsi1_temp + model.beta * model.W(dd,j)^2 * model.m(:,dd) * tempMatrix.B_u{d,j}';%N*M
        
            gPsi3 = gPsi3 - 0.5 * model.beta * model.W(dd,j)^2;

            Tb_u = (1/model.beta) * (tempMatrix.P1_u{d,j}' * tempMatrix.P1_u{d,j});
            Tb_u = Tb_u + (model.W(dd,j)^2 * tempMatrix.B_u{d,j} * tempMatrix.B_u{d,j}');  
    %         gK_uu = gK_uu + 0.5 * (model.invK_uu - Tb_u - (model.beta) * model.invLmT_u * model.C_u{d,j} * model.invLm_u);%30*30
            gK_uu = gK_uu + 0.5 * (model.invK_uu - Tb_u - (model.beta) * model.invLmT_u * tempMatrix.C_u{d,j} * model.invLm_u);%30*30

            PLm_u = tempMatrix.invLatT_u{d,j}*tempMatrix.P_u{d,j}; %30*59
            tmpV_u = model.W(dd,j)^2 * sum(sum(PLm_u.*PLm_u)); %1*1

            gBeta_u = gBeta_u + 0.5*( tempMatrix.TrC_u(d,j) + (-model.k)*sigm - model.W(dd,j)^2*model.Psi3 ...
                + (1/model.beta^2 * (sum(sum(tempMatrix.invLat_u{d,j}.*tempMatrix.invLat_u{d,j}))))...
                + sigm*tmpV_u);

            gW(dd,j) = -trace(model.W(dd,j)*tempMatrix.P1_u{d,j}'*tempMatrix.P1_u{d,j}*model.Psi5) ...
                + model.beta*model.W(dd,j)*model.m(:,dd)'*model.Psi1*tempMatrix.B_u{d,j}...
                -model.beta*model.W(dd,j)^3*trace(tempMatrix.B_u{d,j}*tempMatrix.B_u{d,j}'*model.Psi5)...
                -model.beta*model.W(dd,j)*model.Psi3 + model.beta*trace(model.W(dd,j)*model.Psi5*model.invK_uu);
            
            %approximate gW
%             if d > 1
%                 lastdd = samd(d-1);
%                 nowdd = samd(d);
%                 gW(lastdd+1:nowdd-1,j) = repmat((gW(nowdd,j) + gW(lastdd,j))/2,length(lastdd+1:nowdd-1),1);
%             elseif d == 1
%                 gW(1:dd-1,j) = repmat(gW(dd,j),dd-1,1);
%             end
    end
end

gPsi5 = (model.beta/2) * tempMatrix.T1_u; %30*30


    
   
    
    % gPsi1 = model.beta * model.W^2 * model.m * model.B_u'; %30*359
    gPsi1 = model.d/lend*gPsi1_temp'; % because it is passed to "kernVardistPsi1Gradient" as gPsi1'...
    gPsi3 = model.d/lend*gPsi3;
    gPsi5 = model.d/lend * gPsi5;
    gK_uu = model.d/lend*gK_uu;
    gBeta_u = model.d/lend*gBeta_u;
    % gPsi3 = -0.5 * model.beta * model.W^2 * model.d; %1*1
    

gBeta = gBeta_v + gBeta_u;


%%%%TEMP
%{
    load TEMPbetaGradTrC;
    TEMPbetaGradTrC = [TEMPbetaGradTrC model.d*0.5*model.TrC];
    save 'TEMPbetaGradTrC.mat' 'TEMPbetaGradTrC';

    load TEMPbetaGradNksigm;
    TEMPbetaGradNksigm=[TEMPbetaGradNksigm model.d*0.5*(model.N-model.k)*sigm];
    save 'TEMPbetaGradNksigm.mat' 'TEMPbetaGradNksigm';

    load TEMPbetaGradPsi0;
    TEMPbetaGradPsi0=[TEMPbetaGradPsi0 (-0.5*model.d*model.Psi0)];
    save 'TEMPbetaGradPsi0.mat' 'TEMPbetaGradPsi0';

    load TEMPbetaGradTrPP;
    TEMPbetaGradTrPP=[TEMPbetaGradTrPP 0.5*model.TrPP];
    save 'TEMPbetaGradTrPP.mat' 'TEMPbetaGradTrPP';

    load TEMPbetaGradLat;
    TEMPbetaGradLat=[TEMPbetaGradLat (1/(model.beta^2)) * model.d * sum(sum(model.invLat.*model.invLat))*0.5];
    save 'TEMPbetaGradLat.mat' 'TEMPbetaGradLat';

    load TEMPbetaGradPlm;
    TEMPbetaGradPlm=[TEMPbetaGradPlm sigm*sum(sum(PLm.*PLm))*0.5];
    save 'TEMPbetaGradPlm.mat' 'TEMPbetaGradPlm';
%}
%%%%%


%gBeta = 0.5*(model.d*(model.TrC + (model.N-model.k)*sigm -model.Psi0) ...
%	- model.TrYY + model.TrPP ...
%	+ sigm * sum(sum(model.K_uu .* model.Tb)));

if ~isstruct(model.betaTransform)
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact');
else
    fhandle = str2func([model.betaTransform.type 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact', model.betaTransform.transformsettings);
end


g_Lambda = repmat(-0.5*model.beta*model.d, 1, model.N);

clear tempMatrix;
end
