function [g, model] = varmgplvmPointLogLikeGradient(model, vardistx, y)

% VARGPLVMPOINTLOGLIKEGRADIENT Log-likelihood gradient for of some points of the GP-LVM.
%
%	Description:
%
%	G = VARGPLVMPOINTLOGLIKEGRADIENT(MODEL, X, Y) returns the gradient
%	of the log likelihood with respect to the latent positions, where
%	the log likelihood is conditioned on the training set. The function
%	works even when we are interested only in one point. See
%	vargplvmPointLogLikelihood for more details.
%	 Returns:
%	  G - the gradient of the log likelihood, conditioned on the
%	   training data, with respect to the latent position.
%	 Arguments:
%	  MODEL - the model for which the gradient computation is being
%	   done.
%	  X - the latent positions where the gradient is being computed.
%	  Y - the positions in data space for which the computation is being
%	   done.
%	
%
%	See also
%	VARGPLVMPOINTLOGLIKELIHOOD, VARGPLVMOPTIMISEPOINT, VAGPLVMSEQUENCELOGLIKEGRADIENT


%	Copyright (c) 2009, 2011 Michalis K. Titsias and Neil D. Lawrence
%	Copyright (c) 2011 Andreas C. Damianou
% 	vargplvmPointLogLikeGradient.m SVN version 1408
% 	last update 2011-07-04T19:08:52.114563Z

if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    dynUsed = 1;
else
    dynUsed = 0;
end

indexMissing = find(isnan(y(1,:)));
indexPresent = setdiff(1:model.d, indexMissing);
% yMs = y(:,indexPresent);
% if ~isempty(yMs)
%    myMs = yMs - repmat(model.bias(indexPresent),size(yMs,1),1);  
%    myMs = myMs./repmat(model.scale(indexPresent),size(yMs,1),1);  
% end
map=y;
map(:,indexPresent)=1;%记录哪里丢失数据，有数据的为1，无数据的为0
map(:,indexMissing)=0;%记录哪里丢失数据，有数据的为1，无数据的为0    
model.map=[ones(model.N,model.d);map];
if dynUsed
    % use a special function
    [g, model] = dynPointLogLikeGradient(model, vardistx, y, indexPresent, indexMissing);
    % g = g(:)'; %%% RE-OPT-CODE-REM
    return;
end



% normalize y exactly as model.m is normalized
my = y - repmat(model.bias(indexPresent),size(y,1),1);
my = my./repmat(model.scale(indexPresent),size(y,1),1);

[gPsi0, gPsi1, gPsi2] = vargpCovGrads(model, vardistx, my, indexPresent);

[gKern1, gVarmeans1, gVarcovs1] = kernVardistPsi1Gradient(model.kern, vardistx, model.X_u, gPsi1');
[gKern2, gVarmeans2, gVarcovs2] = kernVardistPsi2Gradient(model.kern, vardistx, model.X_u, gPsi2);
[gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, vardistx, gPsi0);

gVarcovs0 = (gVarcovs0(:).*vardistx.covars(:))';
gVarcovs1 = (gVarcovs1(:).*vardistx.covars(:))';
gVarcovs2 = (gVarcovs2(:).*vardistx.covars(:))';

% KL divergence terms
gVarmeansKL = - vardistx.means(:)';
% !!! the covars are optimized in the log space
gVarcovsKL = 0.5 - 0.5*vardistx.covars(:)';

gVarmeans = gVarmeans0 + gVarmeans1 + gVarmeans2 + gVarmeansKL;
gVarcovs = gVarcovs0 + gVarcovs1 + gVarcovs2 + gVarcovsKL;

g = [gVarmeans gVarcovs];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
function [gPsi0, gPsi1, gPsi2] = vargpCovGrads(model, vardistx, my, indexPresent)
%
%

d = prod(size(indexPresent));

% change the data (by including the new point and taking only the present indices)
model.m = model.m(:,indexPresent);
model.m = [my; model.m];

pointPsi1 = kernVardistPsi1Compute(model.kern, vardistx, model.X_u);
pointPsi2 = kernVardistPsi2Compute(model.kern, vardistx, model.X_u);

model.Psi1 = [pointPsi1; model.Psi1];

model.C = model.invLm * (model.Psi2 + pointPsi2) * model.invLmT;
model.TrC = sum(diag(model.C)); % Tr(C)
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C;
model.Lat = jitChol(model.At)';
model.invLat = model.Lat\eye(size(model.Lat,1));
model.invLatT = model.invLat';
model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
model.P1 = model.invLat * model.invLm; % M x M
model.P = model.P1 * (model.Psi1' * model.m);
model.TrPP = sum(sum(model.P .* model.P));
model.B = model.P1' * model.P;
P1TP1 = (model.P1' * model.P1);
Tb = (1/model.beta) * d * P1TP1;
Tb = Tb + (model.B * model.B');
model.T1 = d * model.invK_uu - Tb;

gPsi2 = (model.beta/2) * model.T1;

gPsi0 = -0.5 * model.beta * d;

gPsi1 = model.beta*(P1TP1*model.Psi1'*model.m*my');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
function [gPointDyn, model] = dynPointLogLikeGradient(model, vardistx, y, indexPresent, indexMissing)
jitter = 1e-6;

%%% RE-OPT-CODE-NEW_
if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise
    % In this case the inducing points and theta_t are optimised and
    % since updateStats is not called, any computations which need X_u
    % or theta_t must be done here.
    model = vargplvmDynamicsUpdateStats(model);
    model.K_uu = kernCompute(model.kern, model.X_u);
    model.K_uu = model.K_uu ...
        + sparseDiag(repmat(jitter, size(model.K_uu, 1), 1));
    model.Lm = jitChol(model.K_uu)';
    model.invLm = model.Lm\eye(model.k);
    model.invLmT = model.invLm';
    model.invK_uu = model.invLmT * model.invLm;
end %%% _RE-OPT-CODE-NEW

%%%%%%%%% Only the missing parts and from the training data only
mOrig = model.m;


N = model.N;
D=model.d;
model.m0=reshape(model.m,[N,D]);
model.y0=reshape(model.y,[N,D]);
Nstar = size(y,1);
Kt = zeros(N+Nstar, N+Nstar);
yMs = y(:,indexPresent);
if ~isempty(yMs)
   myMs = yMs - repmat(model.bias(indexPresent),size(yMs,1),1);  
   myMs = myMs./repmat(model.scale(indexPresent),size(yMs,1),1);  
end
% map=y;
% map(:,indexPresent)=1;%记录哪里丢失数据，有数据的为1，无数据的为0
% map(:,indexMissing)=0;%记录哪里丢失数据，有数据的为1，无数据的为0    


%model.y = model.y(:,indexMissing);
% model.m = model.m(:,indexMissing);
% model.d = prod(size(indexMissing));
% Augment the time vector to include the timestamp of the new point
model.dynamics.t = [model.dynamics.t; model.dynamics.t_star];
% Augment the reparametrized variational parameters mubar and lambda
model.dynamics.vardist.means = [model.dynamics.vardist.means; vardistx.means];
model.dynamics.vardist.covars = [model.dynamics.vardist.covars; vardistx.covars];
model.dynamics.N = model.dynamics.N + Nstar;
model.vardist.numData = model.dynamics.N;
model.dynamics.vardist.numData = model.dynamics.N;
model.vardist.nParams = 2*prod(size(model.dynamics.vardist.means));
model.dynamics.vardist.nParams = 2*prod(size(model.dynamics.vardist.means));
%model.dynamics.seq = model.dynamics.N; %%%%%% Chec
% Faster computation of the Kt matrix
% (all this could have been avoided as Kt is constant... we need to do these
% precomputations outside of this function)
% construct Kt
Kt(1:N,1:N) = model.dynamics.Kt;
% new diagonal block
Kt(N+1:end, N+1:end) = kernCompute(model.dynamics.kern, model.dynamics.t_star);
% cross block
if ~isfield(model.dynamics, 'seq') || isempty(model.dynamics.seq)
    Kt(1:N, N+1:end) = kernCompute(model.dynamics.kern, model.dynamics.t(1:N), model.dynamics.t_star);
    Kt(N+1:end, 1:N) = Kt(1:N, N+1:end)';
end
model.dynamics.Kt = Kt;
model.dynamics.fixedKt = 1;
model = vargpTimeDynamicsUpdateStats(model);
% if ~isfield(model.dynamics, 'onlyTest') || ~model.dynamics.reoptimise
    vardist2 = model.vardist;
    vardist2.means = model.vardist.means(1:end-Nstar,:);
    vardist2.covars = model.vardist.covars(1:end-Nstar,:);
    vardist2.nParams = 2*prod(size(vardist2.means));
    vardist2.numData = size(vardist2.means,1);
    model.Psi0 = kernVardistPsi0Compute(model.kern, vardist2);
    model.Psi1 = kernVardistPsi1Compute(model.kern, vardist2, model.X_u);
    model.Psi2 = kernVardistPsi2Compute(model.kern, vardist2, model.X_u);
% end

vardist2 = model.vardist;
vardist2.means = model.vardist.means(end-Nstar+1:end,:);
vardist2.covars = model.vardist.covars(end-Nstar+1:end,:);
vardist2.nParams = 2*prod(size(vardist2.means));
vardist2.numData = size(vardist2.means,1);
vardist2.map=model.map(end-Nstar+1:end,:);
pointPsi0 = kernVardistPsi0Compute(model.kern, vardist2);
pointPsi1 = kernVardistPsi1Compute(model.kern, vardist2, model.X_u);
pointPsi2 = kernVardistPsi2Compute(model.kern, vardist2, model.X_u);

model.Psi1=splice(model.Psi1,model.N,model.k,model.d,pointPsi1,Nstar);
% model.Psi1 = [model.Psi1; pointPsi1];
model.Psi2 = model.Psi2 + pointPsi2;
model.Psi0 = model.Psi0 + pointPsi0;

model.D=prod(size(model.m0))+prod(size(myMs));% 真实存在的点的 个数N*维度d
mynewMs(:,indexPresent)=myMs;
mynewMs(:,indexMissing)=0;% 缺失点对应的值为0,这样并不影响计算结果
%
model.m0=[model.m0; mynewMs];% 使用0填充的y使之与Psi1对齐，与不填充结果相同
model.m =reshape(model.m0,[],1); 
model.N = model.N + Nstar;
model.dynamics.nParams = model.dynamics.nParams + 2*prod(size(vardistx.means));
model.nParams = model.nParams + 2*prod(size(vardistx.means));

model.C = model.invLm * model.Psi2 * model.invLmT;
%model.TrC = sum(diag(model.C)); % Tr(C)
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
model.Lat = jitChol(model.At)';
model.invLat = model.Lat\eye(size(model.Lat,1));
%model.invLatT = model.invLat';
%model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
model.P1 = model.invLat * model.invLm; % M x M
model.P = model.P1 * (model.Psi1' * model.m);
%model.TrPP = sum(sum(model.P .* model.P));
model.B = model.P1' * model.P;
%model.invK_uu = model.invLmT * model.invLm;
Tb = (1/model.beta) * (model.P1' * model.P1);
Tb = Tb + (model.B * model.B');
model.T1 = model.invK_uu - Tb;

P1TP1 = (model.P1' * model.P1);
model.P = model.P1 * (model.Psi1' * model.m);
model.B = model.P1' * model.P;
Tb = (1/model.beta) * (model.P1' * model.P1);
    Tb = Tb + (model.B * model.B');
model.T1 = model.invK_uu - Tb;
myAdd=reshape(mynewMs,[],1);
gPsi1 = model.beta*(P1TP1*model.Psi1'*model.m*myAdd');
gPsi2 = (model.beta/2) * model.T1;
gPsi0 = -0.5 * model.beta;

%%单独优化测试部分时候需要
if isfield(model.dynamics,'onlyTest') && model.dynamics.onlyTest
    modelTest = model;
    %modelTest.m = [myOb; myMs]; 
    modelTest.dynamics.t = model.dynamics.t_star;
    modelTest.dynamics.vardist = vardistx; 
    modelTest.dynamics.N = size(modelTest.dynamics.vardist.means, 1);
    modelTest.vardist.numData = modelTest.dynamics.N;
    modelTest.vardist.nParams = 2*prod(size(modelTest.dynamics.vardist.means));
    Kt = kernCompute(model.dynamics.kern, model.dynamics.t_star);
    modelTest.dynamics.Kt = Kt;
    modelTest.dynamics.fixedKt = 1;
    % This enables the computation of the variational means and covars
    modelTest.dynamics.seq = []; 
    modelTest.vardist=vardist2;
%     modelTest = vargpTimeDynamicsUpdateStats(modelTest);
end
%%
if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise %%% RE-OPT-CODE-NEW
    [gPointDyn, gDynKern, gInd] = vargplvmLogLikeGradientsVar1(model,1); %%% RE-OPT-CODE-NEW
elseif isfield(model.dynamics,'onlyTest') && model.dynamics.onlyTest
    gPointDyn = vargplvmLogLikeGradientsVarOnlyTest(model,modelTest,1,gPsi0,gPsi1,gPsi2);
else %%% RE-OPT-CODE-NEW
    gPointDyn = vargplvmLogLikeGradientsVar1(model,1);
end %%% RE-OPT-CODE-NEW
if isfield(model.dynamics,'onlyTest') && model.dynamics.onlyTest
    gPointDyn = reshape(gPointDyn, modelTest.dynamics.vardist.numData, model.dynamics.q*2);
else
    gPointDyn = reshape(gPointDyn, model.dynamics.vardist.numData, model.dynamics.q*2);
end

%gPointDyn2 = gPointDyn2(N+1:end,:);

% gPointDyn = gPointDyn1 + gPointDyn2;

%%% RE-OPT-CODE-NEW_
gPointDyn = gPointDyn(:)';
if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise
%     gDynKern = gDynKern1 + gDynKern2;
%     gInd = gInd1 + gInd2;
    if strcmp(model.dynamics.kern.comp{1}.type,'rbf') || strcmp(model.dynamics.kern.comp{1}.type,'matern32')
        if ~isfield(model.dynamics, 'learnVariance') || ~model.dynamics.learnVariance
            gDynKern(2) = 0;
        end
    end
    gPointDyn = [gPointDyn 0*gDynKern gInd]; %%%%%%%%%% TEMP!!!!
end
%%% _RE-OPT-CODE-NEW



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
%function gVar = vargplvmLogLikeGradientsVar1(model, includeKL) %%% RE-OPT-CODE-REM
function [gVar, gDynKern, gInd] = vargplvmLogLikeGradientsVar1(model, includeKL) %%% RE-OPT-CODE-NEW
% Like vargplvmLogLikeGradients but only for the variational parameters (ehm, not exactly). %%% RE-OPT-CODE-MOD


% Likelihood terms (coefficients)
gPsi1 = model.beta * model.m * model.B';
gPsi1 = gPsi1'; % because it is passed to "kernVardistPsi1Gradient" as gPsi1'...
gPsi2 = (model.beta/2) * model.T1;
gPsi0 = -0.5 * model.beta;

model.vardist.map=model.map;
[gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern, model.vardist, model.X_u, gPsi1');
[gKern2, gVarmeans2, gVarcovs2, gInd2] = kernVardistPsi2Gradient(model.kern, model.vardist, model.X_u, gPsi2);
[gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, model.vardist, gPsi0);

gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2;

gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;
if includeKL
    %[gVarmeans gVarcovs gDynKern] = modelPriorReparamGrads(model.dynamics, gVarmeansLik, gVarcovsLik);
    % TODO: This may have to be changed if we have model.fixInducing==1.
    [gVarmeans gVarcovs gDynKern] = modelPriorReparamGrads(model.dynamics, gVarmeansLik, gVarcovsLik,[]);
else
    dynModel = model.dynamics;
    gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
    gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);
    gVarmeans = gVarmeansLik;
    gVarcovs = gVarcovsLik;
    
    gVarcovs = gVarcovs(:)';
    gVarmeans = gVarmeans(:)';
end

gVarcovs = (gVarcovs(:).*model.dynamics.vardist.covars(:))';
gVar = [gVarmeans gVarcovs];

%%% RE-OPT-CODE-NEW_
if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise
    % Compute Gradients with respect to X_u
    gK_uu = 0.5 * (model.T1 - (model.beta ) * model.invLmT * model.C * model.invLm);
    gKX = kernGradX(model.kern, model.X_u, model.X_u);
    gKX = gKX*2;
    dgKX = kernDiagGradX(model.kern, model.X_u);
    for i = 1:model.k
        gKX(i, :, i) = dgKX(i, :);
    end
    gX_u = zeros(model.k, model.q);
    for i = 1:model.k
        for j = 1:model.q
            gX_u(i, j) = gKX(:, j, i)'*gK_uu(:, i);
        end
    end
    gInd = gInd1 + gInd2 + gX_u(:)';
    if ~includeKL
        gDynKern = [];
    end
end
%%% _RE-OPT-CODE-NEW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
%function gVar = vargplvmLogLikeGradientsVar1(model, includeKL) %%% RE-OPT-CODE-REM
function [gVar, gDynKern, gInd] = vargplvmLogLikeGradientsVarOnlyTest(model, modelTest,includeKL,gPsi0,gPsi1,gPsi2) %%% RE-OPT-CODE-NEW
% Like vargplvmLogLikeGradients but only for the variational parameters (ehm, not exactly). %%% RE-OPT-CODE-MOD


% model.vardist.map=model.map;
[gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern, modelTest.vardist, model.X_u, gPsi1');
[gKern2, gVarmeans2, gVarcovs2, gInd2] = kernVardistPsi2Gradient(model.kern, modelTest.vardist, model.X_u, gPsi2);
[gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, modelTest.vardist, gPsi0);

gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2;

gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;
Nstar=modelTest.dynamics.N;
if includeKL
    %[gVarmeans gVarcovs gDynKern] = modelPriorReparamGrads(model.dynamics, gVarmeansLik, gVarcovsLik);
    % TODO: This may have to be changed if we have model.fixInducing==1.
    
    % means
    gVarmeans = modelTest.dynamics.Kt'*(reshape(gVarmeansLik, Nstar, model.q) - modelTest.dynamics.vardist.means); 
    % covars
    gVarcovsLik = reshape(gVarcovsLik, Nstar, model.q);
    gcovTmp = zeros(Nstar,model.dynamics.q);
    for q=1:model.dynamics.q
        LambdaH_q = modelTest.dynamics.vardist.covars(:,q).^0.5;
        Bt_q = eye(Nstar) + LambdaH_q*LambdaH_q'.*modelTest.dynamics.Kt;
        % Invert Bt_q
        Lbt_q = jitChol(Bt_q)';
        G1 = Lbt_q \ diag(LambdaH_q);
        G = G1*modelTest.dynamics.Kt;
        % Find Sq
        Sq = modelTest.dynamics.Kt - G'*G;
        gcovTmp(:,q) = - (Sq .* Sq) * (gVarcovsLik(:,q) + 0.5*modelTest.dynamics.vardist.covars(:,q));
    end
    gVarcovs = gcovTmp.*modelTest.dynamics.vardist.covars; 
    

else
    dynModel = model.dynamics;
    gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
    gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);
    gVarmeans = gVarmeansLik;
    gVarcovs = gVarcovsLik;
    
    gVarcovs = gVarcovs(:)';
    gVarmeans = gVarmeans(:)';
end

% gVarcovs = (gVarcovs(:).*model.dynamics.vardist.covars(:))';
gVar = [gVarmeans gVarcovs];

%%% RE-OPT-CODE-NEW_
if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise
    % Compute Gradients with respect to X_u
    gK_uu = 0.5 * (model.T1 - (model.beta ) * model.invLmT * model.C * model.invLm);
    gKX = kernGradX(model.kern, model.X_u, model.X_u);
    gKX = gKX*2;
    dgKX = kernDiagGradX(model.kern, model.X_u);
    for i = 1:model.k
        gKX(i, :, i) = dgKX(i, :);
    end
    gX_u = zeros(model.k, model.q);
    for i = 1:model.k
        for j = 1:model.q
            gX_u(i, j) = gKX(:, j, i)'*gK_uu(:, i);
        end
    end
    gInd = gInd1 + gInd2 + gX_u(:)';
    if ~includeKL
        gDynKern = [];
    end
end
%%% _RE-OPT-CODE-NEW
function result=splice(A,N,M,D,B,Ns)
    tt=reshape(A',[M,N,D]);
    tt=permute(tt,[2,1,3]);
    tt=reshape(tt,N,[]);
    nn=reshape(B',[M,Ns,D]);
    nn=permute(nn,[2,1,3]);
    nn=reshape(nn,Ns,[]);
    tt=[tt;nn];
    tt=reshape(tt,[N+Ns,M,D]);
    tt=permute(tt,[2,1,3]);
    result=reshape(tt,M,[])';


