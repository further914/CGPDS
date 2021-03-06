function [g, model] = varmgplvmSeqDynLogLikeGradient(model, vardistx, y)

% VARGPLVMSEQDYNLOGLIKEGRADIENT Log-likelihood gradient for of a point of the GP-LVM.
%
%	Description:
%
%	G = VARGPLVMSEQDYNLOGLIKEGRADIENT(MODEL, X, Y) returns the gradient
%	of the log likelihood with respect to the latent position, where the
%	log likelihood is conditioned on the training set.
%	 Returns:
%	  G - the gradient of the log likelihood, conditioned on the
%	   training data, with respect to the latent position.
%	 Arguments:
%	  MODEL - the model for which the gradient computation is being
%	   done.
%	  X - the latent position where the gradient is being computed.
%	  Y - the position in data space for which the computation is being
%	   done.
%	
%
%	See also
%	VARGPLVMPOINTLOGLIKELIHOOD, VARGPLVMOPTIMISEPOINT, VAGPLVMSEQUENCELOGLIKEGRADIENT

%   Combine the situation that there are observations on the first part and no observation on the second part 
%   Copyright (c) 2013 ZhaoJing 
%	Copyright (c) 2011 Michalis K. Titsias and Andreas Damianou
% 	vargplvmSeqDynLogLikeGradient.m SVN version 1312
% 	last update 2011-07-04T19:08:52.494562Z


% y is a new block/sequence 
Nstar = size(y,1);
N = model.N;
D=model.d;
mask = sum(isnan(y),2); 
indexObservedData = find(mask==0)'; 
indexMissingData = setdiff(1:Nstar, indexObservedData);

% Compute fully observed test data points and partially 
% observed data points 
yOb = y(indexObservedData, :);
yMs = y(indexMissingData, :);
mynewMs=[];
map=yMs;%记录哪里丢失数据，有数据的为1，无数据的为0
% Indices of missing dimension in the Missingdata
indexMissing = [];
indexPresent = [1:model.d];
if ~isempty(yMs)
   indexMissing = find(isnan(yMs(1,:)));
   indexPresent = setdiff(1:model.d, indexMissing);
   yMs = yMs(:,indexPresent); 
   map(:,indexPresent)=1;%记录哪里丢失数据，有数据的为1，无数据的为0
   map(:,indexMissing)=0;%记录哪里丢失数据，有数据的为1，无数据的为0    
end
    
% normalize yOb and yMs exactly as model.m is normalized 
myOb = yOb;
if ~isempty(yOb)
  myOb = yOb - repmat(model.bias,size(yOb,1),1); 
  myOb = myOb./repmat(model.scale,size(yOb,1),1); 
end
myMs = yMs;
if ~isempty(yMs)
   myMs = yMs - repmat(model.bias(indexPresent),size(yMs,1),1);  
   myMs = myMs./repmat(model.scale(indexPresent),size(yMs,1),1);  
end
if ~isfield(model,'m0')
model.m0=reshape(model.m,[N,D]);
model.y0=reshape(model.y,[N,D]);
end
mOrig = model.m0;
%modelOrig = model;

% re-order test data so that observed are first and then are the missing 
Order = [indexObservedData, indexMissingData];
model.dynamics.t_star = model.dynamics.t_star(Order);
vardistx.means = vardistx.means(Order,:);
vardistx.covars = vardistx.covars(Order,:);
nObsData = size(indexObservedData,2);
nMissData=size(indexMissingData,2);

% Form the modelTest for the new block which allows to compute the variational
% distribution (given the barmus and lambdas) for this new block. Notice that this distribution
% is independent from the other blocks in the training set. 
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
modelTest = vargpTimeDynamicsUpdateStats(modelTest);


% The fully observed subset of data from the new block must be augmented with all previous
% training data to form the LL1 term in the whole bound 
gVarmeansLik = zeros(1, nObsData*model.dynamics.q);
gVarcovsLik = zeros(1, nObsData*model.dynamics.q);
vardistOb.map=[];
vardistMs.map=[];
if ~isempty(yOb)
   vardistOb = modelTest.vardist; 
   vardistOb.means = modelTest.vardist.means(1:nObsData, :);
   vardistOb.covars = modelTest.vardist.covars(1:nObsData, :);
   vardistOb.nParams = 2*prod(size(vardistOb.means));
   vardistOb.numData = size(vardistOb.means,1);
   vardistOb.map=ones(nObsData,model.d);
   % Psi statistics for the data of the new block/sequence which have a fully
   % observed features/dimensions
   obsPsi0 = kernVardistPsi0Compute(model.kern, vardistOb);
   obsPsi1 = kernVardistPsi1Compute(model.kern, vardistOb, model.X_u);
   obsPsi2 = kernVardistPsi2Compute(model.kern, vardistOb, model.X_u);
   model.Psi1=splice(model.Psi1,N,model.k,model.d,obsPsi1,nObsData);
%    model.Psi1 = [model.Psi1; obsPsi1]; 
   model.Psi2 = model.Psi2 + obsPsi2;
   model.Psi0 = model.Psi0 + obsPsi0;
    model.m0=[mOrig; myOb];
    model.D=prod(size(model.m0));
    model.m =reshape(model.m0,[],1); 
    model.N =  N + nObsData;
    model.d = size(model.m0,2);

%     model.TrYY = sum(sum(model.m .* model.m));
%     model.C = model.invLm * model.Psi2 * model.invLmT;
%     model.TrC = sum(diag(model.C)); % Tr(C)
%     model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
%     model.Lat = jitChol(model.At)';
%     model.invLat = model.Lat\eye(size(model.Lat,1));  
%     model.invLatT = model.invLat';
%     model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
%     model.P1 = model.invLat * model.invLm; % M x M
%     model.P = model.P1 * (model.Psi1' * model.m);
%     model.TrPP = sum(sum(model.P .* model.P));
end
if ~isempty(indexMissing)% to make computing more efficient
    vardistMs = modelTest.vardist; 
    vardistMs.means = modelTest.vardist.means(nObsData+1:end, :);
    vardistMs.covars = modelTest.vardist.covars(nObsData+1:end, :);
    vardistMs.nParams = 2*prod(size(vardistMs.means));
    vardistMs.numData = size(vardistMs.means,1);
    % Psi statistics for the data of the new block/sequence which have partially
    % observed dimensions
    vardistMs.map=map;% 缺失点的映射
    missingPsi0 = kernVardistPsi0Compute(model.kern, vardistMs);
    missingPsi1 = kernVardistPsi1Compute(model.kern, vardistMs, model.X_u);
    missingPsi2 = kernVardistPsi2Compute(model.kern, vardistMs, model.X_u);  

    model.Psi1=splice(model.Psi1,model.N,model.k,model.d,missingPsi1,nMissData);
%     model.Psi1 = [model.Psi1; missingPsi1]; 
    model.Psi2 = model.Psi2 + missingPsi2;
    model.Psi0 = model.Psi0 + missingPsi0;
    model.D=prod(size(model.m0))+prod(size(myMs));% 真实存在的点的 个数N*维度d
    mynewMs(:,indexPresent)=myMs;
    mynewMs(:,indexMissing)=0;% 缺失点对应的值为0,这样并不影响计算结果
    model.m0=[model.m0; mynewMs];% 使用0填充的y使之与Psi1对齐，与不填充结果相同
    model.m =reshape(model.m0,[],1); 
    model.N =  model.N + nMissData;
    
    model.d = size(model.m0,2);
    
%    model.C = model.invLm * model.Psi2 * model.invLmT;
%    model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
%    model.Lat = jitChol(model.At)';
%    model.invLat = model.Lat\eye(size(model.Lat,1));  
%    model.P1 = model.invLat * model.invLm; % M x M
%    P1TP1 = (model.P1' * model.P1);
%    model.P = model.P1 * (model.Psi1' * model.m);
%    model.B = model.P1' * model.P;
%    Tb = (1/model.beta) * model.d * (model.P1' * model.P1);
%         Tb = Tb + (model.B * model.B');
%    model.T1 = model.invK_uu - Tb;

   % Precompuations for the gradients
  
end
model.C = model.invLm * model.Psi2 * model.invLmT;
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
model.Lat = jitChol(model.At)';
model.invLat = model.Lat\eye(size(model.Lat,1));  
model.P1 = model.invLat * model.invLm; % M x M
P1TP1 = (model.P1' * model.P1);
model.P = model.P1 * (model.Psi1' * model.m);
model.B = model.P1' * model.P;
Tb = (1/model.beta) * (model.P1' * model.P1);
    Tb = Tb + (model.B * model.B');
model.T1 = model.invK_uu - Tb;
myAdd=reshape([myOb;mynewMs],[],1);
gPsi1 = model.beta*(P1TP1*model.Psi1'*model.m*myAdd');
%gPsi1 = model.beta * model.m(:,indexMissing) * model.B';
%gPsi1 = gPsi1'; % because it is passed to "kernVardistPsi1Gradient" as gPsi1'...
gPsi2 = (model.beta/2) * model.T1;
gPsi0 = -0.5 * model.beta;
modelTest.vardist.map=[vardistOb.map;vardistMs.map];
[gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern, modelTest.vardist, model.X_u, gPsi1');
[gKern2, gVarmeans2, gVarcovs2, gInd2] = kernVardistPsi2Gradient(model.kern, modelTest.vardist, model.X_u, gPsi2);
[gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, modelTest.vardist, gPsi0);
gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2;
gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;   



 
%

gPointDyn = zeros(Nstar, model.dynamics.q*2); 
% means
gPointDyn(:,1:model.dynamics.q) = modelTest.dynamics.Kt'*(reshape(gVarmeansLik, Nstar, model.q) - modelTest.dynamics.vardist.means); 
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
gPointDyn(:,(model.dynamics.q+1):(model.dynamics.q*2)) = gcovTmp.*modelTest.dynamics.vardist.covars; 


% applying the inverse Order to give the gradeitn in the 
%orginal order of the data
T = eye(Nstar);
T = T(Order,:);
% inverse permutation
T = T';
% there must be a better way to do this
InvOrder = zeros(1,Nstar);
for i=1:Nstar
    InvOrder(i) = find(T(i,:)==1); 
end
gPointDyn = gPointDyn(InvOrder, :);

g = gPointDyn(:)';
end
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
end