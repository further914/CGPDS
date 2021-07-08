function [ll, model] = varmgplvmSeqDynLogLikelihood(model, vardistx, y)

% VARGPLVMSEQDYNLOGLIKELIHOOD Log-likelihood of a point for the GP-LVM.
%
%	Description:
%
%	VARGPLVMSEQDYNLOGLIKELIHOOD(MODEL, VARDISTX, Y, INDEXPRESENT)
%	returns the log likelihood of a latent point and an observed data
%	point for the posterior prediction of the GP-LVM model.
%	 Arguments:
%	  MODEL - the model for which the point prediction will be made.
%	  VARDISTX - the variational distribution over latent point for
%	   which the posterior distribution will be evaluated. It contains
%	   the mean and the diagonal covarriance
%	  Y - the observed data point for which the posterior is evaluated
%	  INDEXPRESENT - indicates which indices from the observed vector
%	   are present (e.g. when when all indices are present, then y will
%	   d-dimensional and indexPresent = 1:D)
%	
%
%	See also
%	VARGPLVMCREATE, VARGPLVMOPTIMISEPOINT, VARGPLVMPOINTOBJECTIVE


%	Copyright (c) 2011 Michalis K. Titsias and Andreas Damianou
% 	vargplvmSeqDynLogLikelihood.m SVN version 1312
% 	last update 2011-04-18T11:54:45.000000Z


% !!!!!! this function can become faster with precomputations stored in the
% structure model !!!!! 


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
modelTest.dynamics.N =  size(modelTest.dynamics.vardist.means, 1);
modelTest.vardist.numData = modelTest.dynamics.N;
modelTest.vardist.nParams = 2*prod(size(modelTest.dynamics.vardist.means));
Kt = kernCompute(model.dynamics.kern, model.dynamics.t_star);
modelTest.dynamics.Kt = Kt;
modelTest.dynamics.fixedKt = 1;
% This enables the computation of the variational means and covars
modelTest.dynamics.seq = []; 
modelTest = vargpTimeDynamicsUpdateStats(modelTest);

ll=0;
% The fully observed subset of data from the new block must be augmented with all previous
% training data to form the LL1 term in the whole bound 
if ~isempty(yOb) % 计算全部都有观测的点的似然
   vardistOb = modelTest.vardist; 
   vardistOb.means = modelTest.vardist.means(1:nObsData, :);
   vardistOb.covars = modelTest.vardist.covars(1:nObsData, :);
   vardistOb.nParams = 2*prod(size(vardistOb.means));
   vardistOb.numData = size(vardistOb.means,1);
 
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
    model.N =  model.N + nObsData;
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
%     
    
    
end 





% The data in the test block/seq with missing values are processed to get their
% psi statisitcs. This is needed to form the LL2 term that is the part 
% of the bound corresponding to dimensions observed everywhere (in all
% training and test points)

if ~isempty(indexMissing)% 计算存在缺失点的似然
    vardistMs = modelTest.vardist; 
    vardistMs.means = modelTest.vardist.means(nObsData+1:end, :);
    vardistMs.covars = modelTest.vardist.covars(nObsData+1:end, :);
    vardistMs.nParams = 2*prod(size(vardistMs.means));
    vardistMs.numData = size(vardistMs.means,1);
    vardistMs.map=map;% 缺失点的映射
    % Psi statistics for the data of the new block/sequence which have partially
    % observed dimensions
    
    missingPsi0 = kernVardistPsi0Compute(model.kern, vardistMs);
    missingPsi1 = kernVardistPsi1Compute(model.kern, vardistMs, model.X_u);% 缺失点对应的值为0
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

    
end
model.TrYY = sum(sum(model.m .* model.m));
model.C = model.invLm * model.Psi2 * model.invLmT;
model.TrC = sum(diag(model.C)); % Tr(C)
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
model.Lat = jitChol(model.At)';
model.invLat = model.Lat\eye(size(model.Lat,1));  
model.invLatT = model.invLat';
model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
model.P1 = model.invLat * model.invLm; % M x M
model.P = model.P1 * (model.Psi1' * model.m);
model.TrPP = sum(sum(model.P .* model.P));


ll = -0.5*(-(model.D-model.k)*log(model.beta) ...
				  + model.logDetAt ...
	      - (model.TrPP ...
	      - model.TrYY)*model.beta);
ll = ll - 0.5*model.beta*model.Psi0 + 0.5*model.beta*model.TrC;

ll = ll-model.D/2*log(2*pi);
% KL TERM 
%
model.dynamics.t = [model.dynamics.t; model.dynamics.t_star];
% Augment the reparametrized variational parameters mubar and lambda
model.dynamics.vardist.means = [model.dynamics.vardist.means; vardistx.means];
model.dynamics.vardist.covars = [model.dynamics.vardist.covars; vardistx.covars]; 
model.dynamics.N = N + Nstar;
model.vardist.numData = model.dynamics.N;
Kt = zeros(N+Nstar,N+Nstar);
Kt(1:N,1:N) = model.dynamics.Kt; 
Kt(N+1:end, N+1:end) = modelTest.dynamics.Kt; 
model.dynamics.Kt = Kt;
model.dynamics.fixedKt = 1;
model.dynamics.seq = [model.dynamics.seq, (N+Nstar)];
KLdiv = modelVarPriorBound(model);

% sum all terms
ll = ll + KLdiv; 
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