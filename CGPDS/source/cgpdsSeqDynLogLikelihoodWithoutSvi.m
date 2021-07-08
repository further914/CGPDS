function [ll, model] = cgpdsSeqDynLogLikelihoodWithoutSvi(model, vardistx, y,samd,samr)

% cgpdsSeqDynLogLikelihoodWithoutSvi Log-likelihood of a point for CGPDS

% !!!!!! this function can become faster with precomputations stored in the
% structure model !!!!! 


% y is a new block/sequence 
Nstar = size(y,1);
N = model.N;

mask = sum(isnan(y),2); 
indexObservedData = find(mask==0)'; 
indexMissingData = setdiff(1:Nstar, indexObservedData);

% Compute fully observed test data points and partially 
% observed data points 
yOb = y(indexObservedData, :);%0*59
yMs = y(indexMissingData, :);%43*59

% Indices of missing dimension in the Missingdata
indexMissing = [];
indexPresent = [1:model.d];
if ~isempty(yMs)
   indexMissing = find(isnan(yMs(1,:)));
   indexPresent = setdiff(1:model.d, indexMissing);
   yMs = yMs(:,indexPresent);   
end
    
% normalize yOb and yMs exactly as model.m is normalized 
myOb = yOb;
if ~isempty(yOb)
  myOb = yOb - repmat(model.bias,size(yOb,1),1); 
  myOb = myOb./repmat(model.scale,size(yOb,1),1); 
end


myMs = yMs;%43*52
if ~isempty(yMs)
   myMs = yMs - repmat(model.bias(indexPresent),size(yMs,1),1);  
   myMs = myMs./repmat(model.scale(indexPresent),size(yMs,1),1);  
end
mOrig = model.m;%353*59

% re-order test data so that observed are first and then are the missing 
Order = [indexObservedData, indexMissingData];
model.dynamics.t_star = model.dynamics.t_star(Order);
vardistx.means = vardistx.means(Order,:);
vardistx.covars = vardistx.covars(Order,:);
nObsData = size(indexObservedData,2);


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


% The fully observed subset of data from the new block must be augmented with all previous
% training data to form the LL1 term in the whole bound 
% if ~isempty(yOb)
%    vardistOb = modelTest.vardist; 
%    vardistOb.means = modelTest.vardist.means(1:nObsData, :);
%    vardistOb.covars = modelTest.vardist.covars(1:nObsData, :);
%    vardistOb.nParams = 2*prod(size(vardistOb.means));
%    vardistOb.numData = size(vardistOb.means,1);
%  
%    % Psi statistics for the data of the new block/sequence which have a fully
%    % observed features/dimensions
%    obsPsi0 = kernVardistPsi1Compute(model.kern_v, vardistOb, model.X_v);
%    obsPsi1 = kernVardistPsi1Compute(model.kern_u, vardistOb, model.X_u);
%    
%    obsPsi2 = kernVardistPsi0Compute(model.kern_v, vardistOb);
%    obsPsi3 = kernVardistPsi0Compute(model.kern_u, vardistOb);
%    
%    obsPsi4 = kernVardistPsi2Compute(model.kern_v, vardistOb, model.X_v);
%    obsPsi5 = kernVardistPsi2Compute(model.kern_u, vardistOb, model.X_u);
%    
%    model.Psi0 = [model.Psi0; obsPsi0];
%    model.Psi1 = [model.Psi1; obsPsi1];
%    
%    model.Psi2 = model.Psi2 + obsPsi2;
%    model.Psi3 = model.Psi3 + obsPsi3;
%    
%    model.Psi4 = model.Psi4 + obsPsi4;
%    model.Psi5 = model.Psi5 + obsPsi5;
%    
%    if ~isempty(indexMissing)
%         model.N = model.N + size(yOb,1);
%         % Is model.m= [mOrig(:,indexMissing); myOb(:, indexMissing)] also OK? Then model.m(:,indexMissing) can be replaced by model.m
%         model.m = [model.m; myOb];
%        
%         model.C_v = model.invLm_v * model.Psi4 * model.invLmT_v;
%         model.TrC_v = sum(diag(model.C_v)); % Tr(C)
%         model.At_v = (1/model.beta) * eye(size(model.C_v,1)) + model.C_v; % At = beta^{-1} I + C
%         model.Lat_v = jitChol(model.At_v)';
%         model.invLat_v = model.Lat_v\eye(size(model.Lat_v,1));  
%         model.invLatT_v = model.invLat_v';
%         model.logDetAt_v = 2*(sum(log(diag(model.Lat_v)))); % log |At|
%         model.P1_v = model.invLat_v * model.invLm_v; % M x M
%         % once in the calculations (P1: MxM, Psi1':MxN, Y: NxD)
%         model.P_v = model.P1_v * (model.Psi0' * model.m);
%         model.TrPP_v = sum(sum(model.P_v .* model.P_v));
%        
%         lend = length(indexPresent);
%         model.C_u = cell(lend,model.J);
%         model.TrC_u = zeros(lend,model.J);
%         model.At_u = cell(lend,model.J);
%         model.Lat_u = cell(lend,model.J);
%         model.invLat_u = cell(lend,model.J);
%         model.invLatT_u = cell(lend,model.J);
%         model.logDetAt_u = zeros(lend,model.J);
%         model.P1_u = cell(lend,model.J);
%         model.P_u = cell(lend,model.J);
%         model.TrPP_u = 0;
%         model.B_u = cell(lend,model.J);
%         model.T1_u = zeros(model.k,model.k);
% 
%         model.C_u_temp = model.invLm_u * model.Psi5 * model.invLmT_u;  
% 
%         for d = 1:length(indexPresent)
%             index = indexPresent(d);
%             for j = 1:model.J
% 
%                 model.C_u{d,j} = model.W(index,j)^2 * model.C_u_temp;
%                 model.TrC_u(d,j) = sum(diag(model.C_u{d,j})); % Tr(C)
%                 % Matrix At replaces the matrix A of the old implementation; At is more stable
%                 % since it has a much smaller condition number than A=sigma^2 K_uu + Psi2
%                 model.At_u{d,j} = (1/model.beta) * eye(size(model.C_u{d,j},1)) + model.C_u{d,j}; % At = beta^{-1} I + C
%                 model.Lat_u{d,j} = jitChol(model.At_u{d,j})';%lower bound
%                 model.invLat_u{d,j} = model.Lat_u{d,j}\eye(size(model.Lat_u{d,j},1));  
%                 model.invLatT_u{d,j} = model.invLat_u{d,j}';
%                 model.logDetAt_u(d,j) = 2*(sum(log(diag(model.Lat_u{d,j})))); % log |At|
% 
%                 model.P1_u{d,j} = model.invLat_u{d,j} * model.invLm_u; % M x M
% 
%                 % First multiply the two last factors; so, the large N is only involved
%                 % once in the calculations (P1: MxM, Psi1':MxN, Y: NxD)
%                 model.P_u{d,j} = model.P1_u{d,j} * (model.Psi1' * model.m(:,d));%M*1  
% 
%                 % Needed for both, the bound's and the derivs. calculations.
%         %         model.TrPP_u = sum(sum(model.P_u .* model.P_u));
%                 model.TrPP_u = model.TrPP_u + model.W(index,j)^2*model.P_u{d,j}'*model.P_u{d,j};
%                 model.B_u{d,j} = model.P1_u{d,j}' * model.P_u{d,j}; %M*1
% 
%             end
%         end
%        
%    end
% else
%    P_v = model.P1_v * (model.Psi0' * model.m(:,indexMissing));%30*7
%    TrPP_v = sum(sum(P_v .* P_v));
%    
%    TrPP_u = 0;
%    for i = 1:length(indexMissing)
%        index = indexMissing(i);
%        for j = 1 : model.J
%            P_u = model.P1_u{index,j} * (model.Psi1' * model.m(:,index));%30*1
%            TrPP_u = TrPP_u + model.W(index,j)^2 * P_u' * P_u;
%        end
%    end
   
% end 



% LL1 TERM
%
ll1 = 0;
if ~isempty(indexMissing)
   dmis = prod(size(indexMissing));
   
   % Precompute again the parts that contain Y
   TrYY = sum(sum(model.m(:,indexMissing) .* model.m(:,indexMissing)));
   
   ll1 = -0.5*(dmis*(-(model.N-model.k)*log(model.beta) ...
				  + model.logDetAt_v) ...
	      - (TrPP_v ...
	      - TrYY)*model.beta);
   ll1 = ll1 - 0.5*model.beta*dmis*model.Psi2 + 0.5*dmis*model.beta*model.TrC_v;
   ll1 = ll1-dmis*model.N/2*log(2*pi);
   
   for i = 1:length(indexMissing)
       index = indexMissing(i);
       for j = 1:model.J
                ll1 = ll1 - 0.5*(model.logDetAt_u(index,j) - (-model.k)*log(model.beta));
%                 %if strcmp(model.approx, 'dtcvar')
                ll1 = ll1 - 0.5*model.beta*model.W(index,j)^2*model.Psi3 + 0.5*model.beta*model.TrC_u(index,j);
       end
   end
   ll1 = ll1 + 0.5*model.beta*TrPP_u;
   
end


% The data in the test block/seq with missing values are processed to get their
% psi statisitcs. This is needed to form the LL2 term that is the part 
% of the bound corresponding to dimensions observed everywhere (in all
% training and test points)
if ~isempty(indexMissing)
    vardistMs = modelTest.vardist; 
    vardistMs.means = modelTest.vardist.means(nObsData+1:end, :);
    vardistMs.covars = modelTest.vardist.covars(nObsData+1:end, :);
    vardistMs.nParams = 2*prod(size(vardistMs.means));
    vardistMs.numData = size(vardistMs.means,1);
  
    % Psi statistics for the data of the new block/sequence which have partially
    % observed dimensions
    missingPsi0 = kernVardistPsi1Compute(model.kern_v, vardistMs, model.X_v);
    missingPsi1 = kernVardistPsi1Compute(model.kern_u, vardistMs, model.X_u);
    
    missingPsi2 = kernVardistPsi0Compute(model.kern_v, vardistMs);
    missingPsi3 = kernVardistPsi0Compute(model.kern_u, vardistMs);
    
    missingPsi4 = kernVardistPsi2Compute(model.kern_v, vardistMs, model.X_v);
    missingPsi5 = kernVardistPsi2Compute(model.kern_u, vardistMs, model.X_u);
    
    model.Psi0 = [model.Psi0; missingPsi0];
    model.Psi1 = [model.Psi1; missingPsi1];  
    
    model.Psi2 = model.Psi2 + missingPsi2;
    model.Psi3 = model.Psi3 + missingPsi3;
    
    model.Psi4 = model.Psi4 + missingPsi4;
    model.Psi5 = model.Psi5 + missingPsi5;
end



model.m = [mOrig(:,indexPresent); myOb(:, indexPresent); myMs]; %(353+43)*52
model.N =  N + Nstar;
model.d = prod(size(indexPresent));

model.TrYY = sum(sum(model.m .* model.m));
model.C_v = model.invLm_v * model.Psi4 * model.invLmT_v;
model.TrC_v = sum(diag(model.C_v)); % Tr(C)
model.At_v = (1/model.beta) * eye(size(model.C_v,1)) + model.C_v; % At = beta^{-1} I + C
model.Lat_v = jitChol(model.At_v)';
model.invLat_v = model.Lat_v\eye(size(model.Lat_v,1));  
model.invLatT_v = model.invLat_v';
model.logDetAt_v = 2*(sum(log(diag(model.Lat_v)))); % log |At|
model.P1_v = model.invLat_v * model.invLm_v; % M x M
% once in the calculations (P1: MxM, Psi1':MxN, Y: NxD)
model.P_v = model.P1_v * (model.Psi0' * model.m);
model.TrPP_v = sum(sum(model.P_v .* model.P_v));

fprintf('model.TrPP_v is %f\n',model.TrPP_v);
    

model.C_u = cell(model.d,model.J);
model.TrC_u = zeros(model.d,model.J);
model.At_u = cell(model.d,model.J);
model.Lat_u = cell(model.d,model.J);
model.invLat_u = cell(model.d,model.J);
model.invLatT_u = cell(model.d,model.J);
model.logDetAt_u = zeros(model.d,model.J);
model.P1_u = cell(model.d,model.J);
model.P_u = cell(model.d,model.J);
model.TrPP_u = 0;
model.B_u = cell(model.d,model.J);
model.T1_u = zeros(model.k,model.k);

model.C_u_temp = model.invLm_u * model.Psi5 * model.invLmT_u;  

for d = 1:length(indexPresent)
    index = indexPresent(d);
    for j = 1:model.J
        
        model.C_u{d,j} = model.W(index,j)^2 * model.C_u_temp;
        model.TrC_u(d,j) = sum(diag(model.C_u{d,j})); % Tr(C)
        % Matrix At replaces the matrix A of the old implementation; At is more stable
        % since it has a much smaller condition number than A=sigma^2 K_uu + Psi2
        model.At_u{d,j} = (1/model.beta) * eye(size(model.C_u{d,j},1)) + model.C_u{d,j}; % At = beta^{-1} I + C
        model.Lat_u{d,j} = jitChol(model.At_u{d,j})';%lower bound
        model.invLat_u{d,j} = model.Lat_u{d,j}\eye(size(model.Lat_u{d,j},1));  
        model.invLatT_u{d,j} = model.invLat_u{d,j}';
        model.logDetAt_u(d,j) = 2*(sum(log(diag(model.Lat_u{d,j})))); % log |At|

        model.P1_u{d,j} = model.invLat_u{d,j} * model.invLm_u; % M x M

        % First multiply the two last factors; so, the large N is only involved
        % once in the calculations (P1: MxM, Psi1':MxN, Y: NxD)
        model.P_u{d,j} = model.P1_u{d,j} * (model.Psi1' * model.m(:,d));%M*1  

        % Needed for both, the bound's and the derivs. calculations.
%         model.TrPP_u = sum(sum(model.P_u .* model.P_u));
        model.TrPP_u = model.TrPP_u + model.W(index,j)^2*model.P_u{d,j}'*model.P_u{d,j};
        model.B_u{d,j} = model.P1_u{d,j}' * model.P_u{d,j}; %M*1
     
    end
end
% fprintf('model.TrPP_u is %f\n',model.TrPP_u);

% model.TrPP_u1=0;
% for d = 1:length(indexPresent)
%     index = indexPresent(d);
%     for j = 1:model.J
%         model.TrPP_u1 = model.TrPP_u1 + trace(model.W(index,j)^2*model.m(:,d)'*model.Psi1*inv(model.W(index,j)^2*model.Psi5+model.K_uu/model.beta)...
%             *model.Psi1'*model.m(:,d));
%     end
% end
% 
% fprintf('model.TrPP_u is %f\n',model.TrPP_u1);


% LL2 TERM 
%
ll2 = 0;
if ~isempty(indexPresent)   
     
    ll2 = -0.5*(model.d*(-(model.N-model.k)*log(model.beta) ...
				  + model.logDetAt_v) ...
	      - (model.TrPP_v ...
	      - model.TrYY)*model.beta);
    ll2 = ll2 - 0.5*model.beta*model.d*model.Psi2 + 0.5*model.d*model.beta*model.TrC_v;

    ll2 = ll2-model.d*model.N/2*log(2*pi);
    
    for d = 1:length(indexPresent)
       index = indexPresent(d);
       for j = 1:model.J
                ll2 = ll2 - 0.5*(model.logDetAt_u(d,j) - (-model.k)*log(model.beta));
%                 %if strcmp(model.approx, 'dtcvar')
                ll2 = ll2 - 0.5*model.beta*model.W(index,j)^2*model.Psi3 + 0.5*model.beta*model.TrC_u(d,j);
       end
    end
    ll2 = ll2 + 0.5*model.beta*model.TrPP_u;
    
%    
end


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
ll = ll1 + ll2 + KLdiv; 
