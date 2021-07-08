function [g, model] = cgpdsSeqDynLogLikeGradient(model, vardistx, y,samd,samr)

% cgpdsSeqDynLogLikeGradient Log-likelihood gradient of CGPDS for testing data. 

%   Combine the situation that there are observations on the first part and no observation on the second part 
%   Copyright (c) 2013 ZhaoJing 
%	Copyright (c) 2011 Michalis K. Titsias and Andreas Damianou
% 	vargplvmSeqDynLogLikeGradient.m SVN version 1312
% 	last update 2011-07-04T19:08:52.494562Z


% y is a new block/sequence 
Nstar = size(y,1);
N = model.N;

mask = sum(isnan(y),2); 
indexObservedData = find(mask==0)'; %%%empty
indexMissingData = setdiff(1:Nstar, indexObservedData);

% Compute fully observed test data points and partially 
% observed data points 
yOb = y(indexObservedData, :);%%empty
yMs = y(indexMissingData, :);

% Indices of missing dimension in the Missingdata
indexMissing = [];
indexPresent = [1:model.d];
if ~isempty(yMs)
   indexMissing = find(isnan(yMs(1,:)));
   indexPresent = setdiff(1:model.d, indexMissing);
%    yMs = yMs(:,indexPresent);   
end
    
% normalize yOb and yMs exactly as model.partm is normalized 
myOb = yOb;
if ~isempty(yOb)
  myOb = yOb - repmat(model.bias,size(yOb,1),1); 
  myOb = myOb./repmat(model.scale,size(yOb,1),1); 
end


myMs = yMs(:,indexPresent);
if ~isempty(myMs)
   myMs = myMs - repmat(model.bias(indexPresent),size(myMs,1),1);  
   myMs = myMs./repmat(model.scale(indexPresent),size(myMs,1),1);  
end
mOrig = model.m;
%modelOrig = model;

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
% if ~isempty(yOb)
%    vardistOb = modelTest.vardist; 
%    vardistOb.means = modelTest.vardist.means(1:nObsData, :);
%    vardistOb.covars = modelTest.vardist.covars(1:nObsData, :);
%    vardistOb.nParams = 2*prod(size(vardistOb.means));
%    vardistOb.numData = size(vardistOb.means,1);
%  
%    % Psi statistics for the data of the new block/sequence which have a fully
%    % observed features/dimensions
%    obsPsi0 = kernVardistPsi0Compute(model.kern, vardistOb);
%    obsPsi1 = kernVardistPsi1Compute(model.kern, vardistOb, model.X_u);
%    obsPsi2 = kernVardistPsi2Compute(model.kern, vardistOb, model.X_u);
%    model.Psi1 = [model.Psi1; obsPsi1]; 
%    model.Psi2 = model.Psi2 + obsPsi2;
%    model.Psi0 = model.Psi0 + obsPsi0;
%    if ~isempty(indexMissing)% to make computing more efficient
%        model.N = N + size(yOb,1);
%        model.partm = [model.partm; myOb];
%        model.d = prod(size(indexMissing));
% 
%        model.C = model.invLm * model.Psi2 * model.invLmT;
%        model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
%        model.Lat = jitChol(model.At)';
%        model.invLat = model.Lat\eye(size(model.Lat,1));  
%        model.P1 = model.invLat * model.invLm; % M x M
%        P1TP1 = (model.P1' * model.P1);
%        model.P = model.P1 * (model.Psi1' * model.partm(:,indexMissing));
%        model.B = model.P1' * model.P;
%        Tb = (1/model.beta) * model.d * (model.P1' * model.P1);
%             Tb = Tb + (model.B * model.B');
%        model.T1 = model.d * model.invK_uu - Tb;
% 
%        % Precompuations for the gradients
%        gPsi1 = model.beta*(P1TP1*model.Psi1'*model.partm(:,indexMissing)*myOb(:,indexMissing)');
%        %gPsi1 = model.beta * model.partm(:,indexMissing) * model.B';
%        %gPsi1 = gPsi1'; % because it is passed to "kernVardistPsi1Gradient" as gPsi1'...
%        gPsi2 = (model.beta/2) * model.T1;
%        gPsi0 = -0.5 * model.beta * model.d;
%        [gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern, vardistOb, model.X_u, gPsi1');
%        [gKern2, gVarmeans2, gVarcovs2, gInd2] = kernVardistPsi2Gradient(model.kern, vardistOb, model.X_u, gPsi2);
%        [gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, vardistOb, gPsi0);
%        gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2;
%        gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;   
%    end
% end 


% GRADIENT FOR THE LL1 TERM
%
gPointDyn1 = zeros(Nstar, model.dynamics.q*2);%43*18
if ~isempty(indexMissing)
%    
   % means
   gPointDyn1(:,1:model.dynamics.q) = modelTest.dynamics.Kt(1:nObsData,:)'*reshape(gVarmeansLik, nObsData, model.q);
   % covars
   gVarcovsLik = reshape(gVarcovsLik, nObsData, model.q);
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
       Sq = - (Sq .* Sq);
       % only the cross matrix is needed as in barmu case
       gcovTmp(:,q) = Sq(1:nObsData,:)'*gVarcovsLik(:,q);
   end
   gPointDyn1(:,(model.dynamics.q+1):(model.dynamics.q*2)) = gcovTmp.*modelTest.dynamics.vardist.covars;% zero matrix????
%   
end


% GRADIENT FOR THE LL2 TERM PLUS KL TERM
%
% The data in the test block/seq with missing values are processed to get their
% psi statisitcs. This is needed to form the LL2 term that is the part 
% of the bound corresponding to dimensions observed everywhere (in all
% training and test points)
if ~isempty(indexMissing)% to make computing more efficient
    vardistMs = modelTest.vardist;%43*9 
    vardistMs.means = modelTest.vardist.means(nObsData+1:end, :);
    vardistMs.covars = modelTest.vardist.covars(nObsData+1:end, :);
    vardistMs.nParams = 2*prod(size(vardistMs.means));
    vardistMs.numData = size(vardistMs.means,1);
    % Psi statistics for the data of the new block/sequence which have partially


    missingPsi0 = kernVardistPsi1Compute(model.kern_v, vardistMs, model.X_v);%43*30
    missingPsi1 = kernVardistPsi1Compute(model.kern_u, vardistMs, model.X_u);
    
    missingPsi2 = kernVardistPsi0Compute(model.kern_v, vardistMs);
    missingPsi3 = kernVardistPsi0Compute(model.kern_u, vardistMs);
    
    missingPsi4 = kernVardistPsi2Compute(model.kern_v, vardistMs, model.X_v);%30*30
    missingPsi5 = kernVardistPsi2Compute(model.kern_u, vardistMs, model.X_u);
    
    model.Psi0 = [model.Psi0; missingPsi0];
    model.Psi1 = [model.Psi1; missingPsi1];  
    
    model.Psi2 = model.Psi2 + missingPsi2;
    model.Psi3 = model.Psi3 + missingPsi3;
    
    model.Psi4 = model.Psi4 + missingPsi4;
    model.Psi5 = model.Psi5 + missingPsi5;
end

model.N = N + Nstar;
model.m = [mOrig(:,indexPresent); myOb(:, indexPresent); myMs]; %(353+43)*52
model.d = prod(size(indexPresent));
model.dynamics.nParams = model.dynamics.nParams + 2*prod(size(vardistx.means));
model.nParams = model.nParams + 2*prod(size(vardistx.means));



% model.T1_v = - model.B_v * model.B_v';

% lenIndexPresent = length(indexPresent);
% model.C_u = cell(lenIndexPresent,model.J);
% model.TrC_u = zeros(lenIndexPresent,model.J);
% model.At_u = cell(lenIndexPresent,model.J);
% model.Lat_u = cell(lenIndexPresent,model.J);
% model.invLat_u = cell(lenIndexPresent,model.J);
% model.invLatT_u = cell(lenIndexPresent,model.J);
% model.logDetAt_u = zeros(lenIndexPresent,model.J);
% model.P1_u = cell(lenIndexPresent,model.J);
% model.P_u = cell(lenIndexPresent,model.J);
% model.TrPP_u = 0;
% model.B_u = cell(lenIndexPresent,model.J);
% model.T1_u = zeros(model.k,model.k);
% P1TP1_u = cell(lenIndexPresent,model.J);
% 
% for d = 1:length(indexPresent)
%     index = indexPresent(d);
%     for j = 1:model.J
%         
%         model.C_u{d,j} = model.W(index,j)^2 * model.C_u_temp;
%         model.TrC_u(d,j) = sum(diag(model.C_u{d,j})); % Tr(C)
%         % Matrix At replaces the matrix A of the old implementation; At is more stable
%         % since it has a much smaller condition number than A=sigma^2 K_uu + Psi2
%         model.At_u{d,j} = (1/model.beta) * eye(size(model.C_u{d,j},1)) + model.C_u{d,j}; % At = beta^{-1} I + C
%         model.Lat_u{d,j} = jitChol(model.At_u{d,j})';%lower bound
%         model.invLat_u{d,j} = model.Lat_u{d,j}\eye(size(model.Lat_u{d,j},1));  
%         model.invLatT_u{d,j} = model.invLat_u{d,j}';
%         model.logDetAt_u(d,j) = 2*(sum(log(diag(model.Lat_u{d,j})))); % log |At|
% 
%         model.P1_u{d,j} = model.invLat_u{d,j} * model.invLm_u; % M x M
%         P1TP1_u{d,j} = model.P1_u{d,j}'*model.P1_u{d,j};
%         % First multiply the two last factors; so, the large N is only involved
%         % once in the calculations (P1: MxM, Psi1':MxN, Y: NxD)
%         model.P_u{d,j} = model.P1_u{d,j} * (model.Psi1' * model.partm(:,d));%M*1  
% 
%         % Needed for both, the bound's and the derivs. calculations.
% %         model.TrPP_u = sum(sum(model.P_u .* model.P_u));
% %         model.TrPP_u = model.TrPP_u + model.W(index,j)^2*model.P_u{d,j}'*model.P_u{d,j};
%         model.B_u{d,j} = model.P1_u{d,j}' * model.P_u{d,j}; %M*1
%         
%         Tb_u = (1/model.beta) * model.W(index,j)^2 * (model.P1_u{d,j}' * model.P1_u{d,j})...
%         	 + model.W(index,j)^4 * model.B_u{d,j} * model.B_u{d,j}';
%          
% %          Tb_u = (1/model.beta) * model.W(index,j)^2 * (model.P1_u{d,j}' * model.P1_u{d,j});
% %         Tb_u = model.W(index,j)^4 * model.B_u{d,j} * model.B_u{d,j}';
% %         model.T1_u = model.T1_u - Tb_u;
% 
%         model.T1_u = model.T1_u + model.W(index,j)^2 * model.invK_uu - Tb_u; % sum w.r.t all D and all J
%         
%     end
% end

% Precompuations for the gradients
% gPsi0 = model.beta*(P1TP1_v*model.Psi0'*model.partm*[myOb(:,indexPresent); myMs]');


model.C_v = model.invLm_v * model.Psi4 * model.invLmT_v;
model.At_v = (1/model.beta) * eye(size(model.C_v,1)) + model.C_v; % At = beta^{-1} I + C
model.Lat_v = jitChol(model.At_v)';
model.invLat_v = model.Lat_v\eye(size(model.Lat_v,1));  
model.P1_v = model.invLat_v * model.invLm_v; % M x M
P1TP1_v = (model.P1_v' * model.P1_v);
model.P_v = model.P1_v * (model.Psi0' * model.m);
model.B_v = model.P1_v' * model.P_v;
Tb_v = (1/model.beta) * model.d * (model.P1_v' * model.P1_v);
     Tb_v = Tb_v + (model.B_v * model.B_v');
model.T1_v = model.d * model.invK_vv - Tb_v;

gPsi0_full = model.beta*(P1TP1_v*model.Psi0'*model.m*model.m');
gPsi0 = gPsi0_full(:,N+1:N+Nstar);
gPsi2 = -0.5 * model.beta * model.d;
gPsi4 = (model.beta/2) * model.T1_v;

samdPresent = intersect(samd,indexPresent);
lendPresent = length(samdPresent);
dpre = length(indexPresent);
% 


partYMs = yMs(:, samdPresent);
partYMs = partYMs - repmat(model.bias(:,samdPresent),size(partYMs,1),1);
partYMs = partYMs./repmat(model.scale(:,samdPresent),size(partYMs,1),1);

model.partm = [mOrig(:,samdPresent); myOb(:, samdPresent); partYMs];
% tempMatrix = calculateMatrixVwithPartM( model,samdPresent,samr );

gPsi1_full = zeros(N+Nstar,model.k);
gPsi3 = 0;
gPsi5 = zeros(model.k,model.k);

for d = 1:length(samdPresent)
    index = samdPresent(d);
    for j = 1:model.J
%         gPsi1_full = gPsi1_full + model.beta * model.W(index,j)^2 * (P1TP1_u{d,j}*model.Psi1'*model.partm(:,d)*model.partm(:,d)');%N*M
        gPsi1_full = gPsi1_full + model.beta*model.W(index,j)^2*model.partm(:,d)*model.partm(:,d)'*model.Psi1...
            *inv(model.W(index,j)^2*model.Psi5+model.K_uu/model.beta);
        
        gPsi3 = gPsi3 - 0.5 * model.beta * model.W(index,j)^2;
        
        gPsi5 = gPsi5 + 0.5*model.beta*model.W(index,j)^2*model.invK_uu...
            -0.5*model.beta*model.W(index,j)^4*inv(model.W(index,j)^2*model.Psi5+model.K_uu/model.beta)*...
            model.Psi1'*model.partm(:,d)*model.partm(:,d)'*model.Psi1*inv(model.W(index,j)^2*model.Psi5+model.K_uu/model.beta)...
            -0.5*model.W(index,j)^2*inv(model.W(index,j)^2*model.Psi5+model.K_uu/model.beta);
    end
end
gPsi1 = gPsi1_full(N+1:N+Nstar,:)';

gPsi1 = dpre/lendPresent * gPsi1;
gPsi3 = dpre/lendPresent * gPsi3;
gPsi5 = dpre/lendPresent * gPsi5;

[gKern0, gVarmeans0, gVarcovs0, gInd0] = kernVardistPsi1Gradient(model.kern_v, modelTest.vardist, model.X_v, gPsi0');
[gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern_u, modelTest.vardist, model.X_u, gPsi1');

[gKern2, gVarmeans2, gVarcovs2] = kernVardistPsi0Gradient(model.kern_v, modelTest.vardist, gPsi2);
[gKern3, gVarmeans3, gVarcovs3] = kernVardistPsi0Gradient(model.kern_u, modelTest.vardist, gPsi3);

[gKern4, gVarmeans4, gVarcovs4, gInd4] = kernVardistPsi2Gradient(model.kern_v, modelTest.vardist, model.X_v, gPsi4);
[gKern5, gVarmeans5, gVarcovs5, gInd5] = kernVardistPsi2Gradient(model.kern_u, modelTest.vardist, model.X_u, gPsi5);


gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2 + gVarmeans3 + gVarmeans4 + gVarmeans5;
gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2 + gVarcovs3 + gVarcovs4 + gVarcovs5; 
%

gPointDyn2 = zeros(Nstar, model.dynamics.q*2); 
% means
gPointDyn2(:,1:model.dynamics.q) = modelTest.dynamics.Kt'*(reshape(gVarmeansLik, Nstar, model.q) - modelTest.dynamics.vardist.means); 
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
gPointDyn2(:,(model.dynamics.q+1):(model.dynamics.q*2)) = gcovTmp.*modelTest.dynamics.vardist.covars; 

gPointDyn = gPointDyn1 + gPointDyn2;

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
