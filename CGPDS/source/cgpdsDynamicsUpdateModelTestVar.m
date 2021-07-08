function [muX, varX, model]  = cgpdsDynamicsUpdateModelTestVar(model, x, varx, y)        

% cgpdsDynamicsUpdateModelTestVar return the original variational means and
%
%	Description:
%	variances for the test data and update the model to be prepared for
%	prediction
%	DESC in the dynamical vargplvm the original variational parameters are
%	reparametrized so that the optimiser optimises the new free ones and then
%	these are mapped back to the original ones. This function takes the new
%	free parameters which are optimised and returns the original ones. Also,
%	because there is coupling in the variational parameters (because the
%	original ones are obtained via a nonlinear operation on the free ones
%	also involving Kt) after introducing test points (which have their own
%	latent positions) the whole variational distribution is changed.
%	Therefore, the updated model is also returned.
%	ARG model : The cgpds model which contains dynamics.
%	ARG x : the optimised variational means (the free parameters, not the ones) of the test points
%	ARG varx : the optimised variational variances (the free parameters, not the
%	ones) of the test points
%	ARG y: the test points
%	RETURN muX : the original,true variational means corresponding to x
%	RETURN varX : the original,true variational variances corresponding to varx
%	RETURN model : the updated model, prepared for prediction
% 	vargplvmDynamicsUpdateModelTestVar.m SVN version 1413
% 	last update 2011-07-04T19:08:52.084563Z
%
% COPYRIGHT : Michalis K. Titsias, 2009-2011
% COPYRIGHT : Neil D. Lawrence, 2009-2011
%
  
Nstar = size(x,1);    
N = model.N;
mask = sum(isnan(y),2); 
indexObservedData = find(mask==0)'; 
indexMissingData = setdiff(1:Nstar, indexObservedData);


% re-order test data so that observed are first and then are the missing 
Order = [indexObservedData, indexMissingData];
model.dynamics.t_star = model.dynamics.t_star(Order);
if isfield(model.dynamics, 'seq') && ~isempty(model.dynamics.seq) 
   mask = sum(isnan(y),2); 
   indexObservedData = find(mask==0)'; 
   yOb = y(indexObservedData, :);
   
   modelTest = model;
   modelTest.dynamics.t = model.dynamics.t_star;
   modelTest.dynamics.vardist.means = x;
   modelTest.dynamics.vardist.covars = varx; 
   modelTest.dynamics.N = size(modelTest.dynamics.vardist.means, 1);
   modelTest.vardist.numData = modelTest.dynamics.N;
   modelTest.vardist.nParams = 2*prod(size(modelTest.dynamics.vardist.means));
   Kt = kernCompute(model.dynamics.kern, model.dynamics.t_star);
   modelTest.dynamics.Kt = Kt;
   modelTest.dynamics.fixedKt = 1;
   % This enables the computation of the variational means and covars
   modelTest.dynamics.seq = []; 
   modelTest = vargpTimeDynamicsUpdateStats(modelTest);
   
   % Variational means and covars for the test block data
   % latent variables
   muX = modelTest.vardist.means;
   varX = modelTest.vardist.covars;
   
   if ~isempty(yOb)     
      % Since the test block also might have a fully observed subblock, you
      % need to update the model structure (this for the subseqeucent prediction
      % using the PostMeanVar function)
      vardistOb = modelTest.vardist; 
      vardistOb.means = modelTest.vardist.means(indexObservedData, :);
      vardistOb.covars = modelTest.vardist.covars(indexObservedData, :);
      vardistOb.nParams = 2*prod(size(vardistOb.means));
      vardistOb.numData = size(vardistOb.means,1);
      % Psi statistics for the data of the new block/sequence which have a fully
      % observed features/dimensions
      obsPsi0 = kernVardistPsi1Compute(model.kern_v, vardistOb, model.X_v);
      obsPsi1 = kernVardistPsi1Compute(model.kern_u, vardistOb, model.X_u);
      obsPsi2 = kernVardistPsi0Compute(model.kern_v, vardistOb);
      obsPsi3 = kernVardistPsi0Compute(model.kern_u, vardistOb);      
      obsPsi4 = kernVardistPsi2Compute(model.kern_v, vardistOb, model.X_v);
      obsPsi5 = kernVardistPsi2Compute(model.kern_u, vardistOb, model.X_u);
      
      model.N = model.N + size(yOb,1); 
      myOb = yOb - repmat(model.bias,size(yOb,1),1); 
      myOb = myOb./repmat(model.scale,size(yOb,1),1); 
      model.m = [model.m; myOb];
          
      model.dynamics.t = [model.dynamics.t; model.dynamics.t_star(indexObservedData)];
      Kt = zeros(N+size(indexObservedData,2),N+size(indexObservedData,2));
      Kt(1:N,1:N) = model.dynamics.Kt; 
      tmpKt = modelTest.dynamics.Kt(indexObservedData,:); 
      tmpKt = tmpKt(:,indexObservedData);
      Kt(N+1:end, N+1:end) = tmpKt; 
      model.dynamics.Kt = Kt;
      model.dynamics.seq = [model.dynamics.seq, (N+size(indexObservedData,2))];
      model.Psi0 = [model.Psi0; obsPsi0]; 
      model.Psi1 = [model.Psi1; obsPsi1]; 
      model.Psi2 = model.Psi2 + obsPsi2;
      model.Psi3 = model.Psi3 + obsPsi3;
      model.Psi4 = model.Psi4 + obsPsi4;
      model.Psi5 = model.Psi5 + obsPsi5;
      
      % dynamics (reparametrized) variational parameters (barmu and lambdas)
      model.dynamics.vardist.means = [model.dynamics.vardist.means; x(indexObservedData,:)];
      model.dynamics.vardist.covars = [model.dynamics.vardist.covars; varx(indexObservedData,:)]; 
      model.dynamics.N = size(model.dynamics.vardist.means, 1); 
      model.dynamics.vardist.transforms.index = (size(model.dynamics.vardist.means,2)*model.N+1):size(model.dynamics.vardist.means,2)*model.N*2;
      model.dynamics.vardist.numData =  size(model.dynamics.vardist.means, 1); 
      model.dynamics.vardist.nParams = 2*prod(size(model.dynamics.vardist.means));
      model.dynamics.nParams=model.dynamics.vardist.nParams +model.dynamics.kern.nParams;
      
      % the actual variational distrbution 
      model.vardist.means = [model.vardist.means; vardistOb.means];
      model.vardist.covars = [model.vardist.covars; vardistOb.covars]; 
      model.vardist.transforms.type = model.vardist.transforms.type;
      model.vardist.transforms.index = (size(model.vardist.means,2)*model.N+1):size(model.vardist.means,2)*model.N*2;
      model.vardist.numData =  size(model.vardist.means, 1); 
      model.vardist.nParams = 2*prod(size(model.vardist.means));
      
      % update in the model structure useful for prediction
      model.C_v = model.invLm_v * model.Psi2 * model.invLmT_v;
      model.TrC_v = sum(diag(model.C_v)); % Tr(C)
      model.At_v = (1/model.beta) * eye(size(model.C_v,1)) + model.C_v; % At = beta^{-1} I + C 
      model.Lat_v = jitChol(model.At_v)';
      model.invLat_v = model.Lat_v\eye(size(model.Lat_v,1));  
      model.invLat_vT = model.invLat_v';
      model.logDetAt_v = 2*(sum(log(diag(model.Lat_v)))); % log |At|
      model.P1_v = model.invLat_v * model.invLm_v; % M x M
      model.P_v = model.P1_v * (model.Psi1' * model.m);
      model.TrPP_v = sum(sum(model.P_v .* model.P_v));
   end
else
    
   %for all data(including training data and testing data) 
   % Augment the time vector to include the timestamp of the new point
   model.dynamics.t = [model.dynamics.t; model.dynamics.t_star];
   % Augment the reparametrized variational parameters mubar and lambda
   model.dynamics.vardist.means = [model.dynamics.vardist.means; x];
   model.dynamics.vardist.covars = [model.dynamics.vardist.covars; varx]; 
   model.dynamics.N =  model.dynamics.N + Nstar;
   model.vardist.numData = model.dynamics.N;
   model.dynamics.vardist.numData = model.dynamics.N;
   model.vardist.nParams = 2*prod(size(model.dynamics.vardist.means));
   model.dynamics.vardist.nParams = 2*prod(size(model.dynamics.vardist.means));
   model.dynamics.seq = model.dynamics.N; %%%%%% Chec
   model = vargpTimeDynamicsUpdateStats(model);    
   
   %for training data
   vardist2 = model.vardist;
   vardist2.means = model.vardist.means(1:end-Nstar,:);%for training data
   vardist2.covars = model.vardist.covars(1:end-Nstar,:);
   vardist2.nParams = 2*prod(size(vardist2.means));
   vardist2.numData = size(vardist2.means,1);
   
   %for training data
   model.Psi0 = kernVardistPsi1Compute(model.kern_v, vardist2, model.X_v);
   model.Psi1 = kernVardistPsi1Compute(model.kern_u, vardist2, model.X_u);
   model.Psi2 = kernVardistPsi0Compute(model.kern_v, vardist2);
   model.Psi3 = kernVardistPsi0Compute(model.kern_u, vardist2);
   model.Psi4 = kernVardistPsi2Compute(model.kern_v, vardist2, model.X_v);
   model.Psi5 = kernVardistPsi2Compute(model.kern_u, vardist2, model.X_u);
   
   %for training data
   model.C_v = model.invLm_v * model.Psi2 * model.invLmT_v;
   model.TrC_v = sum(diag(model.C_v)); % Tr(C)
   model.At_v = (1/model.beta) * eye(size(model.C_v,1)) + model.C_v; % At = beta^{-1} I + C 
   model.Lat_v = jitChol(model.At_v)';
   model.invLat_v = model.Lat_v\eye(size(model.Lat_v,1));  
   model.invLat_vT = model.invLat_v';
   model.logDetAt_v = 2*(sum(log(diag(model.Lat_v)))); % log |At|
   model.P1_v = model.invLat_v * model.invLm_v; % M x M
   model.P_v = model.P1_v * (model.Psi1' * model.m);
   model.TrPP_v = sum(sum(model.P_v .* model.P_v));
   
   % Variational means and covars for the test block data
   % latent variables
   %for testing data
   muX = model.vardist.means(N+1:end,:);
   varX = model.vardist.covars(N+1:end,:);   
   
end
   
model.X = model.vardist.means;

% function result=splice(A,N,M,D,B,Ns)
%     tt=reshape(A',[M,N,D]);
%     tt=permute(tt,[2,1,3]);
%     tt=reshape(tt,N,[]);
%     nn=reshape(B',[M,Ns,D]);
%     nn=permute(nn,[2,1,3]);
%     nn=reshape(nn,Ns,[]);
%     tt=[tt;nn];
%     tt=reshape(tt,[N+Ns,M,D]);
%     tt=permute(tt,[2,1,3]);
%     result=reshape(tt,M,[])';
