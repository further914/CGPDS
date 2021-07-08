function [gKern1, gKern2, gVarmeans, gVarcovars, gInd] = convolveconbiasVardistPsi2Gradient(convolveKern, conbiasKern, vardist, Z, covGrad)

% RBFARD2BIASVARDISTPSI2GRADIENT description.
%
%	Description:
%	[gKern1, gKern2, gVarmeans, gVarcovars, gInd] = convolveconbiasVardistPsi2Gradient(convolvedKern, conbiasKern, vardist, Z, covGrad)
%% 	rbfard2biasVardistPsi2Gradient.m SVN version 583
% 	last update 2011-07-04T19:08:51.164564Z
  
% variational means
N = size(vardist.means,1);
%  inducing variables 
[M Q] = size(Z); 
D=convolveKern.outputDimension;
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
[Psi2, Pnobias, Psi1, Psi2novar] = convolveconbiasVardistPsi2Compute(convolveKern, conbiasKern, vardist, Z);



% gradient wrt variance of the rbfard2 kernel 
% 这里要用三维矩阵转换
covGradtmp=repmat(covGrad,[1,1,D]);%先把系数扩展成3-d matrix
gS_d=reshape(sum(sum(Psi2novar.*covGradtmp,1),2),[1,D]);


%%
% gradient for the bias parameter  
gKern2 = sum(sum(Pnobias.*covGrad)); 
Bnm = conbiasKern.variance*ones(size(Psi1)); 
BPsi1Covg = Psi1.*(Bnm*covGrad); 

% compute the gradient wrt lengthscales, variational means and variational variances  
gVarmeans=zeros(N,Q);
gVarcovars=zeros(N,Q);
% gKernlengcs=zeros(1,Q);
% compute the gradient wrt lengthscales, variational means and variational variances  
for q=1:vardist.latentDimension
%
    S_q = vardist.covars(:,q);  
    Mu_q = vardist.means(:,q); 
    Z_q = Z(:,q)'; 
%     P_q=repmat(convolveKern.P_d(q),[N,1]);
%     L_q=repmat(convolveKern.Lambda_k(q),[N,1]);
    Aq=kron(convolveKern.P_d(:,q)+convolveKern.Lambda_k(q),ones(N,1))+repmat(S_q,[D,1]);% Aq=P+L+S
    
    B_q = repmat(repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1]),[D,1])./repmat(Aq, [1 M]);% B_q=(mu-z)/(P+L+S)
      
    % tmp is the nonGaussian term 这里的区别是我们是D维扩展矩阵
    tmp = B_q.*BPsi1Covg;
    % 乘以缺失点的映射关系
    tmp = tmp.*repmat(reshape(vardist.map,[],1),[1,M]);
    % variational means: you sum out the columns (see report)
    gVarmeanstmp(:,q) = -sum(tmp,2); 
    
    % 对扩展的D维进行分组求和
    gVarmeans(:,q)=sum(reshape(gVarmeanstmp(:,q),[N,D]),2);
    
    % inducing inputs: you sum out the rows 
    gInd(:,q) = sum(tmp,1)'; 
    
    % 
    %
    B_q = B_q.*repmat(repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1]),[D,1]);
   
    tmp=BPsi1Covg./repmat(Aq, [1 M]).*(B_q - 1);
    % 乘以缺失点的映射关系
    tmp = tmp.*repmat(reshape(vardist.map,[],1),[1,M]);
    % sum on M
    gVarcovarstmp(:,q) = sum(tmp,2);
    
    
    gVarcovars(:,q)=0.5*sum(reshape(gVarcovarstmp(:,q),[N,D]),2);% sum on D
    gP_d(:,q)=0.5*sum(reshape(gVarcovarstmp(:,q),[N,D]),1)';
    % gradient wrt L has a addtional term whch is L^-1 on nonGaussian term 
    tmp2=BPsi1Covg.*(repmat(Aq.^-1, [1 M]).*(B_q - 1)+convolveKern.Lambda_k(q)^-1);
    
    gLambda_k(:,q)=0.5*sum(sum(tmp2,2),1);
    %
end
%
gP_d=gP_d(:)';

gKern1 = [gS_d 2*gP_d 2*gLambda_k];

% gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarmeans = gVarmeans'; 
gVarmeans = 2*gVarmeans(:)'; 

% gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarcovars = gVarcovars'; 

gVarcovars = 2*gVarcovars(:)';

% gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
% this will unfold this matrix column-wise 
%gInd = gInd'; 
gInd = 2*gInd(:)'; 




