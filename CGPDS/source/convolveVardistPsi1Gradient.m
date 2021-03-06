function [gKern, gVarmeans, gVarcovars, gInd] = convolveVardistPsi1Gradient(convolveKern, vardist, Z, covGrad)

% all the matrix with the size of ND*M are as [Nd_1',Nd_2',...,Nd_D']'
% the different idea with rbfard2VardistPsi1Gradient is we must sum on
% D and some different details are the sigma are different and we have no
% kernlengcs while we have |sigma|^(-1/2) |Lambda|^(1/2)
%%
 
% variational means
[N,Q] = size(vardist.means);
%  inducing variables 
M = size(Z,1); 
D=convolveKern.outputDimension;
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
% evaluate the kernel matrix 
[K_fu Knovar] = convolveVardistPsi1Compute(convolveKern, vardist, Z);



% gradient wrt variance of the kernel 

gS_d = reshape(sum(sum(reshape((Knovar.*covGrad)',[M,N,D]),1),2),[1,D]); 
% multiply covGrad by psi1 
KfuCovGrad = K_fu.*covGrad;

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
    tmp = B_q.*KfuCovGrad;
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
   
    tmp=KfuCovGrad./repmat(Aq, [1 M]).*(B_q - 1);
    % 乘以缺失点的映射关系
    tmp = tmp.*repmat(reshape(vardist.map,[],1),[1,M]);
    % sum on M
    gVarcovarstmp(:,q) = sum(tmp,2);
    
    
    gVarcovars(:,q)=0.5*sum(reshape(gVarcovarstmp(:,q),[N,D]),2);% sum on D
    gP_d(:,q)=0.5*sum(reshape(gVarcovarstmp(:,q),[N,D]),1)';
    % gradient wrt L has a addtional term whch is L^-1 on nonGaussian term 
    tmp2=KfuCovGrad.*(repmat(Aq.^-1, [1 M]).*(B_q - 1)+convolveKern.Lambda_k(q)^-1);
    
    gLambda_k(:,q)=0.5*sum(sum(tmp2,2),1);
    %
end
%
gP_d=gP_d(:)';

gKern = [gS_d gP_d gLambda_k];

% gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarmeans = gVarmeans'; 
gVarmeans = gVarmeans(:)'; 

% gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarcovars = gVarcovars'; 

gVarcovars = gVarcovars(:)';

% gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
% this will unfold this matrix column-wise 
%gInd = gInd'; 
gInd = gInd(:)'; 

