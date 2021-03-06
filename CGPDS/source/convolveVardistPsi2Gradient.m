function [gKern, gVarmeans, gVarcovars, gInd] = convolveVardistPsi2Gradient(convolveKern, vardist, Z, covGrad)

% copyright ZhaoJing
% last update 2013-09-17 
% variational means
N = size(vardist.means,1);
%  inducing variables 
[M Q] = size(Z); 
D = convolveKern.outputDimension;
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
% evaluate the kernel matrix 
[K, outKern, sumKern, Kgvar] = convolveVardistPsi2Compute(convolveKern, vardist, Z);

% inverse variances
% A = 1./convolveKern.Lambda_k;

% gradient wrt variance of the kernel 
gKernvar = 2*convolveKern.S.*reshape(sum(sum(Kgvar.*repmat(covGrad,[1,1,D]),1),2),[D,1]);  


% 1) line compute 0.5*(z_mq + z_m'q) for any q and store the result in a "M x Q x M" 
%  matrix where M is the number of inducing points and Q the latent dimension
% 2) line compute the z_mq - z_m'q, for any q
ZmZm  = zeros(M,Q,M);
ZmDZm = zeros(M,Q,M);
for q=1:size(Z,2)
  ZmZm(:,q,:) = 0.5*(repmat(Z(:,q),[1 1 M]) + repmat(reshape(Z(:,q),[1 1 M]),[M 1 1]));
  ZmDZm(:,q,:) = repmat(Z(:,q),[1 1 M]) - repmat(reshape(Z(:,q),[1 1 M]),[M 1 1]);
end
gInd = zeros(M,Q);
gLambda_k = zeros(1,Q);
gVarcovars = zeros(N,Q); 
gVarmeans = zeros(N,Q); 
partP2 = zeros(D,Q);

for d=1:D
    sigma_d=2*(convolveKern.P_d(d,:)+convolveKern.Lambda_k);
    covGrad_d = convolveKern.S(d)^2*(covGrad.*outKern(:,:,d));
    covGrad_d = reshape(covGrad_d,[M 1 M]);
    sumKern_d = reshape(sumKern(:,:,d),[M 1 M]);

    Amq = repmat(sigma_d.^-1,[M 1]);
    partInd1 = - 2*Amq.*sum(ZmDZm.*repmat(sumKern_d.*covGrad_d,[1 Q 1]),3);
    partInd2 = zeros(M,Q);
    ZmDZmA= bsxfun(@rdivide, ZmDZm.^2, sigma_d);% (z-z')^2 / 2(P+L);
    partP1(d,:) = sum(sum((ZmDZmA-1).*repmat(sumKern_d.*covGrad_d,[1 Q 1]),3),1).*sigma_d.^-1; 
    partL1 = partP1(d,:)+sum(sum(sumKern_d.*covGrad_d))./convolveKern.Lambda_k;
   gVarcovars_d = zeros(N,Q); 
   gVarmeans_d = zeros(N,Q);

    % Compute the gradient wrt lengthscales, variational means and variational variances  
    % For loop over training points  
    for n=find(vardist.map(:,d))'
        %
        %  
        mu_n = vardist.means(n,:); 
        s2_n = vardist.covars(n,:); 
        sA_n = 0.5*(convolveKern.P_d(d,:)+convolveKern.Lambda_k)+s2_n;% (p+L)/2+s_n
        MunZmZm = bsxfun(@minus,mu_n,ZmZm);% mu-(z+z')/2
        MunZmZmA =  bsxfun(@rdivide, MunZmZm, sA_n);% mu-(z+z')/2 / (p+L)/2+s_n
        k2Kern_n = sum(  bsxfun(@times, MunZmZmA.^2,sA_n),2);% (mu-(z+z')/2)^2 / (p+L)/2+s_n
        k2Kern_n = exp(-0.5*k2Kern_n)/prod(sqrt(sA_n));% element at index d,n with size of M*Q*M
        k2ncovG = repmat(k2Kern_n.*covGrad_d,[1 Q 1]); 
        tmp = MunZmZmA.*k2ncovG;
        tmp = sum(tmp,3);
        gVarmeans_d(n,:) = -(sum(tmp,1));
        partInd2 = partInd2 + tmp;
        MunZmZmA = MunZmZmA.*MunZmZm; 
        gVarcovars_d(n,:) = sum(sum( bsxfun(@times, (MunZmZmA - 1).*k2ncovG, 0.5./sA_n),1),3);

    end
    partP2(d,:)=0.5*sum(gVarcovars_d,1);%(p+L)/2+s_n
    gLambda_k=gLambda_k+partL1+partP2(d,:);
    gVarmeans=gVarmeans+gVarmeans_d;
    gVarcovars=gVarcovars+gVarcovars_d;
    gInd=gInd+partInd1 + partInd2; 
end
gP_d = partP1 + partP2; 

gKern = [gKernvar' gP_d(:)' gLambda_k];

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



    