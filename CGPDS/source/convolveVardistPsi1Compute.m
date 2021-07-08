function [K,Knovar]=convolveVardistPsi1Compute(convolveKern, vardist, Z)

%   last update 2013-07-03 
%   last update 2013-11-05
%   Copyright (c) ZhaoJing

%% 为了计算每一维的独立方差
if isfield(vardist,'chooseD')&& vardist.chooseD 
    convolveKern.outputDimension=1;
%     convolveKern.Lambda_k=convolveKern.Lambda_k(vardist.chooseD,:);
    convolveKern.P_d=convolveKern.P_d(vardist.chooseD,:);
    convolveKern.S=convolveKern.S(vardist.chooseD,:);
end
%%

M = size(Z,1); 
D = convolveKern.outputDimension;
% for d=1:D
    N  = size(vardist.means,1);
    if ~isfield(vardist,'map')|| isempty(vardist.map)
        vardist.map=ones(N,D);
    end
    argExp = zeros(N*D,M); 
    normfactor = ones(N*D,1);
    for q=1:vardist.latentDimension
    %
        S_q = repmat(vardist.covars(:,q),[D,1]); 
        sigma_q=kron(convolveKern.Lambda_k(q)+convolveKern.P_d(:,q),ones(N,1))+S_q;
        normfactor = normfactor.*sigma_q;
        Mu_q = vardist.means(:,q); 
        Z_q = Z(:,q)';
        distan = (repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1])).^2;
        argExp = argExp + repmat(sigma_q.^-1, [1 M]).*repmat(distan,[D,1]);
    %
    end
    normfactor = sqrt(prod(convolveKern.Lambda_k))./normfactor.^0.5;%|L|^0.5/|P+L+S|^0.5

    Knovar = repmat(normfactor,[1 M]).*exp(-0.5*argExp); 
    K = kron(convolveKern.S,ones(N,M)).*Knovar;
%     if isfield(vardist,'map')&& ~isempty(vardist.map)
        Knovar = repmat(reshape(vardist.map,[],1),1,M).*Knovar;
        K = repmat(reshape(vardist.map,[],1),1,M).*K;
%     end
% end
