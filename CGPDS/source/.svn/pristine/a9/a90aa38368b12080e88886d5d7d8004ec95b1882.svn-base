function [k0 k0part]=convolveVardistPsi0Compute(convolveKern, vardist)

% CONVOLVEPSI0COMPUTE description
% 
% Description


%   last update 2013-07-05 
%   last update 2013-10-29 
%   last update 2013-11-05
%   Copyright (c) ZhaoJing

%fprintf('#Calculating Psi0......\n');
%% 为了计算每一维的独立方差
if isfield(vardist,'chooseD')&& vardist.chooseD 
    convolveKern.outputDimension=1;
%     convolveKern.Lambda_k=convolveKern.Lambda_k(vardist.chooseD,:);
    convolveKern.P_d=convolveKern.P_d(vardist.chooseD,:);
    convolveKern.S=convolveKern.S(vardist.chooseD,:);
end
%%
N = size(vardist.means,1);
D = convolveKern.outputDimension;
sigma=repmat(convolveKern.Lambda_k,[D,1])+2*convolveKern.P_d;

%% 加入缺失数据情况
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
    Nvec=sum(vardist.map,1)';
    k0part=Nvec.*convolveKern.S.^2*sqrt(prod(convolveKern.Lambda_k))./sqrt(prod(sigma,2));% k0part is a vector composing of D elements
    k0=sum(k0part);   
% else
% %%
%     k0part=N*convolveKern.S.^2*sqrt(prod(convolveKern.Lambda_k))./sqrt(prod(sigma,2));% k0part is a vector composing of D elements
%     k0=sum(k0part);
% end

