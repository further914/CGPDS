function [K, outKern, sumKern, Kgvar]=convolveVardistPsi2Compute(convolveKern, vardist, Z)

% CONVOLVEVARDISTPSI2COMPUTE description
%
% Description
% sumkern is the parts sensitive to n
% outkern is the parts with no relation to n
% K is a matrix

%%
%   last update 2013-07-03 
%   last update 2013-11-05
%   Copyright (c) ZhaoJing
% fprintf('#Calculating Psi2......\n');
%% 为了计算每一维的独立方差
if isfield(vardist,'chooseD')&& vardist.chooseD 
    convolveKern.outputDimension=1;
%     convolveKern.D_l=1;
%     convolveKern.Lambda_k=convolveKern.Lambda_k(vardist.chooseD,:);
    convolveKern.P_d=convolveKern.P_d(vardist.chooseD,:);
    convolveKern.S=convolveKern.S(vardist.chooseD,:);
end
%%

N= size(vardist.means,1);

M = size(Z,1); 

D = convolveKern.outputDimension;


if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
    

sumKern = zeros(M,M,D);
for d=1:D
    Ad = convolveKern.Lambda_k+convolveKern.P_d(d,:);
    for n=find(vardist.map(:,d))'
        %    
        AS_n = (0.5*Ad+vardist.covars(n,:));  %(L+P)/2+S
        normfactor =1/prod(sqrt(AS_n));%1/|(L+P)/2+S|^0.5

        Z_n = bsxfun(@minus, vardist.means(n,:), Z)*0.5;% (mu-z)/2

        Z_n = bsxfun(@times, Z_n, sqrt(AS_n.^-1)); 
        distZ = dist2(Z_n,-Z_n); % (mu-(z+z')/2)^2 / ((P+L)/2+S)
        % sumKern is a three-dimensional matrix 
        sumKern(:,:,d) = sumKern(:,:,d) + normfactor*exp(-0.5*distZ);  
        %
    end
    
ZZ =  bsxfun(@times, Z, sqrt(2*Ad).^-1);
distZZ = dist2(ZZ,ZZ);
normfactor =  prod(convolveKern.Lambda_k)/prod(sqrt(2*Ad));%|L|/|2(P+L)|^0.5
outKern(:,:,d) = normfactor*exp(-0.5*distZZ);
end
Kgvar = outKern.*sumKern; 
K = sum(repmat(reshape(convolveKern.S.^2,[1,1,D]),[M,M,1]).*Kgvar,3);




