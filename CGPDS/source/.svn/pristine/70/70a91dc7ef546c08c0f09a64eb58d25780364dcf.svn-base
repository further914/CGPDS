function [Psi2, Pnobias, Psi1, Psi2novar] = convolveconbiasVardistPsi2Compute(convolveKern, conbiasKern, vardist, Z)

% CONVOLVECONBIASVARDISTPSI2COMPUTE description.
%
%	Description:
%	[Psi2, Pnobias, Psi1] = convolveconbiasVardistPsi2Compute(rbfardKern, biasKern, vardist, Z)
%% 	convolveconbiasVardistPsi2Compute.m SVN version 583
% 	last update 2011-07-04T19:08:51.124563Z
%% 为了计算每一维的独立方差
if isfield(vardist,'chooseD')&& vardist.chooseD 
    conbiasKern.outputDimension=1;
end
%%
D = conbiasKern.outputDimension;
M = size(Z,1);
N  = size(vardist.means,1);
%% 加入缺失数据情况
if ~isfield(vardist,'map')|| isempty(vardist.map)
    vardist.map=ones(N,D);
end
[Psi1, Knovar]= convolveVardistPsi1Compute(convolveKern, vardist, Z); 
%% 为了后面计算S_d梯度时使用
Knovar=reshape(Knovar',M,N,D);
Knovar=permute(Knovar,[2,1,3]);
sumKnovar=sum(Knovar,1);
Psi2novar=repmat(sumKnovar,[M,1,1]);
Pnovarnobias=Psi2novar+permute(Psi2novar,[2,1,3]);
Psi2novar=conbiasKern.variance*Pnovarnobias; %3-d matrix with size of M*M*D

%%
sumPsi1 = sum(Psi1,1); 

Psi2 = ones(size(Z,1),1)*sumPsi1; 

Pnobias = Psi2 + Psi2'; 
Psi2 = conbiasKern.variance*Pnobias; 

% second naive way
%Psi22 = zeros(size(Z,1),size(Z,1));
%for j=1:size(Z,1)
%    for i=1:size(Z,1)
%        Psi22(j,i) = biasKern.variance*(sum(Psi1(:,j)) + sum(Psi1(:,i)));
%    end
%end
%sum(sum(abs(Psi2 - Psi22))) 

