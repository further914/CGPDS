function K = sdlfmvXsdlfmvKernComputeBlock(lfmKern1, lfmKern2, t1, t2, ...
    kyy, kyv, kvy, kvv, i, j, generalConst)

% SDLFMVXSDLFMVKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
%
%	Description:
%
%	K = SDLFMVXSDLFMVKERNCOMPUTEBLOCK(LFMKERN1, LFMKERN2, T1, T2, KYY,
%	KYV, KVY, KVV, I, J, GENERALCONSTANT) computes the kernel matrix for
%	the SDLFM kernel function in the block specified at indeces i,j. It
%	assumes the computation for functions that describe velocities
%	(velocity 1 and velocity 2).
%	 Returns:
%	  K - the kernel matrix portion of block i,j
%	 Arguments:
%	  LFMKERN1 - structure containing parameters for the system 1
%	  LFMKERN2 - structure containing parameters for the system 2
%	  T1 - times at which the system 1 is evaluated
%	  T2 - times at which the system 2 is evaluated
%	  KYY - covariance for the initial conditions between position 1 and
%	   position 2 at block i,j
%	  KYV - covariance for the initial conditions between position 1 and
%	   velocity 2 at block i,j
%	  KVY - covariance for the initial conditions between velocity 1 and
%	   position 2 at block i,j
%	  KVV - covariance for the initial conditions between velocity 1 and
%	   velocity 2 at block i,j
%	  I - interval to be evaluated for system 1
%	  J - interval to be evaluated for system 2
%	  GENERALCONSTANT - constants evaluated with
%	   sdlfmKernComputeConstant.m


%	Copyright (c) 2010. Mauricio A. Alvarez
% 	sdlfmvXsdlfmvKernComputeBlock.m SVN version 820
% 	last update 2011-06-16T07:23:44.000000Z

if nargin<11
    j = i;
    generalConst = [];
end

g1 = sdlfmvMeanCompute(lfmKern1(1), t1, 'Pos');
h1 = sdlfmvMeanCompute(lfmKern1(1), t1, 'Vel');
g2 = sdlfmvMeanCompute(lfmKern2(1), t2, 'Pos');
h2 = sdlfmvMeanCompute(lfmKern2(1), t2, 'Vel');

K = kyy*g1*g2.' + kyv*g1*h2.' + kvy*h1*g2.' + kvv*h1*h2.';

if i==j
    for k=1:length(lfmKern1)
        K  = K + lfmvXlfmvKernCompute(lfmKern1(k), lfmKern2(k), t1, t2);
    end
else    
    if i>j
        PosVel = zeros(1, length(t2));
        VelVel = zeros(1, length(t2));
        for k=1:length(lfmKern1)
            PosVel = PosVel + lfmvXlfmKernCompute(lfmKern2(k), lfmKern1(k), t2, lfmKern2(k).limit).'; 
            VelVel = VelVel + lfmvXlfmvKernCompute(lfmKern1(k), lfmKern2(k), lfmKern2(k).limit, t2);
        end
        if isempty(generalConst{i,j})
            K = K + g1*PosVel + h1*VelVel;        
        else
            K = K + (generalConst{i,j}(1,1)*g1 + generalConst{i,j}(2,1)*h1)*PosVel + ...
                (generalConst{i,j}(1,2)*g1 + generalConst{i,j}(2,2)*h1)*VelVel;           
        end 
    else
        PosVel = zeros(length(t1),1);
        VelVel = zeros(length(t1),1);
        for k =1:length(lfmKern1)
            PosVel = PosVel + lfmvXlfmKernCompute(lfmKern1(k), lfmKern2(k), t1, lfmKern1(k).limit);
            VelVel = VelVel + lfmvXlfmvKernCompute(lfmKern1(k), lfmKern2(k), t1, lfmKern1(k).limit);
        end
        if isempty(generalConst{i,j})
            K = K + PosVel*g2.' + VelVel*h2.';
        else
            K = K + PosVel*(generalConst{i,j}(1,1)*g2.' + generalConst{i,j}(2,1)*h2.') + ...
                VelVel*(generalConst{i,j}(1,2)*g2.' + generalConst{i,j}(2,2)*h2.');
        end
    end
end

