function K = sdlfmXsdlfmKernComputeBlock(lfmKern1, lfmKern2, t1, t2, ...
    kyy, kyv, kvy, kvv, i, j, generalConst)

% SDLFMXSDLFMKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
%
%	Description:
%
%	K = SDLFMXSDLFMKERNCOMPUTEBLOCK(LFMKERN1, LFMKERN2, T1, T2, KYY,
%	KYV, KVY, KVV, I, J, GENERALCONSTANT) computes the kernel matrix for
%	the SDLFM kernel function in the block specified at indeces i,j. It
%	assumes the computation for functions that describe positions
%	(position 1 and position 2).
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
% 	sdlfmXsdlfmKernComputeBlock.m SVN version 820
% 	last update 2010-05-28T06:01:33.000000Z

if nargin<11
    j = i;
    generalConst = [];
end

c1 = sdlfmMeanCompute(lfmKern1(1), t1, 'Pos');
e1 = sdlfmMeanCompute(lfmKern1(1), t1, 'Vel');
c2 = sdlfmMeanCompute(lfmKern2(1), t2, 'Pos');
e2 = sdlfmMeanCompute(lfmKern2(1), t2, 'Vel');

K = kyy*c1*c2.' + kyv*c1*e2.' + kvy*e1*c2.' + kvv*e1*e2.';

if i==j
    for k=1:length(lfmKern1)
        K  = K + lfmXlfmKernCompute(lfmKern1(k), lfmKern2(k), t1, t2);
    end
else
    if i>j
        PosPos = zeros(1,length(t2));
        VelPos = zeros(1,length(t2));
        for k=1:length(lfmKern1)
            PosPos = PosPos + lfmXlfmKernCompute(lfmKern1(k), lfmKern2(k), ...
                lfmKern2(k).limit, t2);
            VelPos = VelPos + lfmvXlfmKernCompute(lfmKern1(k), lfmKern2(k), ...
                lfmKern2(k).limit, t2);
        end
        if isempty(generalConst{i,j})
            K = K + c1*PosPos + e1*VelPos;
        else
            K = K + (generalConst{i,j}(1,1)*c1 + generalConst{i,j}(2,1)*e1)*PosPos + ...
                (generalConst{i,j}(1,2)*c1 + generalConst{i,j}(2,2)*e1)*VelPos;
        end 
    else
        PosPos = zeros(length(t1),1);
        PosVel = zeros(length(t1),1);
        for k=1:length(lfmKern1)
            PosPos = PosPos + lfmXlfmKernCompute(lfmKern1(k), lfmKern2(k), t1, lfmKern1(k).limit);
            PosVel = PosVel + lfmvXlfmKernCompute(lfmKern2(k), lfmKern1(k), lfmKern1(k).limit, t1).';
        end        
         if isempty(generalConst{i,j})
             K = K + PosPos*c2.' + PosVel*e2.';               
         else
             K = K + PosPos*(generalConst{i,j}(1,1)*c2.' + generalConst{i,j}(2,1)*e2.') + ...
                 PosVel*(generalConst{i,j}(1,2)*c2.' + generalConst{i,j}(2,2)*e2.');
         end
    end
end

