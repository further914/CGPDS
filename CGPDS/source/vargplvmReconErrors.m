function [errStruct,Ypred,Ytrue] = vargplvmPredErrors(model, Y, Ytest, startInd, origBias, origScale,...
    missingInd, name, expNo,preday)

% VARGPLVMTAYLORANGLEERRORS Helper function for computing angle errors for CMU 35 data using
%
%	Description:
%	GPLVM with dynamics.
%	
%	Description:
%	[errStruct] = vargplvmTaylorAngleErrors(model, Y, Ytest, startInd, origBias, origScale,...
%	missingInd, name, expNo);
%	Based on fgplvmTaylorAngleErrors.m for FGPLVM.


%	Copyright (c) Andreas Damianou, Michalis Titsias 2010-2011 Neil Lawrence
% 	vargplvmTaylorAngleErrors.m SVN version 1767
% 	last update 2011-11-21T23:01:18.331219Z
%%%NEW
Ygplvm=Y;
YtestGplvm=Ytest;
YtestGplvmOrig = Ytest;

Y = Y - repmat(origBias, size(Y, 1), 1);
Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
Y = Y.*repmat(origScale, size(Y, 1), 1);
Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);
%%%


origYtest = Ytest;
YtrueTest = Ytest;

capName = name;
capName(1) = upper(capName(1));

temp = ones(1, size(Ytest, 2));
temp(missingInd) = 0;
presentInd = find(temp);

if ~isempty(presentInd)
    
    YtrainNn = Y(:, presentInd);
    YtestNn = origYtest(startInd:end, presentInd);



    disp('# NN prediction...');

    % First nearest neighbour
     dists = dist2(YtrainNn./(repmat(origScale(presentInd), size(YtrainNn, 1), 1)), ...
        YtestNn./(repmat(origScale(presentInd), size(YtestNn, 1), 1)));
    [void, bestIndNn] = min(dists);
    lenVal = size(Ytest, 1);
    err = (YtrueTest(startInd:end, missingInd) - Y(bestIndNn, missingInd))...
        ./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);
    err = err*180/pi;
    errStruct.angleErrorNnScaled = sqrt(mean(mean((err.*err))));
    load cmu35TaylorScaleBias
    err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
    errStruct.taylorErrorNnScaled = sum(sum(err.*err))/length(missingInd);

    % Second nearest neighbour
    dists = dist2(YtrainNn, ...
        YtestNn);
    [void, bestIndNn] = min(dists);
    lenVal = size(Ytest, 1);
    err = (YtrueTest(startInd:end, missingInd) - Y(bestIndNn, missingInd))...
        ./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);
    err = err*180/pi;
    errStruct.angleErrorNn = sqrt(mean(mean((err.*err))));
    load cmu35TaylorScaleBias
    err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
    errStruct.taylorErrorNn = sum(sum(err.*err))/length(missingInd);
    
    Ytest(startInd:end, missingInd) = NaN;


    %%%%%%%%%%%%%  VAR  GPLVM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    YtestGplvm(startInd:end, missingInd) = NaN;

    disp('# Initializing latent points...');
    indexPresent = setdiff(1:model.d, missingInd);
else
    YtrainNn = Y;
    YtestNn = origYtest(startInd:end, :);
    Ntrain=size(Y,1);
    index=[1:Ntrain];
    cycle=2;%一个周期的长度
    hcycle=cycle/2;
    number=Ntrain/cycle;%一共多少个周期
    mask=repmat([ones(1,hcycle),zeros(1,hcycle)],1,number);
    Ibefore=find(mask);
    Inext=setdiff(index,Ibefore);%前一天的下标
    Ybefore=YtrainNn(Ibefore,:);%前一天的数据
    Ynext=YtrainNn(Inext,:);
    Ybefore=reshape(Ybefore',[],number);%96*7;
    YtestNn=reshape(YtestNn',[],1);
    dists=dist2(Ybefore',YtestNn');
    [void, bestIndNn] = min(dists);
    mini=(bestIndNn-1)*cycle+1:(bestIndNn-1)*cycle+hcycle;%训练集中与测试集最近的点下标
    Ypred=Y(mini+hcycle,:);
    err = (YtrueTest(startInd:end, :) - Ypred)/origScale;
    err=err(1);                                                                                                                                                                                                 
    errStruct.ErrorNn = sqrt(mean(mean((err.*err))));
    
    
    Ytest = Ytest(1:startInd-1,:);
    %%%%%%%%%%%%%  VAR  GPLVM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    YtestGplvm = YtestGplvm(1:startInd-1,:);
    disp('# Initializing latent points...');
    indexPresent = [1:model.d];
    timeStampsTest=model.dynamics.t_star(end-startInd+1:end,:);
    model.dynamics.t_star=model.dynamics.t_star(1:end-startInd+1);
end



%iters=2000;% default: 1000
iters = model.reconstrIters;
fprintf(1,'# Optimising for %d iterations...\n',iters)

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
   if ~exist('mini')
       for i=1:size(YtestGplvm,1)
        % initialize the latent points using the nearest neighbour from the training data
        dst = dist2(YtestGplvm(i,indexPresent), Ygplvm(:,indexPresent));
        [mind, mini(i)] = min(dst);
       end
    end
    
    
    
    
    indexMissingData = startInd:size(YtestGplvm);
    
    vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = 0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
%     load('model.mat');
   
%     iters=1;
% preday=2;
modelType=model.type;
modelType(1) = upper(modelType(1));
if exist('preday') & preday
    fileName=['dem' modelType,name, num2str(expNo),'#',num2str(preday)];
else
    fileName=['dem' modelType,name, num2str(expNo)];
end
try 
    load(fileName);
catch  
        [x, varx] = vargplvmOptimiseSeqDyn(model, vardistx, YtestGplvm, 1, iters);
        % keep the optimized variational parameters
        barmu = x;
        lambda = varx;

        % Get the variational means and variacnes for the new test sequcen and
        % update the model to be prepared for prediction
        [x, varx, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(model, barmu, lambda, YtestGplvm);
        save(fileName,'modelUpdated','x','varx');
 
end
    % Latent variables corresponding to the data  with missing dimensions
    if ~isempty(indexMissingData)
        Testmeans = x(indexMissingData, :);
        Testcovars = varx(indexMissingData, :);
    else
        
        
        
        % model = vargplvmRestorePrunedModel(prunedModel, Ytr);
        % clear modelTr %%
        % model.y = Ytr;
        % model.m= gpComputeM(model); %%%
        % model.y=[];
        % load('model.mat');
%         Nstar = size(YtsOriginal,1);
        %%% NOTE: If you are reloading and the restoring the "modelUpdated" you will need to:
        %  modelUpdated.m = modelUpdated.mOrig;
        %  modelUpdated.P = modelUpdated.P1 * (modelUpdated.Psi1' * modelUpdated.m);
        %%
        
        % Prediction using the only information in the test time points
%         modelUpdated=model;
        [Testmeans Testcovars] = vargplvmPredictPoint(modelUpdated.dynamics, timeStampsTest);
        %% 下面是在隐空间上采用最近邻的方法
%         XtrainNn = modelUpdated.vardist.means(1:end-startInd+1, :);
%         XtestNn = x;        
%         Xbefore=YtrainNn(Ibefore,:);
%         Ynext=YtrainNn(Inext,:);
%         Xbefore=reshape(Xbefore',[],number);%96*7;
%         XtestNn=reshape(XtestNn',[],1);
%         dists=dist2(Xbefore',XtestNn');
%         [void, bestIndNn] = min(dists);
%         mini=(bestIndNn-1)*cycle+1:(bestIndNn-1)*cycle+hcycle;%训练集中与测试集最近的点坐标
        

%         Testmeans=[x;modelUpdated.vardist.means(mini+hcycle,:)];
%         Testcovars=[varx;modelUpdated.vardist.covars(mini+hcycle,:)];
        
       
%         Testmeans(97:192,:)=modelUpdated.vardist.means(mini+96,:);
%         Testcovars(97:192,:)=modelUpdated.vardist.covars(mini+96,:);
    end
    
    %  YtsOrMs = YtsOriginal(indexMissingData,:);
    disp('# Making the predictions for Y...')
    % Reconstruct all dimensions (missing and observed) for the test data
    mu = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
    Ypred = mu(startInd:end,:);
else
    %... STatic...
    % Static case doesnt capture data dependencies. Hence, we reconstruct
    % only the missing inputs:
    YtestGplvm = YtestGplvm(startInd:end, :);

    
   for i=1:size(YtestGplvm,1)
        % initialize the latent points using the nearest neighbour from the training data
        dst = dist2(YtestGplvm(i,indexPresent), Ygplvm(:,indexPresent));
        [mind, mini(i)] = min(dst);
    end

   
    vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = 0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
    % mu=[];
    % sigma=[];
    % for i=1:size(YtestGplvm,1)
    % optimize mean and vars of the latent point
    model.vardistx = vardistx;
    
    [x, varx] = vargplvmOptimisePoint(model, vardistx, YtestGplvm, 1, iters); %
    
    % reconstruct the missing outputs
    mu = vargplvmPosteriorMeanVar(model, x, varx);
    Ypred = mu;
    %    mu(i,:) = mu_i;
    %    sigma(i,:) = sigma_i;
    %end
    %save 'temp_temp.mat' %%%%
end

%{
[Xpred, varx] = vargplvmOptimiseSeqDyn(model, vardistx, YtestGplvm, 1, iters); % default: 1000
save('TEMPtaylorVargplvm2.mat','model','vardistx');
disp('# Getting the original variational parameters...');
barmu = vardistx.means;
lambda = vardistx.covars;
modelTest = vargplvmDynamicsUpdateModelTestVar(model, barmu, lambda);
Nstar = size(YtestGplvm,1);
overL=0;
Xpred = modelTest.vardist.means(end-(Nstar-1+overL):end,:);
varx= modelTest.vardist.covars(end-(Nstar-1+overL):end,:);
%}


%%%%%%%%%%%%%
% fName=['TEMPTaylorVargplvm' num2str(iters) 'it' num2str(expNo) '.mat'];
% fprintf(1,'# Saving the workspace in file %s...\n',fName);
% save(fName);

%%%%%%%%%


%%% Normalise predictions as in Fgplvm
% Ypred = Ypred - repmat(origBias, size(Ypred, 1), 1);
% Ypred = Ypred.*repmat(origScale, size(Ypred, 1), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




disp('# Final calculations...')

YtrueTest=YtestGplvmOrig;

% % De-normalise the predictions.
% err = (YtrueTest(startInd:end, missingInd) - Ypred(:, missingInd))./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);

err = YtrueTest(startInd:end, missingInd) - Ypred(:, missingInd);
err=err(1);
errStruct.ErrorGplvm = sqrt(mean(mean((err.*err))));
fprintf('errStruct.ErrorNn=%f\n',errStruct.ErrorNn);
fprintf('errStruct.ErrorGplvm=%f\n',errStruct.ErrorGplvm);
Ypred=Ypred(:, missingInd);
Ytrue= YtrueTest(startInd:end, missingInd);

%%%
