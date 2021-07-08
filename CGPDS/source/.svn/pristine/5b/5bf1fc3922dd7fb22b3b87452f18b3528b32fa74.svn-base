
% DEMBRENDANVARGPLVM3 Run variational GPLVM on Brendan face data.
%


%	Description:
%   Copyright (c) 2013 ZhaoJing
%	% 	demBrendanVargplvm3.m SVN version 1509


% 	last update 2011-08-03T14:12:49.214877Z
%   last update 2013-5-28 

% Fix seeds
fprintf('==============================================================\n');
fprintf('GP-LVM\n');
fprintf('--------------------------------------------------------------\n');
disp('Run variational GPLVM on Brendan face data!');
fprintf('AISTATS\n');
fprintf('Zhao Jing\n');
fprintf('==============================================================\n');
fprintf('\n');



randn('seed', 1e6);
rand('seed', 1e6);

dataSetName = 'brendan';
experimentNo = 3;
printDiagram = 1;



% load data
disp('Loading dataset......');
[Y, lbls] = lvmLoadData(dataSetName);
disp('Loading finished!')
% training and test sets
Ntr = 1000; 
perm = randperm(size(Y,1)); 
Ytr = Y(perm(1:Ntr),:);      %lblsTr = lbls(perm(1:Ntr),:);
Yts = Y(perm(Ntr+1:end),:);  %lblsTs = lbls(perm(Ntr+1:end),:);
reTrainModels=1;
rePredict=0;
reDraw=0;
if reTrainModels
    % Set up model
    options = vargplvmOptions('dtcvar');
    options.kern = {'rbfard2', 'white'};
    options.numActive = 50; 
    options.scale2var1 = 1; % scale data to have variance 1
    %options.tieParam = 'tied';  

    options.optimiser = 'scg';
    latentDim = 30;
    d = size(Y, 2);

    disp('Create Model!');
    model = vargplvmCreate(latentDim, d, Ytr, options);
    %
    disp('Init Model!');
    model = vargplvmParamInit(model, model.m, model.X); 

    iters = 100;
    display = 1;
    disp('Start Optimising!');
    fprintf('Total cycles are %d !\n ',iters);
    model = vargplvmOptimise(model, display, iters);
    capName = dataSetName;
    capName(1) = upper(capName(1));
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    save(['mat_files\dem' capName modelType num2str(experimentNo) '.mat'], 'model');
    disp('Training finished!');
end 
if rePredict
    
    iters = 1;
    display = 1;

   load('demBrendanVargplvm3.mat');
    %%

    % 50% missing outputs from the each test point
    numIndPresent = round(0.5*model.d);
    indexP = [];
    Init = [];
    Testmeans = []; 
    Testcovars = [];
    Varmu = [];
    Varsigma = [];
    YtsOriginal = Yts; %%%NEW
    % patrial reconstruction of test points
    for i=1:size(Yts,1)
        %
        % randomly choose which outputs are present
        permi = randperm(model.d);
        indexPresent =  permi(1:numIndPresent);
        indexP(i,:) = indexPresent;
        indexMissing = setdiff(1:model.d, indexPresent); %
        Yts(i,indexMissing) = NaN; %
        % initialize the latent point using the nearest neighbour 
        % from he training data
        dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
        [mind, mini] = min(dst);

        Init(i,:) = model.vardist.means(mini,:);
        % create the variational distribtion for the test latent point
        vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
        vardistx.covars = 0.2*ones(size(vardistx.covars));

        % optimize mean and vars of the latent point 
        model.vardistx = vardistx;
    %    [x, varx] = vargplvmOptimisePoint(model, vardistx, Yts(i, indexPresent), indexPresent, display, iters); %old
         [x, varx] = vargplvmOptimisePoint(model, vardistx, Yts(i, :), display, iters); %
        Testmeans(i,:) = x;
        Testcovars(i,:) = varx;

        % reconstruct the missing outputs  
        [mu, sigma] = vargplvmPosteriorMeanVar(model, x, varx);
        Varmu(i,:) = mu; 
        Varsigma(i,:) = sigma; 
        %
    end
    capName = dataSetName;
    capName(1) = upper(capName(1));
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    save(['mat_files\dem' capName modelType num2str(experimentNo) 'Pred.mat'], 'model', 'perm', 'indexP', 'Varmu', 'Varsigma');
end
if reDraw
    load('demBrendanVargplvm3Pred.mat');
    if exist('printDiagram') & printDiagram
      lvmPrintPlot(model, lbls, dataSetName, experimentNo);
    end

    % Load the results and display dynamically.
    lvmResultsDynamic(model.type, dataSetName, experimentNo, 'image', [20 28], 1, 0, 1)
end 

    

