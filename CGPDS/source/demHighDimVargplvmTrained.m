
% DEMHIGHDIMVARGPLVMTRAINED Perform tasks (predictions, sampling) for an already trained model on high dimensinoal video datasets.
%
%	Description:
%	DESC: Load an already trained model on a video dataset and do sampling or predictions with that.
%	!!! The parameters that concern the dataset and the experimentNo must be exactly the same as those for the trained model.
%	
%
%	See also
%	DEMHIGHDIMVARGPLVM3.M


%	Copyright (c) Michalis K. Titsias, 2010 - 2011 Andreas C. Damianou
% 	demHighDimVargplvmTrained.m SVN version 1570
% 	last update 2011-08-30T14:57:48.000000Z
%---- This file is in essence the same as demHighDimVargplvm3.m with the training phase omitted. ---
%-------------



% Fix seeds
clear timeStamps; % in case it's left from a previous experiment

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404;      end
if ~exist('itNo')         ,  itNo =2000;              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 49;          end     % Default: 49
if ~exist('latentDim')    ,  latentDim = 40;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 1;             end
% Set to 1 to keep only the dimensions modelling the head (missa dataset)
if ~exist('fixedBetaIters'), fixedBetaIters = 0;      end     % DEFAULT: 23
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('fixInd')        ,     fixInd = 0;                             end
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
% if ~exist('mappingKern')   ,     mappingKern = {'convolve', 'bias', 'white'}; end
if ~exist('mappingKern')   ,     mappingKern = {'rbfard2'}; end
if ~exist('reconstrIters') ,     reconstrIters = 1000;                   end
% 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
if ~exist('dataSetSplit')   ,    dataSetSplit = 'everyTwo';              end
if ~exist('blockSize')      ,    blockSize = 8;                          end
% 720x1280
if ~exist('dataSetName')    ,    dataSetName = 'missa';                  end
if ~exist('testReoptimise') ,    testReoptimise = 1;                     end
if ~exist('invWidthMultDyn'),    invWidthMultDyn = 100;                     end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('initX'),     initX ='ppca';   end
if ~exist('regularizeMeans'),  regularizeMeans = 0; end
if ~exist('initVardistIters'),  initVardistIters = 0; end

if ~exist('enableParallelism'), enableParallelism = 0 ; end
if ~exist('groupNumL'), groupNumL=3; end
% modified by ZhaoJing 
% if ~exist('enableParallelism'), enableParallelism = 1; end
% 2013-06-08  


% if strcmp(dataSetName, 'missa') & strcmp(dataSetSplit,'randomBlocks')
%     rand; % TMP (just to make the random seed more convenient! (this results in Ytr close to Yts).
% end
% if strcmp(dataSetName, 'ocean') & strcmp(dataSetSplit,'randomBlocks')
%     for i=1:5, rand; end % TMP (just to make the random seed more convenient! (this results in Ytr close to Yts).
% end

% dataSetName = 'susie';
% load susie_352_240

fprintf(1,'\n#----------------------------------------------------\n');
fprintf(1,'# Dataset: %s\n',dataSetName);
fprintf(1,'# ExperimentNo: %d\n', experimentNo);
fprintf(1,'# Latent dimensions: %d\n',latentDim);
fprintf(1,'# Iterations (with/without fixed Beta): %d / %s\n',fixedBetaIters, num2str(itNo));
fprintf(1,'# Reconstruction iterations: %d\n',reconstrIters);
fprintf(1,'# Tie Inducing points: %d\n',fixInd);
fprintf(1,'# InitX: %s\n',initX);
fprintf(1,'# Reoptimise inducing points (for reconstr): %d \n',testReoptimise);
fprintf(1,'# Dynamics used: %d\n', dynUsed);
if dynUsed
    fprintf(1,'# Dynamics kern: ');
    disp(dynamicKern);
end

fprintf(1,'# VardistCovarsMult: %d \n', vardistCovarsMult);
fprintf(1,'# InvWidthMultDyn: %d \n', invWidthMultDyn);
fprintf(1,'# InvWidthMult: %d \n', invWidthMult);
if exist('dataToKeep')  fprintf(1,'# DataToKeep: %d \n',dataToKeep); end
%fprintf(1,'#----------------------------------------------------\n');

switch dataSetName
    case 'missa' % There is a translation between frames 65...102
            Y = vargplvmLoadData(dataSetName);
            
            width=360;
            height=288;
    case 'missa36000' % There is a translation between frames 65...102
            Y = vargplvmLoadData(dataSetName);
            
            width=360;
            height=288;
    case 'missa-low' % There is a translation between frames 65...102
            data = vargplvmLoadData(dataSetName);
            Y=data.Y;
            width=data.W;
            height=data.H;
    case 'ocean'
        try
            load 'DATA_Ocean'
        catch
            Y=vargplvmLoadData('ocean');
        end
        width=1280; height=720;
    case 'horse'
        try
            load 'DATA_Horse'
        catch
            Y=vargplvmLoadData('horse');
        end
        Y = Y(90:end,:);
        width=249;
        height=187;
    case 'horse2'
        Y=vargplvmLoadData('horse2');
        width=320;
        height=240;
    otherwise
        % Y might have been loaded before this script is run. Otherwise
        % call lvmLoadData
        if ~(exist('Y') && exist('width') && exist('height'))
            try
                [Y, lbls] = vargplvmLoadData(dataSetName);
                height = lbls(1);
                width = lbls(2);
            catch
                load(dataSetName);
            end
        end
end

if exist('dataToKeep')
    Y = Y(1:dataToKeep,:);
end


% Crop video: keep only "window" of each frame
% vargplvmCropVideo



%%%%%%%TMP (for gradchek)
% itNo=2;
% indPoints=4;
% latentDim=3;
% dataSetSplit='halfAndHalf';
% reconstrIters = 2;
% Y = Y(1:10,:);
%%%%%%%%%



dims = size(Y,2);

trainModel = 1;% Set it to 1 to retrain the model. Set it to 0 to load an already trained one.

% Take a downsample version of the video for training and test on the
% remaining frames
N = size(Y,1);

fprintf(1,'# Preparing the dataset...\n');
switch dataSetSplit
    case 'everyTwo'
        % Training: 1,3,5... Test: 2,4,6...
        % assumes the number of data is even
        Nstar = N/2; Nts = N/2; Ntr = N/2;
        indTr = 1:2:N; indTs = 2:2:N;
    case 'halfAndHalf'
        % Training: 1,2,...N/2  Test: N/2+1, ...,N
        Nstar = round(N/2);
        indTr = 1:size(Y,1)-Nstar; indTs = size(Y,1)-Nstar+1:size(Y,1);
    case 'blocks'
        % Training: 1,2,..B, 2*B+1,.. 3*B,...
        % Test:B+1,...2*B,...   i.e. like "everyTwo" but with blocks
        lastBlockSize = mod(N,blockSize*2);
        mask = [ones(1,blockSize) zeros(1,blockSize)];
        mask = repmat(mask, 1, floor(floor(N/blockSize)/2));
        %md=mod(lastBlockSize,blockSize);
        %md2=lastBlockSize-md;
        %if mask(end)
        %    mask = [mask zeros(1,md2) ones(1,md)];
        %else
        %    mask = [mask ones(1, md2) zeros(1,md)];
        %end
        mask = [mask ones(1,lastBlockSize)]; % always end with the tr. set
        indTr = find(mask);
        indTs = find(~mask);
        if exist('msk')
            indTr = sort([indTr msk]);
            indTs=setdiff(1:N,indTr);
        end
        Nstar = size(indTs,2);
    case 'randomBlocks'
        mask = [];
        lastTrPts = 5;
        r=1; % start with tr. set
        while length(mask)<size(Y,1)-lastTrPts %The last lastTrPts will be from YTr necessarily
            blockSize = randperm(8);
            blockSize = blockSize(1);
            pts = min(blockSize, size(Y,1)-lastTrPts - length(mask));
            if r
                mask = [mask ones(1,pts)];
            else
                mask = [mask zeros(1,pts)];
            end
            r = ~r; % alternate between tr. and test set
        end
        mask = [mask ones(1,lastTrPts)];
        indTr = find(mask);
        indTs = find(~mask);
        if sum(sort([indTr indTs]) - (1:size(Y,1)))
            error('Something went wrong in the dataset splitting...');
        end
        % indPoints = length(indTr); %%%%% temp
        Nstar = length(indTs);
    case 'custom' %indTr and indTs must be provided
        Nstar = length(indTs);
end

if indPoints == -1
    indPoints = length(indTr);
end

Ytr = Y(indTr,:); Yts = Y(indTs,:);
t = linspace(0, 2*pi, size(Y, 1)+1)'; t = t(1:end-1, 1);
timeStampsTraining = t(indTr,1); timeStampsTest = t(indTs,1);


%-- For DEBUG
if exist('testOnTrainingTimes'), timeStampsTest = timeStampsTraining; end
if exist('testOnTrainingData'), Yts = Ytr(1:length(timeStampsTest),:); end
if exist('testOnReverseTrData'), Yts = Ytr(end:-1:1,:); timeStampsTest = timeStampsTest(1:size(Yts,1)); end
%--

YtsOriginal = Yts;

fprintf(1,'# Inducing points: %d\n',indPoints);
fprintf(1,'# Dataset size used (train/test) : %d / %d \n', size(Ytr,1), size(Yts,1));
fprintf(1,'# Dataset Split: %s ',dataSetSplit);
if strcmp(dataSetSplit,'blocks'),    fprintf(1,' (blockSize:%d)',blockSize); end
fprintf(1,'\n');
%fprintf(1,'# CropVideo / removeTranslation: %d / %d \n', cropVideo, removeTransl);
fprintf(1,'#----------------------------------------------------\n');

clear Y % Free up some memory

%{
% % Play movie
% for i=1:size(Ytr,1)
%     fr=reshape(Ytr(i,:),height,width);
%     imagesc(fr); colormap('gray');
%     pause(0.08);
% end
% %pause
% for i=1:size(Yts,1)
%     fr=reshape(Yts(i,:),height,width);
%     imagesc(fr); colormap('gray');
%     pause(0.08);
% end
%}

%%

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
options.numActive = indPoints;
if ~isempty(which('scg2'))
    options.optimiser = 'scg2';
else
    options.optimiser = 'scg';
end
% modified by ZhaoJing
DgtN=1; 
% 2013-06-08

if ~exist('DgtN') || ~DgtN
    options.enableDgtN = false;
end
d = size(Ytr, 2);



    % demo using the variational inference method for the gplvm model
    fprintf(1,'# Creating the model...\n');
    % scale = std(Ytr);
    % scale(find(scale==0)) = 1;
    %options.scaleVal = mean(std(Ytr));
    options.scaleVal = sqrt(var(Ytr(:)));
    if fixInd
        options.fixInducing=1;
        options.fixIndices=1:size(Ytr,1);
    end
    model = vargplvmCreate(latentDim, d, Ytr, options);
    
    
    
    % Temporary: in this demo there should always exist the mOrig field
    if ~isfield(model, 'mOrig')
        model.mOrig = model.m;
    end
    
    model = vargplvmParamInit(model, model.mOrig, model.X);
    
    %%% NEW
    model.kern.comp{1}.inputScales = invWidthMult./(((max(model.X)-min(model.X))).^2); % Default 5
    params = vargplvmExtractParam(model);
    model = vargplvmExpandParam(model, params);
    %%%
    
    
    model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
    %model.kern.comp{1}.variance = max(var(Y)); %%%
    
    
    %-------- Add dynamics to the model -----
    if dynUsed
        fprintf(1,'# Adding dynamics to the model...\n');
        optionsDyn.type = 'vargpTime';
        optionsDyn.t=timeStampsTraining;
        optionsDyn.inverseWidth=invWidthMultDyn; % Default: 100
        optionsDyn.testReoptimise = testReoptimise;
        optionsDyn.initX = initX;
        optionsDyn.regularizeMeans = regularizeMeans;
        
        kern = kernCreate(t, dynamicKern); % Default: {'rbf','white','bias'}
           
                
        %-- Initialize each element of the compound kernel (optional but
        % recommended) 
        % ATTENTION: For the gradients we assume that the base kernel (rbf,
        % matern etc) must be the FIRST one and if a second base kernel
        % (not white or bias) exist must be the LAST one!!!!!!!!!!!!!!
        vargplvmInitDynKernel;
        
        
        optionsDyn.kern = kern;
        optionsDyn.vardistCovars = vardistCovarsMult; % 0.23 gives true vardist.covars around 0.5 (DEFAULT: 0.23) for the ocean dataset
        
        % Fill in with default values whatever is not already set
        optionsDyn = vargplvmOptionsDyn(optionsDyn);
        model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
        
        fprintf(1,'# Further calibration of the initial values...\n');
        model = vargplvmInitDynamics(model,optionsDyn);
        
        %___NEW TEMP: to also not learn the last kernel's variance
        if numel(kern.comp) > 1 && exist('learnSecondVariance') && ~learnSecondVariance
            fprintf(1,'# The variance for %s in the dynamics is not learned!\n',kern.comp{end}.type)
            model.dynamics.learnSecondVariance = 0;
            model.dynamics.kern.comp{end}.inverseWidth = model.dynamics.kern.comp{1}.inverseWidth/10; %%% TEMP
        end
        %____
    end
    
    if model.N > 50 && enableParallelism
        fprintf('# Parallel computations w.r.t the datapoints!\n');
         model.vardist.parallel = 1;
    end
    
    %%%
    model.dataSetInfo.dataSetName = dataSetName;
    model.dataSetInfo.dataSetSplit = dataSetSplit;
    if strcmp(dataSetSplit, 'custom')
        model.dataSetInfo.indTr = indTr;
        model.dataSetInfo.indTs = indTs;
    end
    %%%
    
    model.beta=1/(0.01*var(model.mOrig(:)));
    modelInit = model;
    
   % disp(model.vardist.covars)
   % disp(model.kern.comp{1}.inputScales)
    
    
    %%---
    capName = dataSetName;
    capName(1) = upper(capName(1));
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
    %%---
    
    if exist('experimentDesc'), model.expDesc = experimentDesc; end
    
    display = 1;
%%
load(fileToSave);

model = vargplvmRestorePrunedModel(prunedModel, Ytr);
clear modelTr %%
model.y = Ytr;
model.m= gpComputeM(model); %%%
model.y=[];
Nstar = size(YtsOriginal,1);
%%% NOTE: If you are reloading and the restoring the "modelUpdated" you will need to:
%  modelUpdated.m = modelUpdated.mOrig;
%  modelUpdated.P = modelUpdated.P1 * (modelUpdated.Psi1' * modelUpdated.m);
%%

% Prediction using the only information in the test time points
[Testmeans2 Testcovars2] = vargplvmPredictPoint(model.dynamics, timeStampsTest);
Varmu2 = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2);
% Mean absolute error per pixel
errorOnlyTimes = mean(abs(Varmu2(:) - YtsOriginal(:)));

%%%%% CHANGE THE NN PREDICTION, it's not correct like this when predicting
%%%%% ahead (i.e. when the dataSet split is not 'everyTwo'
% 2-nearest (in time) prediction
% NNmu = zeros(Nstar, model.d);
% for i=1:Nstar
%     if i < Nstar
%         NNmu(i, :) = 0.5*(Ytr(i,:) + Ytr(i+1,:));
%     else
%         NNmu(i, :) = Ytr(end,:);
%     end
% end
% Mean absolute error per pixel
%errorNN = mean(abs(NNmu(:) - YtsOriginal(:)));
%errorNN = sum(sum( abs(NNmu - YtsOriginal)) )/prod(size(YtsOriginal)); %equivalent

fprintf(1,'# Error GPLVM: %d\n', errorOnlyTimes);
