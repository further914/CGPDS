%	Copyright (c) Michalis K. Titsias, 2010 - 2011 Andreas C. Damianou
%	Copyright (c) 2011 Andreas C. Damianou 
%	Copyright (c) 2015 Jing Zhao
clear timeStamps; % in case it's left from a previous experiment

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);
%% 修改临时变量区域
% experimentNo = 60;
cut='fix';
indexMissing=[4];
% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 0;      end
if ~exist('constructNo') ,  constructNo = 1;      end
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
if ~exist('mappingKern')   ,     mappingKern = {'convolve'}; end
if ~exist('modelType'),          modelType='varmgplvm'; end
if ~exist('reconstrIters') ,     reconstrIters = 1000;                   end
% 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
if ~exist('dataSetSplit')   ,    dataSetSplit = 'everyTwo';              end
if ~exist('blockSize')      ,    blockSize = 8;                          end
% 720x1280
if ~exist('dataSetName')    ,    dataSetName = 'missa';                  end
if ~exist('testReoptimise') ,    testReoptimise = 0;                     end
if ~exist('invWidthMultDyn'),    invWidthMultDyn = 100;                     end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('initX'),     initX ='ppca';   end
if ~exist('regularizeMeans'),  regularizeMeans = 0; end
if ~exist('initVardistIters'),  initVardistIters = 0; end

if ~exist('enableParallelism'), enableParallelism = 0 ; end
%% 关于训练，重构优化，重构，画图，续优化的相关参数的设置



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
fprintf(1,'# ModelType: %s\n',modelType);
fprintf(1,'# ExperimentNo: %d\n', experimentNo);
fprintf(1,'# ConstructNo: %d\n', constructNo);
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
    case 'missa-low' % There is a translation between frames 65...102
        data = vargplvmLoadData(dataSetName);
        Y=data.Y;
        width=data.W;
        height=data.H;
    case 'sdata' % There is a translation between frames 65...102
        data = vargplvmLoadData(dataSetName);
        Y=data.Y;
        width=data.W;
        height=data.H;
    case 'frey_face' % There is a translation between frames 65...102
        data = vargplvmLoadData(dataSetName);
        Y=data.Y;
        width=data.W;
        height=data.H;
    case 'cmu35WalkJog'
        [Y, lbls] = lvmLoadData(dataSetName);
        seq = cumsum(sum(lbls)) - [1:31];
        dataSetName = 'cmu35gplvm';
        % load data
        [Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
        % Fix times:
        prevSeq = 1;
        timeStampsTraining = [];
        dt=0.05;
        for i=1:length(seq)
            t = ([0:(seq(i)-prevSeq)].*dt)';
            prevSeq = seq(i)+1;
            timeStampsTraining = [timeStampsTraining ;t];
        end;
        timeStampsTest = ([0:size(Ytest,1)-1].*dt)';
    case 'TrafficFlow1D'
        load(dataSetName);
        Ytr=Y;
        Yts=Ytest;
        Y=[Ytr;Yts];
        seq=size(Ytr,1);
        prevSeq = 1;
        timeStampsTraining = [];
        dt=0.05;
        for i=1:length(seq)
            t = ([0:(seq(i)-prevSeq)].*dt)';
            prevSeq = seq(i)+1;
            timeStampsTraining = [timeStampsTraining ;t];
        end;
        timeStampsTest = ([0:size(Ytest,1)-1].*dt)';
    case 'TrafficFlow3DReconstruct'
        load(dataSetName);
        Ytr=Y;
        Yts=Ytest;
        Y=[Ytr;Yts];
        seq=size(Ytr,1);
        prevSeq = 1;
        timeStampsTraining = [];
        dt=0.05;
        for i=1:length(seq)
            t = ([0:(seq(i)-prevSeq)].*dt)';
            prevSeq = seq(i)+1;
            timeStampsTraining = [timeStampsTraining ;t];
        end;
        timeStampsTest = ([0:size(Ytest,1)-1].*dt)';
        width=2;height=2;
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
    case 'toyMultigp2D'
        if experimentNo>10 && experimentNo<21 || experimentNo>60 && experimentNo<71 || experimentNo>40&& experimentNo<51
            load(['toyMultigp2Dregression' num2str(experimentNo-10)]);
        else
            load(['toyMultigp2DReconstruct' num2str(experimentNo)]);
        end
        Ytr=cell2mat(y);
        Yts=cell2mat(yTest);
        %         Yts=Yts(11:15,:);
        %         XTest=XTest(11:15);
        timeStampsTraining=X{1};
        timeStampsTest=XTest';
        Y=[Ytr;Yts];
        width=2;
        height=2;
    case 'sarcos_inv'
        load 'sarcos_inv'
        load sarcos_inv_test
        [Nall,Dall] = size(sarcos_inv);
        randn('seed', experimentNo);
        rand('seed', experimentNo);
        Ni = randperm(Nall);
        Ni = Ni(1:Ntrain);
        Ytr = sarcos_inv(Ni,end-6:end);
        X = sarcos_inv(Ni,1:end-7);
        Yts = sarcos_inv_test(:,end-6:end);
        Xts = sarcos_inv_test(:,1:end-7);
        Ytr = Ytr(:,chosen);
        Yts = Yts(:,chosen);
        Y = [Ytr;Yts];
        timeStampsTraining=X;
        timeStampsTest=Xts;
        % 中心化
        
        %         xmean = mean(timeStampsTraining);
        %         timeStampsTraining = timeStampsTraining- repmat(xmean,[size(timeStampsTraining,1),1]);
        %         timeStampsTest = timeStampsTest- repmat(xmean,[size(timeStampsTest,1),1]);
        
        % 标准化
        [timeStampsTraining,xmean,xstd] = standardize(timeStampsTraining,[],[]);
        timeStampsTest = standardize(timeStampsTest,xmean,xstd);
        %         pca_for_input
        %         [timeStampsTraining, timeStampsTest]=sarcos_all_pca(timeStampsTraining,timeStampsTest,6);
        % %         再标准化
        %         [timeStampsTraining,xmean,xstd] = standardize(timeStampsTraining,[],[]);
        %         timeStampsTest = standardize(timeStampsTest,xmean,xstd);
        clear sarcos_inv sarcos_inv_test
    case 'sarcos_inv_test'
        %          load 'sarcos_inv'
        load sarcos_inv_test
        [Nall,Dall] = size(sarcos_inv_test);
        randn('seed', experimentNo);
        rand('seed', experimentNo);
        Ni = randperm(Nall);
        Ni = Ni(1:Ntrain);
        Nis = setdiff([1:Nall],Ni);
        Ytr = sarcos_inv_test(Ni,end-6:end);
        X = sarcos_inv_test(Ni,1:end-7);
        Yts = sarcos_inv_test(Nis,end-6:end);
        Xts = sarcos_inv_test(Nis,1:end-7);
        Ytr = Ytr(:,chosen);
        Yts = Yts(:,chosen);
        Y = [Ytr;Yts];
        timeStampsTraining=X;
        timeStampsTest=Xts;
        % 中心化
        
        %         xmean = mean(timeStampsTraining);
        %         timeStampsTraining = timeStampsTraining- repmat(xmean,[size(timeStampsTraining,1),1]);
        %         timeStampsTest = timeStampsTest- repmat(xmean,[size(timeStampsTest,1),1]);
        
        % 标准化
        [timeStampsTraining,xmean,xstd] = standardize(timeStampsTraining,[],[]);
        timeStampsTest = standardize(timeStampsTest,xmean,xstd);
        %         pca_for_input
        %         [timeStampsTraining, timeStampsTest]=sarcos_all_pca(timeStampsTraining,timeStampsTest,6);
        % %         再标准化
        %         [timeStampsTraining,xmean,xstd] = standardize(timeStampsTraining,[],[]);
        %         timeStampsTest = standardize(timeStampsTest,xmean,xstd);
        clear sarcos_inv sarcos_inv_test
    case 'smallToyData'
        load('smallToyData');
        Ytr=cell2mat(y);
        Yts=cell2mat(yTest);
        timeStampsTraining=X{1};
        timeStampsTest=XTest';
        
        Y=[Ytr;Yts];
    case 'sinToy'
        try
            load('sin1Dregression');
            Ytr=cell2mat(y);
            Yts=cell2mat(yTest);
            timeStampsTraining=X{1};
            timeStampsTest=XTest';
            Y=[Ytr;Yts];
        catch
            error('no file');
        end
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
% add by ZhaoJing
if exist('timeStampsTraining') & exist('timeStampsTest')
    Nstar=size(Yts,1);
    t=sort([timeStampsTraining;timeStampsTest]);
else
    % 2013-10-09
    Ytr = Y(indTr,:); Yts = Y(indTs,:);
    t = linspace(0, 2*pi, size(Y, 1)+1)'; t = t(1:end-1, 1);
    timeStampsTraining = t(indTr,1); timeStampsTest = t(indTs,1);
end
% add by ZhaoJing
if strcmp(dataSetName,'toyMultigp2D')||strcmp(dataSetName,'smallToyData')||strcmp(dataSetName,'sinToy')
    Nstar=size(Yts,1);
    t=sort([timeStampsTraining;timeStampsTest]);
elseif strcmp(dataSetName,'sarcos_inv')||strcmp(dataSetName,'sarcos_inv_test')%2015-11-02
    Nstar=size(Yts,1);
    dataSetName = (['sarcos_inv' num2str(chosen)]);
    %     t=([timeStampsTraining;timeStampsTest]);
else
    % 2015-11-05
    Ytr = Y(indTr,:); Yts = Y(indTs,:);
    t = linspace(0, 2*pi, size(Y, 1)+1)';
    
    t = t(1:end-1, :);
    timeStampsTraining = t(indTr,:); timeStampsTest = t(indTs,:);
end

%-- For DEBUG
% testOnTrainingTimes=1;testOnTrainingData=1;
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
DgtN=0;
% 2013-06-08

if ~exist('DgtN') || ~DgtN
    options.enableDgtN = false;
end
d = size(Ytr, 2);


if trainModel
    % demo using the variational inference method for the gplvm model
    fprintf(1,'# Creating the model...\n');
    % scale = std(Ytr);
    % scale(find(scale==0)) = 1;
    % options.scaleVal = mean(std(Ytr));
    options.scaleVal = sqrt(var(Ytr));
    if fixInd
        options.fixInducing=1;
        options.fixIndices=1:size(Ytr,1);
    end
    %% added by JingZhao
    if latentDim>d
        options.initX = sarcos_pca(timeStampsTraining, [],1);
    else
        
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
        
        kern = kernCreate(timeStampsTraining, dynamicKern); % Default: {'rbf','white','bias'}
        
        
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
%     model.beta = 150;
    modelInit = model;
    
    % disp(model.vardist.covars)
    % disp(model.kern.comp{1}.inputScales)
    
    
    %%---
    capName = dataSetName;
    capName(1) = upper(capName(1));
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    fileToSave = ['dem' capName modelType num2str(experimentNo) ];
    %%---
    
    if exist('experimentDesc'), model.expDesc = experimentDesc; end
    
    display = 1;
    %     %%%% Optimisation
    %     % do not learn beta for few iterations for intitialization
    
    % modefied by ZhaoJing
    
    if exist('continueTrain')&& continueTrain % 继续上一次的迭代
        capName = dataSetName;
        capName(1) = upper(capName(1));
        %     modelType = model.type;
        modelType(1) = upper(modelType(1));
        fileToLoad = ['dem' capName modelType num2str(experimentNo)];
        load(fileToLoad);
%         fhandle=str2func([model.type 'RestorePrunedModel']);
%         model=fhandle(prunedModel, Ytr);
        itNo=[1000,1000];
        display = 1;
        
    else    % in 2013-06-07
        
        if fixedBetaIters ~=0
            model.learnBeta = 0;
            fprintf(1,'# Intitiliazing the model (fixed beta)...\n');
            model = vargplvmOptimise(model, display, fixedBetaIters); % Default: 20
            fprintf(1,'1/b = %.4d\n',1/model.beta);
            modelBetaFixed = model;
            model.fixedBetaIters=fixedBetaIters;
            model.learnBeta = 1;
        end
        
        if initVardistIters ~=0
            model.initVardist = 1;
            model.learnBeta = 0; model.learnSigmaf = 0; % This should be merged with the initVardist field
            fprintf(1,'# Intitiliazing the model (fixed beta)...\n');
            model = vargplvmOptimise(model, display, initVardistIters); % Default: 20
            fprintf(1,'1/b = %.4d\n',1/model.beta);
            model.learnSigmaf = 1;
            model.initVardist = 0;
            model.learnBeta = 1; 
        end
        if usingEM1~=0
            for i=1:alternateIters1
                model.alternateIters1=alternateIters1;
                model.onlyVariational = 1;
                fprintf(1,'%d #  E Step......\n', i);
                model = vargplvmOptimise(model, display, EIters1); % Default: 20
                model.onlyVariational = 0;
                model.onlyKernel = 1;
                fprintf(1,'%d #  M Step\n',i);
                model = vargplvmOptimise(model, display, MIters1); % Default: 20
                model.onlyKernel = 0;
            end
        end
        if usingEM2~=0
            for i=1:alternateIters2
                model.alternateIters2=alternateIters2;
                model.onlyVariational = 1;
                fprintf(1,'%d #  E Step......\n', i);
                model = vargplvmOptimise(model, display, EIters2); % Default: 20
                model.onlyVariational = 0;
                model.onlyKernel = 1;
                fprintf(1,'%d #  M Step\n',i);
                model = vargplvmOptimise(model, display, MIters2); % Default: 20
                model.onlyKernel = 0;
            end
        end
        
    end
    
    model.learnBeta = 1;
    model.iters = 0;
    
    %new!!!
    prunedModelInit = vargplvmPruneModel(modelInit);
    clear modelInit
    %__
    
    % Optimise the model.
    % modified by ZhaoJing
    alltime=0;
    tic;
    % 2013-06-08
    
    for i=1:length(itNo)
        iters = itNo(i); % default: 2000
        fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
        model = vargplvmOptimise(model, display, iters);
        model.iters = model.iters + iters;
        fprintf(1,'1/b = %.4d\n',1/model.beta);
        modelTr = model;
        
        fprintf(1,'# 1/b=%.4f\n var(m)=%.4f\n',1/model.beta, var(model.mOrig(:)));
        
        % Save model
        fprintf(1,'# Saving %s\n',fileToSave);
        prunedModel = vargplvmPruneModel(model);
        
        % prunedModelTr = vargplvmPruneModel(modelTr);
        save([fileToSave,'prune'], 'prunedModel', 'prunedModelInit');
        save(fileToSave, 'model');
    end
    
else
    % Load pre-trained model
    capName = dataSetName;
    capName(1) = upper(capName(1));
    %     modelType = model.type;
    modelType(1) = upper(modelType(1));
    fileToLoad = ['dem' capName modelType num2str(experimentNo)];
    load(fileToLoad);
    fileToSave = fileToLoad;
end


% Don't make predictions for the non-dynamic case
if ~dynUsed
    return
end



%%%-----------------------   RECONSTRUCTION --------------------%%%

%%
% load(['./N100/' fileToSave]);
if ~exist('model')
    fhandle=str2func([prunedModel.type 'RestorePrunedModel']);
    model=fhandle(prunedModel, Ytr);
end
% model = vargplvmRestorePrunedModel(prunedModel, Ytr);
% clear modelTr %%
% model.y = Ytr;
% model.m= gpComputeM(model); %%%
% model.y=[];
% load('model.mat');
Nstar = size(YtsOriginal,1);
%%% NOTE: If you are reloading and the restoring the "modelUpdated" you will need to:
%  modelUpdated.m = modelUpdated.mOrig;
%  modelUpdated.P = modelUpdated.P1 * (modelUpdated.Psi1' * modelUpdated.m);
%%

% Prediction using the only information in the test time points
[Testmeans2 Testcovars2] = vargplvmPredictPoint(model.dynamics, timeStampsTest);


%% calculate -log poseterior likelihood

% load('VGPDSToyindPosterior','pnew');
% pnew=[pnew;p];
% save('VGPDSToyindPosterior','pnew');
%
[Varmu2,varsigma] = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2);
% Mean absolute error per pixel
% errorOnlyTimes = mean(abs(Varmu2(:) - YtsOriginal(:)))/mean(abs(mean(Ytr(:)) - YtsOriginal(:)).^2);


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

% fprintf(1,'# Error GPLVM: %d\n', errorOnlyTimes);
%fprintf(1,'# Error NN: %d\n', errorNN);

for i=1:size(Varmu2,2)
    %     SMSEy=mean(abs(Varmu2(:,i) - YtsOriginal(:,i)).^2)/mean(abs(mean(YtsOriginal(:,i)) - YtsOriginal(:,i)).^2);
    %     fprintf(1,'# SMSE on y_%d:%d\n',i, SMSEy);
    RMSE(i)=sqrt(mean(abs(Varmu2(:,i) - YtsOriginal(:,i)).^2));
    fprintf(1,'# RMSE on y_%d:%d\n',i, RMSE(i));
    SMSEy(i)=mean(abs(Varmu2(:,i) - YtsOriginal(:,i)).^2)/mean(abs(mean(Ytr(:,i)) - YtsOriginal(:,i)).^2);
    fprintf(1,'# SMSE on y_%d:%d\n',i, SMSEy(i));
end
for i=1:size(YtsOriginal,2)
    p(i)=calculatePosteriorLikelihood(model,Testmeans2,Testcovars2,YtsOriginal,i);
end
% for i=1:size(Varmu2,2)
%     p(i)=calculatePosteriorLikelihood(model,Testmeans2,Testcovars2,YtsOriginal,i);
%     fprintf(1,'LP on y_%d=%.10f\n',i,p(i));
% end
% load ([fileToSave 'result']);

vars = mean(varsigma);
save([fileToSave 'result'],'Varmu2','YtsOriginal','RMSE','p','vars','SMSEy');
if size(Varmu2,2)<5
    % added by ZhaoJing
    
    scaleVal=1;
    if showFigure
        close all
        for i=1:size(Varmu2,2);
            figure
            hold on
            xlim=[timeStampsTest(1),timeStampsTest(end)];
            %         ymax=max()
            
            
            if exist('varsigma')&&~isempty(varsigma)
                f = [(Varmu2(:,i)+2*real(sqrt(varsigma(:,i))))*scaleVal;flipdim((Varmu2(:,i)-2*real(sqrt(varsigma(:,i))))*scaleVal,1)];
                a = fill([timeStampsTest; flipdim(timeStampsTest,1)], f, [6 6 6]/8, 'EdgeColor', [6 6 6]/8);
                
            end
            
            %         a =plot(timeStampsTest, Varmu2(:,i),'k-');
            c =plot(timeStampsTest,YtsOriginal(:,i),'r.');
            b= plot(timeStampsTest, Varmu2(:,i)*scaleVal,'b-');
            %         minimum = min((Varmu2(:,i)-2*real(sqrt(varsigma(:,i)))));
            %         maximum = max((Varmu2(:,i)+2*real(sqrt(varsigma(:,i)))));
            %         ylim = [min(minimum,min(YtsOriginal(:,i))) max(maximum,max(YtsOriginal(:,i)))];
            %         ylim=[0.65,1.05];
            xlabel('Time-t','fontsize',20);
            ylabel('Value-y','fontsize',20);
            %         title('VDM-GPDS','fontsize',12);
            l=legend([b c],'Pred','True','Location','SouthEast');
            set(l,'Fontsize',20);
            %         ylim=[min(YtsOriginal(:,i))-1,max(YtsOriginal(:,i))+0.1];
            set(c,   'markersize', 14);
            set(b,   'lineWidth', 2);
            set(gca, 'fontname', 'arial', 'fontsize', 20, 'box', 'on','xlim', xlim,'ylim',ylim);
        end
    end
    % 2013-10-23
end
% % Visualization of the reconstruction
%{
% for i=1:Nstar
%     subplot(1,2,1);
%     fr=reshape(YtsOriginal(i,:),height,width);
%     imagesc(fr);
%     colormap('gray');
%     subplot(1,2,2);
%     fr=reshape(Varmu2(i,:),height,width);
%     imagesc(fr);
%     colormap('gray');
%     pause(0.2);
% end
%}




%%%------------------ Extra --------------------%%%%%%%

% The full movie!!
%{
% for i=1:size(Varmu2,1)
%     fr=reshape(Ytr(i,:),height,width);
%     imagesc(fr);
%     colormap('gray');
%     pause(0.09)
%     fr=reshape(Varmu2(i,:),height,width);
%     imagesc(fr);
%     colormap('gray');
%     pause(0.001)
% end
%}


%%


% if you like you can also do prediction after training with partially
% observed test images (also good for de-bugging because this should be better
% than the onlyTimes prediction)



if predWithMs
    %
    display = 1;
    indexP = [];
    Init = [];
    
    fprintf(1, '# Partial reconstruction of test points...\n');
    % w=360; % width
    % h=288; % height
    w=width; h=height;
    
    if ~exist('cut')
        cut = 'vertically';
        if strcmp(dataSetName,'missa') || strcmp(dataSetName,'ADN') || strcmp(dataSetName,'claire')
            cut = 'horizontally';
        elseif strcmp(dataSetName,'sdata')|| strcmp(dataSetName,'toyMultigp2D') ||strcmp(dataSetName,'TrafficFlow4DReconstruct')
            cut='fix';
        end
        
        
    end
    
    switch cut
        case 'horizontally'
            if strcmp(dataSetName,'ADN')
                cutPoint=round(h/1.55);
            elseif strcmp(dataSetName,'claire')
                cutPoint=round(h/2);
            elseif strcmp(dataSetName, 'grandma')
                cutPoint=round(h/1.68);
            else %missa
                cutPoint=round(h/2)+13;
            end
            if exist('cutP')    cutPoint = cutP; end
            mask = [ones(1,cutPoint) zeros(1,h-cutPoint)];
            mask=repmat(mask, 1,w);
            indexMissing = find(mask);
        case 'vertically'
            if strcmp(dataSetName,'missa')
                indexMissing=1:round(size(Yts,2)/1.8);
            elseif strcmp(dataSetName,'dog')
                indexMissing = 1:round(size(Yts,2)/1.70);
            elseif strcmp(dataSetName,'ADN')
                indexMissing = 1:round(size(Yts,2)/1.7);
                
            elseif strcmp(dataSetName,'ocean')
                indexMissing = 1:round(size(Yts,2)/1.6);
            elseif strcmp(dataSetName,'head')
                indexMissing=1:round(size(Yts,2)/2.08);
            else
                indexMissing=1:round(size(Yts,2)/2);
            end
        case 'fix'
            
            
            %             YtsOriginal=YtsOriginal(1:50,:);
    end
    indexPresent = setdiff(1:model.d, indexMissing);
    %     YtsOriginal=YtsOriginal(1:10,:);
    %     timeStampsTest=timeStampsTest(1:10,:);
    %
    
    Yts=YtsOriginal;
    timeStampsTestOriginal=timeStampsTest;
    %     YtsBack = YtsOriginal;
    %     timeStampsTestBack=timeStampsTest;
    % See the movie after cutting some dims
    %{
%     Ytemp=Yts;
%     Ytemp(:,indexMissing)=NaN;
%     for i=1:size(Ytemp,1)
%         fr=reshape(Ytemp(i,:),h,w);
%         imagesc(fr);
%         colormap('gray');
%         pause(0.1)
%     end
%     clear Ytemp
    %}
    %     for index=1:30
    %     Yts=YtsBack(index,:);
    %      YtsOriginal=YtsBack(index,:);
    %      timeStampsTest=timeStampsTestBack(index,:);
    % testIndex={[1,2,3,4,5];[6,7,8,9,10];[11,12,13,14,15]};
    testIndex=[1:15];
    Varmu=[];
    
    
    if trainModelReconstruct
        for cycle=1:size(testIndex,1)
            Yts=YtsOriginal(testIndex(cycle,:),:);
            timeStampsTest=timeStampsTestOriginal(testIndex(cycle,:),:);
            mini =[];
            for i=1:size(Yts,1)
                % initialize the latent points using the nearest neighbour
                % from he training data
                dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
                % mind: smaller distance %mini(i): smaller index
                [mind, mini(i)] = min(dst);
            end
            
            %{
        % new barmu
        %jit = 1e-6;
        %%muInit = [model.vardist.means; model.vardist.means(mini,:)];
        %Kt = kernCompute(model.dynamics.kern, timeStampsTest);
        %Lkt = chol(Kt + jit*mean(diag(Kt))*eye(size(Kt,1)))';
        %%barmuInit = Lkt'\(Lkt\muInit);
        %%model.dynamics.vardist.means = barmuInit(1:model.N,:);
        %vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
        %vardistx.means = Lkt'\(Lkt\model.vardist.means(mini,:)) + 0.5*randn(size(vardistx.means));
        %vardistx.covars = 0.2*ones(size(vardistx.covars));
            %}
            vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
            vardistx.covars = model.dynamics.vardist.covars(mini,:);%0.2*ones(size(vardistx.covars));
            model.vardistx = vardistx;
            
            Yts(:,indexMissing) = NaN;
            
            model.dynamics.t_star = timeStampsTest;
            iters=reconstrIters;
            
            %-- NEW
            model.mini = mini;
            clear mini %%%
            model.dynamics.reoptimise=0;
            model.dynamics.onlyTest=0;
            %---
            [x, varx, modelUpdated] = vargplvmOptimisePoint(model, vardistx, Yts, display, iters);
            barmu = x;
            lambda = varx;
            [x, varx, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(modelUpdated, barmu, lambda, Yts);
            Testmeans = x;
            Testcovars = varx;
            [VarmuEvery sigma] = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
            Varmu=[Varmu;VarmuEvery];
            
            %     fhandle=str2func([prunedModelUpdated.type 'RestorePrunedModel']);
            %     modelUpdated=fhandle(prunedModelUpdated, Ytr);
            
        end
        % Find the absolute error
        save([fileToSave 'Pred' num2str(constructNo)], 'Testmeans', 'Testcovars', 'modelUpdated','Varmu');
        fprintf(1,'# Saved %s\n',[fileToSave 'Pred' num2str(constructNo)]);
        
    else
        
        load([fileToSave 'Pred' num2str(constructNo)]);% Mean error per pixel
        %     Varmu = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
    end
    %     SMSEy=mean(abs(Varmu2(:,i) - YtsOriginal(:,i)).^2)/mean(abs(mean(YtsOriginal(:,i)) - YtsOriginal(:,i)).^2);
    %     fprintf(1,'# GPLVM Error (in the missing dims) with missing inputs:%d\n', errorFull);
    p=calculatePosteriorLikelihood(modelUpdated,Testmeans,Testcovars,YtsOriginal,indexMissing);
    save(['ptemp' num2str(experimentNo)],'p');
    [Varmu varsigma] = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
    
    %     Varmu=[Varmu;VarmuEvery];
    errorFull = sum(sum( abs(Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    %     Norm=sum(sum( abs(mean(YtsOriginal(:,indexMissing)) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    %     SMSE=errorFull/Norm;
    fprintf(1,'# GPLVM Error (in the missing dims) with missing inputs:%d\n', errorFull);
    
    errOnlyTimes2 = sum(sum( abs(Varmu2(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# GPLVM Error (in the missing dims) only times:%d\n', errOnlyTimes2);
    
    errorFullPr = sum(sum( abs(Varmu(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
    fprintf(1,'# GPLVM Error (in the present dims) with missing inputs:%d\n', errorFullPr);
    
    errorFullPr2 = sum(sum( abs(Varmu2(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
    fprintf(1,'# GPLVM Error (in the present dims) only times:%d\n', errorFullPr2);
    RMSE=sqrt(mean((Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)).^2));
    fprintf(1,'# RMSE on y_%d :%.10f\n', indexMissing,RMSE);
    %-------- NN  ----------
    fprintf(1,'# NN prediction...\n');
    mini =[];
    sortedInd = zeros(size(Yts,1), size(Ytr,1));
    for i=1:size(Yts,1)
        dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
        % mind: smaller distance %mini(i): smaller index
        %[mind, mini(i)] = min(dst);
        [void ,sortedInd(i,:)] = sort(dst,2);
    end
    clear dst %%%
    
    NNmuPart = zeros(Nstar, size(Ytr,2));
    % Set to 1 to find the two NN's in time space. set to 0 to find in
    % data space. timeNN=1 should be set only for datasets created with
    % the "everyTwo" option for the split.
    timeNN=0;
    k=2; % the k parameter in k-NN
    for i=1:Nstar
        if timeNN
            if i < Nstar
                NNmuPart(i, indexMissing) = 0.5*(Ytr(i,indexMissing) + Ytr(i+1,indexMissing));
            else
                NNmuPart(i, indexMissing) = Ytr(end,indexMissing);
            end
        else
            NNmuPart(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing);
            for n=2:k
                NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
            end
            NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)*(1/k);
        end
        NNmuPart(i, indexPresent) = Yts(i, indexPresent);
    end
    
    % Mean absolute error per pixel
    %errorNNPart = mean(abs(NNmuPart(:) - YtsOriginal(:)));
    errorNNPart = sqrt(mean((NNmuPart(:,indexMissing) - YtsOriginal(:,indexMissing)).^2)) ;
    fprintf(1,'# NN(2) Error (in the missing dims) with missing inputs:%d\n', errorNNPart);
    
    
    %%%%%%  Try different k for k-NN to see which is best
    if ~timeNN
        bestK=2; bestErr = errorNNPart; NNmuPartBest = NNmuPart;
        clear NNmuPart %%%
        for k=[1 3 4 5]; % the k parameter in k-NN, try different values
            NNmuPartK = zeros(Nstar, size(Ytr,2));
            for i=1:Nstar
                NNmuPartK(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing); % first NN
                for n=2:k % more NN's, if k>1
                    NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
                end
                NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)*(1/k); % normalise with the number of NN's
                
                NNmuPartK(i, indexPresent) = Yts(i, indexPresent);
            end
            errorNNPartK = sqrt(mean((NNmuPartK(:,indexMissing) - YtsOriginal(:,indexMissing)).^2));
            if errorNNPartK < bestErr
                bestErr = errorNNPartK;
                bestK=k;
                NNmuPartBest = NNmuPartK;
            end
        end
        clear NNmuPartK
    end
    % Mean absolute error per pixel
    fprintf(1,'# NNbest(%d) Error (in the missing dims) with missing inputs:%d\n',bestK, bestErr);
    %%%%%%%%%%%%%
    
    
    %     for i=1:Nstar
    %         subplot(1,2,1);
    %         fr=reshape(YtsOriginal(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         subplot(1,2,2);
    %         fr=reshape(NNmuPart(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         pause(0.5);
    %     end
    close all;
    figure
    hold on
    xlim=[timeStampsTest(1),timeStampsTest(end)];
    %         ymax=max()
    %         ylim=[1.8,3.4];
    i=indexMissing;
    if exist('varsigma')&&~isempty(varsigma)
        f = [(Varmu(:,i)+2*real(sqrt(varsigma(:,i))))*scaleVal;flipdim((Varmu(:,i)-2*real(sqrt(varsigma(:,i))))*scaleVal,1)];
        a = fill([timeStampsTest; flipdim(timeStampsTest,1)], f, [6 6 6]/8, 'EdgeColor', [6 6 6]/8);
        b= plot(timeStampsTest, Varmu(:,i)*scaleVal,'b-');
    end
    %         a =plot(timeStampsTest, Varmu2(:,i),'k-');
    c =plot(timeStampsTest,YtsOriginal(:,i),'r.');
    minimum = min((Varmu(:,i)-2*real(sqrt(varsigma(:,i)))));
    maximum = max((Varmu(:,i)+2*real(sqrt(varsigma(:,i)))));
    % if isfield(model, 'X_u') && ~isempty(model.X_u);
    %   b = plot(model.X_u, -1.5, 'kx');
    %         set(b, 'linewidth', 2)
    %         set(b, 'markersize', 10);
    % end
    
    ylim = [min(minimum,min(YtsOriginal(:,i))) max(maximum,max(YtsOriginal(:,i)))];
    ylim=[2.7,3.3];
    xlabel('Time-t','fontsize',20);
    ylabel('Value-y','fontsize',20);
    %         title('VDM-GPDS','fontsize',12);
    l=legend([b,c], 'Pred','True','Location','SouthEast');
    set(l,'fontsize',20);
    %         ylim=[min(YtsOriginal(:,i))-1,max(YtsOriginal(:,i))+0.1];
    set(c,   'markersize', 24);
    set([a b],   'lineWidth', 2);
    set(gca, 'fontname', 'arial', 'fontsize', 20, 'box', 'on','xlim', xlim,'ylim',ylim);
    figure;
    % h1=plot(Varmu(:,indexMissing),'r-*');
    % hold on
    % % h2=plot(Varmu2(:,indexMissing));
    % % hold on
    b=plot(timeStampsTest,NNmuPartBest(:,indexMissing),'b-');
    hold on
    c =plot(timeStampsTest,YtsOriginal(:,i),'r.');
    ylim=[2.7,3.3];
    xlabel('Time-t','fontsize',20);
    ylabel('Value-y','fontsize',20);
    %         title('VDM-GPDS','fontsize',12);
    l=legend([b,c], 'Pred','True','Location','SouthEast');
    set(l,'fontsize',20);
    %         ylim=[min(YtsOriginal(:,i))-1,max(YtsOriginal(:,i))+0.1];
    set(c,   'markersize', 24);
    set(b,   'lineWidth', 2);
    set(gca, 'fontname', 'arial', 'fontsize', 20, 'box', 'on','xlim', xlim,'ylim',ylim);
    % % hold on
    % h4=plot(YtsOriginal(:,indexMissing),'k+');
    % legend([h1,h3,h4]','Varmu','NNmuPartBest','YtsOriginal');
    % The following causes OUTOFMEMORY exception
    % lvmVisualise(model, [], 'imageVisualise', 'imageModify', [720 1280],1,0,1);
    save([fileToSave 'Pred' num2str(constructNo) 'result'],'Varmu','YtsOriginal','indexMissing','bestErr');
end



