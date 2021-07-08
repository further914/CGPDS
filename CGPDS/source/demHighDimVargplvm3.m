% DEMHIGHDIMVARGPLVM3 This is the main script to run CGPDS on
% high dimenisional video datasets.
%
% COPYRIGHT: Andreas C. Damianou, Michalis K. Titsias, 2010 - 2011


% Fix seeds
randn('seed',1e5);
rand('seed',1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 30;      end
if ~exist('itNo')         ,  itNo =2000;              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 49;          end     % Default: 49
if ~exist('latentDim')    ,  latentDim = 25;          end
% Set to 1 to use dynamics
if ~exist('dynUsed')      ,  dynUsed = 1;             end
% Set to 1 to keep only the dimensions modelling the head (missa dataset)
% if ~exist('onlyVariationalParam'), optiLocalParam = 200;      end     % DEFAULT: 23
% if ~exist('onlyModelParam'), optiGlobalParam = 200;      end     % DEFAULT: 23
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('mappingKern')   ,     mappingKern = {'rbfard2', 'bias', 'white'}; end
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
if ~ exist('enableParallelism'), enableParallelism = 1; end
% Set it to 1 to retrain the model. Set it to 0 to load an already trained
% one.
if ~exist('isTraining'), isTraining = 1; end

fprintf(1,'\n#----------------------------------------------------\n');
fprintf(1,'# Dataset: %s\n',dataSetName);
fprintf(1,'# ExperimentNo: %f\n', experimentNo);
fprintf(1,'# Latent dimensions: %f\n',latentDim);
fprintf(1,'# Reconstruction iterations: %f\n',reconstrIters);
fprintf(1,'# InitX: %s\n',initX);
fprintf(1,'# Reoptimise inducing points (for reconstr): %f \n',testReoptimise);
fprintf(1,'# Dynamics used: %f\n', dynUsed);
if dynUsed
    fprintf(1,'# Dynamics kern: ');
    disp(dynamicKern);
end
fprintf(1,'# VardistCovarsMult: %f \n', vardistCovarsMult);
fprintf(1,'# InvWidthMultDyn: %f \n', invWidthMultDyn);
fprintf(1,'# InvWidthMult: %f \n', invWidthMult);
if exist('dataToKeep')  fprintf(1,'# DataToKeep: %f \n',dataToKeep); end

switch dataSetName
    case 'missa' % There is a translation between frames 65...102
        Y = vargplvmLoadData('missa');
        height = 288; width = 360;
    case 'dog' % There is a translation between frames 65...102
        Y = vargplvmLoadData('dog');
        height = 480; width = 854;
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
        mask = [mask ones(1,lastBlockSize)]; % always end with the tr. set
        indTr = find(mask);
        indTs = find(~mask);
        if exist('msk')
            indTr = sort([indTr msk]);
            indTs=setdiff(1:N,indTr);
        end
        Nstar = size(indTs,2);
        Ntr = length(indTr);
        Nts = length(indTs);
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
        Ntr = length(indTr);
        Nts = length(indTs);
    case 'custom' %indTr and indTs must be provided
        Ntr = length(indTr);
        Nts = length(indTs);
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

fprintf(1,'# Inducing points: %f\n',indPoints);
fprintf(1,'# Dataset size used (train/test) : %f / %f \n', size(Ytr,1), size(Yts,1));
fprintf(1,'# Dataset Split: %s ',dataSetSplit);
if strcmp(dataSetSplit,'blocks'),    fprintf(1,' (blockSize:%f)',blockSize); end
fprintf(1,'\n');
fprintf(1,'#----------------------------------------------------\n');

clear Y % Free up some memory

%{
 % Play movie
 for i=1:size(Ytr,1)
     fr=reshape(Ytr(i,:),height,width);
     imagesc(fr); colormap('gray'); title(num2str(indTr(i)))
     pause%(0.08);
 end
 pause
 for i=1:size(Yts,1)
     fr=reshape(Yts(i,:),height,width);
     imagesc(fr); colormap('gray'); title(num2str(indTs(i)))
     pause%(0.08);
 end
%}

%%

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
options.numActive = indPoints;
if ~isempty(which('scg2'))
    options.optimiser = 'scg2';
else
%     options.optimiser = 'scg';
    options.optimiser = 'Adam';
end
if ~exist('DgtN') || ~DgtN
    options.enableDgtN = false;
end
if exist('scale2var1')
    options.scale2var1 = scale2var1;
end
d = size(Ytr, 2);

capName = dataSetName;
capName(1) = upper(capName(1));
modelType = 'cgpds';
ModelFile = [pwd filesep 'CGPDSResultsMat' filesep capName 'Results\dem' capName modelType num2str(experimentNo) '.mat'];
PredFile = [pwd filesep 'CGPDSResultsMat' filesep capName 'Results\dem' capName modelType num2str(experimentNo) 'Pred.mat'];
ResultFile = [pwd filesep 'CGPDSResultsMat' filesep capName 'Results\dem' capName modelType num2str(experimentNo) 'Results.mat' ];

try
    load(ModelFile);
    fprintf(1,'# Loading pre-trained model...\n');
catch
    % demo using the variational inference method for the cgpds
    fprintf(1,'# Creating the model...\n');
    % scale = std(Ytr);
    % scale(find(scale==0)) = 1;
    %options.scaleVal = mean(std(Ytr));
    if ~exist('scale2var1')
        if ~exist('scaleVal')
            options.scaleVal = sqrt(var(Ytr(:)));
            options.scale2var1 = 1;
        else
            options.scaleVal = scaleVal;
        end
    end
    
    J = 1;
    Rtr = 1;% the number of seq in trainning data
    Rts = 1;
    Ctr = ceil(Ntr/Rtr);%the length of each seq in trainning data
    Cts = ceil(Nts/Rts);
    seqTrain = 0;
    for r = 1:Rtr
        seqTrain = [seqTrain min(r*Ctr,Ntr)];
    end
    
    seqTest = 0;
    for r = 1:Rts
        seqTest = [seqTest min(r*Cts,Nts)];
    end
    model = cgpdsCreate(latentDim, d, J, Ytr, options);
    model.Ntr = Ntr;
    model.Nts = Nts;
    model.Rtr = Rtr;
    model.Rts = Rts;
    model.Ctr = Ctr;
    model.Cts = Cts;
    model.seqTrain = seqTrain;
    model.seqTest = seqTest;
    
    % Temporary: in this demo there should always exist the mOrig field
    if ~isfield(model, 'mOrig')
        model.mOrig = model.m;
    end
    
    dims = [1:model.d];
    seqs = [1:model.Rtr];
    model = cgpdsUpdateStats(model, model.X_v ,model.X_u, dims, seqs);
    model.kern_v.comp{1}.inputScales = invWidthMult./(((max(model.X)-min(model.X))).^2); % Default 5
    model.kern_u.comp{1}.inputScales = invWidthMult./(((max(model.X)-min(model.X))).^2); % Default 5
    params = cgpdsExtractParam(model,dims,seqs);
    model = cgpdsExpandParam(model, params,dims,seqs);
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
        model = cgpdsAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq,dims,seqs);
        
        fprintf(1,'# Further calibration of the initial values...\n');
        model = cgpdsInitDynamics(model,optionsDyn,dims,seqs);
        
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
end
if exist('experimentDesc'), model.expDesc = experimentDesc; end
    
% train the cgpds using adam algorithm
if isTraining  
    display = 1;
    stepSize = 1e-1;
    model.iters = 0;
    for i=1:length(itNo)
        iters = itNo(i); % default: 2000
        fprintf(1,'\n# Optimising all parameters for %f iterations (session %f)...\n',iters,i);
        model.optiLocalParam = 0;
        model.optiGlobalParam = 0;
        dims = [1:model.d];
        seqs = [1:1];
        model.optimiser = 'Adam';
        model.fixW = 1;
        model = cgpdsOptimise(model, display, iters, stepSize, dims ,seqs); % Default: 20
%         save(ModelFile, 'model');
        model.iters = model.iters + iters;
    end
else
    % Load pre-trained model
    disp(['# Loading pre-trained model number ' num2str(experimentNo) '...']);
    load(ModelFile);
    if exist('testReoptimise'), model.dynamics.reoptimise = testReoptimise; end
end

% Don't make predictions for the non-dynamic case
if ~dynUsed
    return
end



%%%-----------------------   RECONSTRUCTION --------------------%%%

% %%
% model = modelTr;
% clear modelTr %%
% model.y = Ytr;
% if isfield(model, 'mOrig')
%     model.m = model.mOrig;
% else
%     model.m= gpComputeM(model); %%%
% end
% model.y=[];
% Nstar = size(YtsOriginal,1);


% Prediction using the only information in the test time points
fprintf(1, '# Predicting only with the test time points...\n');
[Testmeans2 Testcovars2] = cgpdsPredictPoint(model.dynamics, timeStampsTest);
[Varmu2 Varcovars2] = cgpdsPosteriorMeanVar(model, Testmeans2, Testcovars2);

Varmu2 = Varmu2 .* repmat(model.scale,size(Yts,1),1);
Varmu2 = Varmu2 + repmat(model.bias,size(Yts,1),1);
Varcovars2 = Varcovars2.*repmat(model.scale.*model.scale,size(Yts,1), 1);

% % Mean absolute error per pixel
errorOnlyTimes = mean(abs(Varmu2(:) - YtsOriginal(:)));

%%%%% CHANGE THE NN PREDICTION, 
% 2-nearest (in time) prediction
NNmu = zeros(Nts, model.d);
for i=1:Nts
    if i < Nts
        NNmu(i, :) = 0.5*(Ytr(i,:) + Ytr(i+1,:));
    else
        NNmu(i, :) = Ytr(end,:);
    end
end
% Mean absolute error per pixel
errorNN = mean(abs(NNmu(:) - YtsOriginal(:)));
errorNN = sum(sum( abs(NNmu - YtsOriginal)) )/prod(size(YtsOriginal)); %equivalent

fprintf(1,'# Error CGPDS only time: %f\n', errorOnlyTimes);
fprintf(1,'# Error NN: %f\n', errorNN);

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
if ~exist('predWithMs'), predWithMs = 1; end

if predWithMs
    %%
    display = 1;
    indexP = [];
    Init = [];
    fprintf(1, '# Partial reconstruction of test points...\n');
    demHighDimPrepareTestData
    
    mini =[];
    for i=1:size(Yts,1)
        % initialize the latent points using the nearest neighbour
        % from he training data
        dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
        % mind: smaller distance %mini(i): smaller index
        [mind, mini(i)] = min(dst);
    end
    
    vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = model.dynamics.vardist.covars(mini,:);%0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
    Yts(:,indexMissing) = NaN;
    model.dynamics.t_star = timeStampsTest;
    iters=reconstrIters;
    model.mini = mini;
    clear mini %%%
    
    try
        load(PredFile);
    catch
        samd = [1:model.d];
        samr = [1:model.Rtr];
        model.optimiser = 'Adam';
        [x, varx] = cgpdsOptimiseSeqDyn(model, vardistx, Yts, 1, reconstrIters,samd,samr);
        % keep the optimized variational parameters
        barmu = x;
        lambda = varx;
        % Get the variational means and variacnes for the new test sequcen and
        % update the model to be prepared for prediction
        [x, varx, modelUpdated] = cgpdsDynamicsUpdateModelTestVar(model, barmu, lambda, Yts);
        Testmeans = x;
        Testcovars = varx;
    end
    
    [Varmu, Varcovars] = cgpdsPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
    Varmu = Varmu .* repmat(model.scale,size(Yts,1),1);
    Varmu = Varmu + repmat(model.bias,size(Yts,1),1);
    Varcovars = Varcovars.*repmat(model.scale.*model.scale,size(Yts,1), 1);
    
    ErrorPartObserMissingDim = sum(sum( abs(Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# CGPDS Error (in the missing dims) with missing inputs:%f\n', ErrorPartObserMissingDim);
    
    ErrorPartObserPresentDim = sum(sum( abs(Varmu(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
    fprintf(1,'# CGPDS Error (in the present dims) with missing inputs:%f\n', ErrorPartObserPresentDim);
    
    ErrorOnlyTimes2MissingDim = sum(sum( abs(Varmu2(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# CGPDS Error (in the missing dims) only times:%f\n', ErrorOnlyTimes2MissingDim);
    
    ErrorOnlyTimes2PresentDim = sum(sum( abs(Varmu2(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
    fprintf(1,'# CGPDS Error (in the present dims) only times:%f\n', ErrorOnlyTimes2PresentDim);
    
    likehood = 0;
    for n = 1:size(YtsOriginal,1)
        likehood = likehood - model.d/2*log(2*pi)-0.5*sum(log((Varcovars2(n,:))))-0.5*sum((YtsOriginal(n,:)-Varmu2(n,:)).*...
            (1./Varcovars2(n,:)).*(YtsOriginal(n,:)-Varmu2(n,:)));
    end
    MSLLOnlyTimeCGPDS = - likehood/prod(size(YtsOriginal));
    fprintf(1,'# CGPDS MSLL only times:%f\n', MSLLOnlyTimeCGPDS);
    
    likehood = 0;
    for n = 1:size(YtsOriginal,1)
        likehood = likehood - length(indexMissing)/2*log(2*pi)-0.5*sum(log((Varcovars(n,indexMissing))))-0.5*sum((YtsOriginal(n,indexMissing)-Varmu(n,indexMissing)).*...
            (1./Varcovars(n,indexMissing)).*(YtsOriginal(n,indexMissing)-Varmu(n,indexMissing)));
    end
    MSLLPartObserMissingDimCGPDS = - likehood/(size(YtsOriginal,1)*length(indexMissing));
    fprintf(1,'# CGPDS MSLL (in the missing dims) with missing inputs:%f\n', MSLLPartObserMissingDimCGPDS);
    
    %                 save(ResultFile,'Testmeans','Testcovars','Varmu','Varcovars',...
    %                     'ErrorPartObserMissingDim','ErrorOnlyTimes2MissingDim','ErrorPartObserPresentDim','ErrorOnlyTimes2PresentDim',...
    %                     'Testmeans2','Testcovars2','Varmu2','Varcovars2','errorOnlyTimes','-append');

            
    
    



   
     % Visualization of the reconstruction
%          for i=1:Nstar
%              subplot(1,2,1);
%              fr=reshape(YtsOriginal(i,:),height,width);
%              imagesc(fr);
%              colormap('gray');
%              subplot(1,2,2);
%              fr=reshape(Varmu(i,:),height,width);
%              imagesc(fr);
%              colormap('gray');
%              pause(0.5);
%          end
%          fr=reshape(Varmu(46,:),height,width); imagesc(fr); colormap('gray'); % VGPDS


    %% -------- NN  ----------
    meanFrame = repmat(mean(Ytr),size(YtsOriginal,1),1);
    errMeanFrame = sum(sum(abs(meanFrame(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# Mean Frame Error (in the missing dims) with missing inputs: %f\n', errMeanFrame);

%     fprintf('# Sequential NN...\n')
%     try
%         [NNpredSeq, errorsNNSeq] = NNseq(Ytr, YtsOriginal, Yts, 1:6);
%         [minEr, indMinEr] = min(errorsNNSeq);
%         fprintf(1,'# NN Seq.(%f) Error (in the missing dims) with missing inputs:%f\n', indMinEr, minEr);
%     catch e
%     end
%     
%     fprintf(1,'# Classic NN prediction...\n');
% 
%     mini =[];
%     sortedInd = zeros(size(Yts,1), size(Ytr,1));
%     for i=1:size(Yts,1)
%         dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
%         % mind: smaller distance %mini(i): smaller index
%         %[mind, mini(i)] = min(dst);
%         [void ,sortedInd(i,:)] = sort(dst,2);
%     end
%     clear dst %%%
%     
%     NNmuPart = zeros(Nstar, size(Ytr,2));
%     % Set to 1 to find the two NN's in time space. set to 0 to find in
%     % data space. timeNN=1 should be set only for datasets created with
%     % the "everyTwo" option for the split.
%     timeNN=0;
%     k=2; % the k parameter in k-NN
%     for i=1:Nstar
%         if timeNN
%             if i < Nstar
%                 NNmuPart(i, indexMissing) = 0.5*(Ytr(i,indexMissing) + Ytr(i+1,indexMissing));
%             else
%                 NNmuPart(i, indexMissing) = Ytr(end,indexMissing);
%             end
%         else
%             NNmuPart(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing);
%             for n=2:k
%                 NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
%             end
%             NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)*(1/k);
%         end
%         NNmuPart(i, indexPresent) = Yts(i, indexPresent);
%     end
%     
%     % Mean absolute error per pixel
%     %errorNNPart = mean(abs(NNmuPart(:) - YtsOriginal(:)));
%     errorNNPart = sum(sum( abs(NNmuPart(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
%     fprintf(1,'# NN(2) Error (in the missing dims) with missing inputs:%f\n', errorNNPart);
%     
%     
%     %%%%%%  Try different k for k-NN to see which is best
%     if ~timeNN
%         bestK=2; bestNNErr = errorNNPart; NNmuPartBest = NNmuPart;
%         clear NNmuPart %%%
%         for k=[1 3 4 5]; % the k parameter in k-NN, try different values
%             NNmuPartK = zeros(Nstar, size(Ytr,2));
%             for i=1:Nstar
%                 NNmuPartK(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing); % first NN
%                 for n=2:k % more NN's, if k>1
%                     NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
%                 end
%                 NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)*(1/k); % normalise with the number of NN's
%                 
%                 NNmuPartK(i, indexPresent) = Yts(i, indexPresent);
%             end
%             errorNNPartK = sum(sum( abs(NNmuPartK(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
%             if errorNNPartK < bestNNErr
%                 bestNNErr = errorNNPartK;
%                 bestK=k;
%                 NNmuPartBest = NNmuPartK;
%             end
%         end
%         clear NNmuPartK
%     end
%     % Mean absolute error per pixel
%     fprintf(1,'# NNbest(%f) Error (in the missing dims) with missing inputs:%f\n',bestK, bestNNErr);
%     %%%%%%%%%%%%%
%     save(ResultFile,'errorNNPart','bestK','bestNNErr','-append');
    
    
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
end

%% Movies
% playMov(h,w,[],Varmu, Yts)
% playMov(h,w,[],Varmu,NNmuPartBest)

% The following causes OUTOFMEMORY exception
% lvmVisualise(model, [], 'imageVisualise', 'imageModify', [720 1280],1,0,1);

%%
%{
plot(modelUpdated.X(:,1), modelUpdated.X(:,2), 'x'); hold on; plot(modelUpdated.X_u(:,1), modelUpdated.X_u(:,2), 'ro'); legend('mUpdX','mUpdXu')
figure
plot(model.X(:,1), model.X(:,2), '+'); hold on; plot(model.X_u(:,1), model.X_u(:,2),'ro'); legend('modelX', 'modelX_u')
%}