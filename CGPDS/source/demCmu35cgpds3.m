% DEMCMU35CGPDS3 Run CGPDS on CMU35 data.

% Fix seeds
clear
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 1;      end
if ~exist('isTraining') ,  isTraining = 1;      end
if ~exist('indPoints')    ,  indPoints = 30;          end     % Default:
if ~exist('latentDim')    ,  latentDim = 9;          end
% Set to 1 to use dynamics
if ~exist('dynUsed')      ,  dynUsed = 1;             end
if ~exist('onlyVariationalParamIter'), onlyVariationalParamIter = 2000;      end
if ~exist('onlyModelParamIter'), onlyModelParamIter = 2000;      end
% if ~exist('onlyBeta'), onlyBeta = 50;      end
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('fixInd')        ,     fixInd = 0;                             end
if ~exist('baseKern')   ,     baseKern = {'rbfard2'}; end % rbfard2 and convolve to be chosen
if ~exist('modelType'),          modelType='cgpds'; end %varmgplvm and vargplvm to be chosen
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'bias', 'white'}; end
if ~exist('reconstrIters') ,     reconstrIters = 5000;                   end
if ~exist('learnVariance'),     learnVariance =0;   end
%if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
if ~exist('initX'),     initX ='ppca';   end
if ~exist('doReconstr'),     doReconstr=1;   end
if ~exist('onlyJog'),     onlyJog=1;   end
if ~exist('onlyWalk'),     onlyWalk=0;   end
if ~exist('continueTrain'),     continueTrain=1;   end
if ~exist('subset'),     subset=9;   end

% Get the sequence numbers.
[Y, lbls] = lvmLoadData('cmu35WalkJog');
seq = cumsum(sum(lbls)) - [1:31];

dataSetName = 'cmu35gplvm';

% load data
% 1-16 for walk
% 17-25 for jog
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);

if exist('onlyJog')&& onlyJog
    Y=Y(seq(16)+1:seq(25),:);
    seq=seq-seq(16);
    seq=seq(17:17+8);
elseif exist('onlyWalk') && onlyWalk
    Y=Y(1:seq(15),:);
%     seq=seq-seq(16);
    seq=seq(1:15);
end
refY=Y;
Y=refY(1:seq(subset),:);
sNum=size(seq,2);
for i=1:sNum-subset
    Ytests{i}=refY(seq(i-1+subset)+1:seq(i+subset),:);
end
% Ytests{i+1}=Ytest;
seq=seq(1:subset);

% save([dataSetName 'JogSub'],'Y','Ytests');
fprintf(1,'\n#----------------------------------------------------\n');
fprintf(1,'# Dataset: %s\n',dataSetName);
fprintf(1,'# ExperimentNo: %d\n', experimentNo);
fprintf(1,'# Inducing points: %d\n',indPoints);
fprintf(1,'# Latent dimensions: %d\n',latentDim);
fprintf(1,'# Reconstruction iterations: %d\n', reconstrIters);
fprintf(1,'# Dataset size used (train) : %d \n', size(Y,1));
fprintf(1,'# Tie Inducing points: %d\n',fixInd);
fprintf(1,'# InitX: %s\n',initX);
fprintf(1,'# Base kern: ');
disp(baseKern);
fprintf(1,'# Dynamics used: %d\n', dynUsed);
if dynUsed
    fprintf(1,'# Dynamics kern: ');
    disp(dynamicKern);
end

fprintf(1,'#----------------------------------------------------\n');
fileToSave = [pwd filesep 'CGPDSResultsMat' filesep 'cmu35Result' filesep 'demCmu35cgpds' num2str(experimentNo) '.mat'];
LastModel = [pwd filesep 'CGPDSResultsMat' filesep 'cmu35Result' filesep 'demCmu35cgpds' num2str(experimentNo-1) '.mat'];
% Fix times:
prevSeq = 1;
timeStampsTraining = [];
dt=0.05;
for i=1:length(seq)
    t = ([0:(seq(i)-prevSeq)].*dt)';
    prevSeq = seq(i)+1;
    timeStampsTraining = [timeStampsTraining ;t];
end
timeStampsTest = ([0:size(Ytest,1)-1].*dt)';

% Set up model
options = vargplvmOptions('dtcvar');
%options.kern = {'rbfard2', 'bias', 'white'};
options.kern = baseKern;
options.numActive = indPoints; % DEfault: 100

options.optimiser = 'scg';

d = size(Y, 2);
fprintf(1,'# Creating the model...\n');
if fixInd
    options.fixInducing=1;
    options.fixIndices=1:indPoints;
end

J=1;
dataSetName(1) = upper(dataSetName(1));

model = cgpdsCreate(latentDim, d, J, Y, options);
% samd and samr are used for stochastic optimization, but here we do not
% use stochastic optimization, so samd is the collection of all dimensions.
samd = [1:model.d];
samr = [1:1];

try
    load(LastModel);
    fprintf('load pre-trained model...\n');
catch
    model = cgpdsParamInit(model, model.m, model.X, samd, samr);
    model.beta=1/(0.01*var(model.m(:)));
    model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
    model.type = modelType;
end


%-------- Add dynamics to the model -----
if dynUsed
    optionsDyn.type = 'vargpTime';
    optionsDyn.t=timeStampsTraining;
    optionsDyn.inverseWidth=30;
    %   optionsDyn.vardistCovars = vardistCovarsMult;
    optionsDyn.seq = seq;
    optionsDyn.learnVariance = learnVariance;
    optionsDyn.initX = initX;
    optionsDyn.regularizeMeans = 0;
    
    % Dynamic kernel: the kernel of Ktt
    kern = kernCreate(t, dynamicKern); % Default: {'rbf','white','bias'}
    
    if strcmp(kern.comp{2}.type, 'white')
        kern.comp{2}.variance = 1e-2; % Usual values: 1e-1, 1e-3
    end
    
    if strcmp(kern.comp{2}.type, 'whitefixed')
        if ~exist('whiteVar')
            whiteVar = 1e-6;
        end
        kern.comp{2}.variance = whiteVar;
        fprintf(1,'# fixedwhite variance: %d\n',whiteVar);
    end
    
    if strcmp(kern.comp{1}.type, 'rbfperiodic')
        if exist('periodicPeriod')
            kern.comp{1}.period = periodicPeriod;
        end
        fprintf(1,'# periodic period: %d\n',kern.comp{1}.period);
    end
    
    % The following is related to the expected number of
    % zero-crossings.(larger inv.width numerator, rougher func)
    if ~strcmp(kern.comp{1}.type,'ou')
        kern.comp{1}.inverseWidth = optionsDyn.inverseWidth./(((max(t)-min(t))).^2);
        kern.comp{1}.variance = 1;
    end
    optionsDyn.kern = kern;
    
    if exist('vardistCovarsMult')
        optionsDyn.vardistCovars = vardistCovarsMult;
    end
    
    % Fill in with default values whatever is not already set
    optionsDyn = vargplvmOptionsDyn(optionsDyn);
    model = cgpdsAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq,samd,samr);
    
    fprintf(1,'# Further calibration of the initial parameters...\n');
    model = cgpdsInitDynamics11(model,optionsDyn,samd,samr);
end
modelInit = model;
model.reconstrIters = reconstrIters;


stepSize = 0.01;
if isTraining == 1
    % Optimise the model.
    display = 1;
    model.learnBeta = 1;
    model.iters=0;
    
    for i = 1:20
        if onlyVariationalParamIter > 0
            display = 1;
            fprintf(1,'# Optimizing the variational parameters %d iterations...\n',onlyVariationalParamIter);
            model.onlyVariationalParam = 1;
            model.onlyModelParam = 0;
            model = cgpdsOptimise(model, display, onlyVariationalParamIter,stepSize,samd,samr);
            disp('# Saving model after optimising variational parameters...');
            save(fileToSave,'model');
        end
        
        if onlyModelParamIter > 0
            display = 1;
            fprintf(1,'# Optimizing the model parameters %d iterations...\n',onlyModelParamIter);
            model.onlyModelParam = 1;
            model.onlyVariationalParam = 0;
            model = cgpdsOptimise(model, display, onlyModelParamIter,stepSize,samd,samr);
            disp('# Saving model after optimising model parameters...');
            save(fileToSave,'model');
        end
    end
else
    load(fileToSave);
end

% Reconstruction can also be done separately calling demCmu35cgpdsReconstructTaylor
if doReconstr
    %---- RECONSTRUCTION ---%
    fprintf('# Taylor Reconstruction for expNo:%d\n', experimentNo);
    demCmu35cgpdsReconstructTaylor
end
