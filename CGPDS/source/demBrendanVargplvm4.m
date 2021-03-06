
% DEMBRENDANVARGPLVM4 Run variational GPLVM on Brendan face data with missing (test) data.
%
%	Description:
%	Missing data are indicated with a NaN value in the corresponding cell of
%	the test matrix.
% 	demBrendanVargplvm4.m SVN version 1252
% 	last update 2011-03-04T09:54:39.000000Z

% Fix seeds
randn('seed', 1e6);
rand('seed', 1e6);

dataSetName = 'brendan';
experimentNo = 4;
printDiagram = 1;

fprintf(1, '# Preparing the dataset...\n');

% load data
[Y, lbls] = lvmLoadData(dataSetName);

dims = size(Y,2);

% For a subset of the data
N=size(Y,1); % Set N=size(Y,1) for the full dataset
Y = Y(1:N,:); 

% training and test sets
Ntr = round(size(Y,1)/2);

% Create Permutation matrix to split randomly training and test data
dataPerm = randperm(N); 
perm = dataPerm;

Ytr = Y(perm(1:Ntr),:);      %lblsTr = lbls(perm(1:Ntr),:);
Yts = Y(perm(Ntr+1:end),:);  %lblsTs = lbls(perm(Ntr+1:end),:);

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'white'};
options.numActive = 50; % Default 50
options.scale2var1 = 1; % scale data to have variance 1
%options.tieParam = 'tied';  

options.optimiser = 'scg';
latentDim = 30; % Default 30
d = size(Y, 2);

fprintf(1, '# Creating the model...\n');
model = vargplvmCreate(latentDim, d, Ytr, options);
%
model = vargplvmParamInit(model, model.m, model.X); 

iters = 500; % default 1500
display = 1;

fprintf(1, '# Training the model (%d iterations)...\n',iters);
model = vargplvmOptimise(model, display, iters);

iters = 100; %Default: 100
display = 0;

fprintf(1, '# Removing some outputs randomly from test data...\n');

% randomly choose which outputs are present.
% Here we assume 50% missing outputs from each test point
numIndPresent = round(0.5*dims);
indicesPresent = zeros(N,numIndPresent);
indicesMissing = zeros(N, dims-numIndPresent);
for i=1:N
    permutat = randperm(dims);
    indicesPresent(i,:) =  permutat(1:numIndPresent);
    indicesMissing(i,:) = setdiff(1:dims, indicesPresent(i,:));
end

% Missing data are indicated with NaN values. This demo creates randomly
% such data.
for i=1:size(Yts,1)
    Yts(i,indicesMissing(i,:))=NaN;
end

indexP = [];
Init = [];
Testmeans = []; 
Testcovars = [];
Varmu = [];
Varsigma = [];

fprintf(1, '# Partial reconstruction of test points...\n');

% patrial reconstruction of test points
for i=1:size(Yts,1)
    indexPresent = indicesPresent(i,:);
    indexP(i,:) = indexPresent;
    
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
    [x, varx] = vargplvmOptimisePoint(model, vardistx, Yts(i, :), display, iters);

    Testmeans(i,:) = x;
    Testcovars(i,:) = varx;
    
    % reconstruct the missing outputs  
    [mu, sigma] = vargplvmPosteriorMeanVar(model, x, varx);
    Varmu(i,:) = mu; 
    Varsigma(i,:) = sigma; 
    
    %
end

% Save the results.
modelWriteResult(model, dataSetName, experimentNo);

if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, lbls, dataSetName, experimentNo);
end

% Load the results and display dynamically.
lvmResultsDynamic(model.type, dataSetName, experimentNo, 'image', [20 28], 1, 0, 1)

