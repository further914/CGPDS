% Reconstruct data in cmu16
clear
if ~exist('experimentNo') experimentNo = 1; end
if ~exist('index') index =2; end
randn('seed', 1e5);
rand('seed', 1e5);
[~, lbls] = lvmLoadData('cmu16Jog');
seq = cumsum(sum(lbls)) - [1:5];
dataSetName = 'cmu16gplvm';

% load test data
Y = lvmLoadData(dataSetName);
Ytests{1} = Y(1:seq(1),:);
for i =2:5
    Ytests{i} = Y(seq(i-1)+1:seq(i),:);
end
save('cmu16gpdmTs','Ytests');
Ytest = Ytests{index};% for one test
% save('cmu16JogTs','Ytest');
dt=0.05;
timeStampsTest = ([0:size(Ytest,1)-1].*dt)';

% load data
dataSetName = 'cmu35gplvm';
load([dataSetName 'JogSub']);
origBias = mean(Y);
origScale = 1./sqrt(var(Y));

% load trained model
fileToSave = [pwd filesep 'CGPDSResultsMat' filesep 'cmu16Result' filesep 'demCmu35cgpds' num2str(experimentNo) '.mat'];

load(fileToSave); % From demCmu35gplvmcgpds1.m
fprintf(1,'# Loaded %s\n',fileToSave);


% Indices associated with left leg.
legInd = [8:14];
% Indices associated with upper body.
bodyInd = [21:50];

fprintf(1,'# Dynamics kernel:');
kernDisplay(model.dynamics.kern);
model.dynamics.t_star = timeStampsTest;
kernName = model.dynamics.kern.comp{1}.type;
fprintf(1,'# experimentNo=%d\n',experimentNo);
fprintf(1,'# Reconstruction Iters=%d\n',model.reconstrIters);

model.reconstrIters = 5000; %%%TEMP

startInd = 1;


fileOfmodelUpdated = [pwd filesep 'CGPDSResultsMat' filesep 'cmu16Result' filesep 'demCmu16cgpdsModelUpdated' num2str(experimentNo) '.mat'];


% To compute LegErrs, just replace the parameters 'bodyInd' with 'legInd' and replace 'Body' with 'Leg' in function vargplvmTaylorAngleErrors1
BodyErrs = vargplvmTaylorAngleErrors1(model, Y, Ytest, startInd,  bodyInd, origBias, origScale, [dataSetName ...
    'Body'],  experimentNo, fileOfmodelUpdated);
BodyErrs






