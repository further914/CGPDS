
% DEMCMU35CGPDSRECONSTRUCTTAYLOR Reconstruct right leg and body of CMU 35.
%
%	Description:
%	
%	DESC
%	Load a model learned with cgpds and reconstruct the missing indices
%	using the Taylor framework developed for cgpds (angle errors)

randn('seed', 1e5);
rand('seed', 1e5);
onlyJog = 1;
onlyWalk = 0;
dynUsed=1;
dataSetName = 'cmu35gplvm';
capName = dataSetName;
capName(1) = upper(capName(1));

fileToSave = [pwd filesep 'CGPDSResultsMat' filesep 'cmu35Result' filesep 'demCmu35cgpds' num2str(experimentNo) '.mat'];

if exist('onlyJog')&& onlyJog
    load([dataSetName 'JogSub']);
elseif exist('onlyWalk') && onlyWalk
    load([dataSetName 'Walk']);
end

origBias = mean(Y);
origScale = 1./sqrt(var(Y));

% REMOVE LEG TEST
% Indices associated with right leg.
legInd = [8:14];
% Indices associated with upper body.
bodyInd = [21:50];

dt=0.05;
Ytest = Ytests{1};
timeStampsTest = ([0:size(Ytest,1)-1].*dt)';

% Load saved model.
load(fileToSave); % From demCmu35gplvmVargplvm1.m
fprintf(1,'# Loaded %s\n',fileToSave);

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    fprintf(1,'# Dynamics kernel:');
    kernDisplay(model.dynamics.kern);
    model.dynamics.t_star = timeStampsTest;
    kernName = model.dynamics.kern.comp{1}.type;
else
    kernName = 'static';
end
fprintf(1,'# experimentNo=%d\n',experimentNo);
fprintf(1,'# Reconstruction Iters=%d\n',model.reconstrIters);

if exist('onlyJog')&& onlyJog
    startInd=1;
end
if exist('onlyWalk')&& onlyWalk
    startInd=1;
end

fileOfmodelUpdated = [pwd filesep 'CGPDSResultsMat' filesep 'cmu35Result' filesep 'demCmu35cgpdsModelUpdated' num2str(experimentNo) '.mat'];

%for leg error, just replace paramters 'bodyInd' with 'legInd' and 'Body'
%with 'Leg'
BodyErrs = vargplvmTaylorAngleErrors1(model, Y, Ytest, startInd, bodyInd, origBias, origScale, [dataSetName ...
    'Body'], experimentNo, fileOfmodelUpdated);
BodyErrs





