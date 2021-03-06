
% DEMUSPSVARGPLVM1 Demonstrate variational GPLVM on USPS data.
%
%	Description:
%	% 	demUspsVargplvm1.m SVN version 583
% 	last update 2009-11-08T13:07:33.000000Z

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'usps';
experimentNo = 1;
printDiagram = 1;

% load data
[YTrain, lblsTrain, YTest, lblsTest] = lvmLoadData(dataSetName);

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'white'};
options.numActive = 50; 
%options.tieParam = 'tied';

options.optimiser = 'scg';
latentDim = 10;
d = size(YTrain, 2);

iters = 1000;
display = 1;

% create a separate vargplvm for each digit
for i=1:10
    %
    Y = YTrain(lblsTrain(:,i)==1,:);
 
    model = vargplvmCreate(latentDim, d, Y, options);
    %
    model = vargplvmParamInit(model, model.m, model.X); 

    model = vargplvmOptimise(model, display, iters);
    
    varmodel{i} = model;
    % 
end  

iters = 100;
display = 0;

capName = dataSetName;
capName(1) = upper(capName(1));
modelType = varmodel{1}.type;
modelType(1) = upper(modelType(1));
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'varmodel');

% measure performance on test data 
indexPresent = 1:size(YTest,2);
TestError = 0;
for n=1:size(YTest,1)
    %
    % compute the approximate class conditional density for each digit
    for i=1:10
       prob(n,i) = vargplvmProbabilityCompute(varmodel{i}, YTest(n,:), indexPresent);
    end
    [maxP C] = max(prob(n,:));
    if lblsTest(n,C) == 0
        TestError = TestError + 1; 
    end
    %
end

%
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'varmodel', 'prob', 'TestError');
