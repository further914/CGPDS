function [errStruct] = vargplvmTaylorAngleErrors1(model, Y, Ytest, startInd,...
    missingInd, origBias, origScale, name, expNo, fileOfmodelUpdated)
% VARGPLVMTAYLORANGLEERRORS Helper function for computing angle errors for CMU 35 data using
%	GPLVM with dynamics.
%
%	Description:
%	[errStruct] = vargplvmTaylorAngleErrors(model, Y, Ytest, startInd, origBias, origScale,...
%                                              missingInd, name, expNo);
% 	Based on fgplvmTaylorAngleErrors.m for FGPLVM.
% COPYRIGHT: Neil Lawrence, Andreas Damianou, Michalis Titsias 2010-2011

%%%NEW

if ~exist('experimentNo') experimentNo = expNo; end

Ygplvm=Y;%353*9
YtestGplvm=Ytest;%cmu35 43*59 unnormlized cmu16 47*59
YtestGplvmOrig = Ytest;

Y = Y - repmat(origBias, size(Y, 1), 1);
Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
Y = Y.*repmat(origScale, size(Y, 1), 1);
Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);
%%%

origYtest = Ytest;
YtrueTest = Ytest;%after normlising

capName = name;
capName(1) = upper(capName(1));

temp = ones(1, size(Ytest, 2));
temp(missingInd) = 0;
presentInd = find(temp);

YtrainNn = Y(:, presentInd);%353*52
YtestNn = origYtest(startInd:end, presentInd);%43*52
YtsTmp = YtestGplvmOrig(startInd:end, :);

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
errStruct.taylorRMSENnScaled = sqrt(mean(mean((err.*err))));

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
errStruct.taylorRMSENn = sqrt(mean(mean((err.*err))));

Ytest(startInd:end, missingInd) = NaN;

%%%%%%%%%%%%%  VAR  GPLVM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YtestGplvm(startInd:end, missingInd) = NaN;

disp('# Initializing latent points...');
indexPresent = setdiff(1:model.d, missingInd);

%iters=2000;% default: 1000
iters = model.reconstrIters;
fprintf(1,'# Optimising for %d iterations...\n',iters);

for i=1:size(YtestGplvm,1)
    % initialize the latent points using the nearest neighbour from the training data
    dst = dist2(YtestGplvm(i,indexPresent), Ygplvm(:,indexPresent));
    [mind, mini(i)] = min(dst);
end
indexMissingData = startInd:size(YtestGplvm);
vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
vardistx.covars = 0.2*ones(size(vardistx.covars));
model.vardistx = vardistx;


try
    load(fileOfmodelUpdated);
catch
    samd = [1:model.d];
    samr = [1:1];
    % Do also reconstruction in test data
    [x, varx] = cgpdsOptimiseSeqDyn(model, vardistx, YtestGplvm, 1, iters, samd, samr);
    
    % keep the optimized variational parameters
    barmu = x;
    lambda = varx;
    % Get the variational means and variacnes for the new test sequcen and
    % update the model to be prepared for prediction
    [x, varx, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(model, barmu, lambda, YtestGplvm);
    modelUpdated.origBias = origBias;
    modelUpdated.origScale = origScale;
    save(fileOfmodelUpdated,'x','varx','modelUpdated');
end

% Latent variables corresponding to the data  with missing dimensions
Testmeans = x(indexMissingData, :);
Testcovars = varx(indexMissingData, :);


%  YtsOrMs = YtsOriginal(indexMissingData,:);
disp('# Making the predictions for Y...')
% Reconstruct all dimensions (missing and observed) for the test data
[mu, sigma] = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
Ypred = mu;
YpredVar = sigma;

disp('# Final calculations...')

YtrueTest=YtestGplvmOrig;%unnormalised
% YtrueTest=origYtest;%unnormalised
% % De-normalise the predictions.
Ypred = Ypred./repmat(origScale, size(Ypred, 1), 1);
Ypred = Ypred + repmat(origBias, size(Ypred, 1), 1);


% err = (YtrueTest(startInd:end, missingInd) - Ypred(:, missingInd))./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);
err = YtrueTest(startInd:end, missingInd) - Ypred(:, missingInd);
errStruct.ErrorCGpds = sqrt(mean(mean((err.*err))));  

% Convert to degrees.
err = err*180/pi;
% Compute average of mean square error.
errStruct.angleErrorCGpds = sqrt(mean(mean((err.*err))));
load cmu35TaylorScaleBias
err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
errStruct.taylorErrorCGpds = sum(sum(err.*err))/length(missingInd);
errStruct.taylorRMSECGpds = sqrt(mean(mean((err.*err))));


likehood = 0;
pdf4y = zeros(size(YtestGplvm,1),1);

for n = 1:size(YtestGplvm,1)
    likehood = likehood - length(missingInd)/2*log(2*pi)-0.5*log(det(diag(YpredVar(n,missingInd))))-0.5*(YtrueTest(n,missingInd)-Ypred(n,missingInd))*...
        diag(YpredVar(n,missingInd))^-1*(YtrueTest(n,missingInd)-Ypred(n,missingInd))';
    
    pdf4y(n)=mvnpdf(YtrueTest(n,missingInd),Ypred(n,missingInd),diag(YpredVar(n,missingInd)));
end
errStruct.p = - likehood/size(YtestGplvm,1);
errStruct.p1 = - mean(pdf4y,1);

% errStruct.p=calculatePosteriorLikelihoodCMGPDS(model,Testmeans,Testcovars,YtsTmp,missingInd,ones(1,length(missingInd)),origBias,origScale);
% errStruct.pa=calculatePosteriorLikelihoodCMGPDS(model,Testmeans,Testcovars,YtsTmp,missingInd,180/pi*ones(1,length(missingInd)),origBias,origScale);
% errStruct.ps=calculatePosteriorLikelihoodCMGPDS(model,Testmeans,Testcovars,YtsTmp,missingInd,180/pi*scale(1, missingInd+9),origBias,origScale);
% if isinf(errStruct.pa) || isinf(errStruct.ps)||isinf(errStruct.p)
% errStruct.p=calculatePosteriorLikelihoodCMGPDS2(model,Testmeans,Testcovars,YtsTmp,missingInd,ones(1,length(missingInd)));
% errStruct.pa2=calculatePosteriorLikelihoodCMGPDS2(model,Testmeans,Testcovars,YtsTmp,missingInd,180/pi*ones(1,length(missingInd)));
% errStruct.ps2=calculatePosteriorLikelihoodCMGPDS2(model,Testmeans,Testcovars,YtsTmp,missingInd,180/pi*scale(1, missingInd+9));
% end


% plotRange = missingInd;
% colordef white   
% for plotNo = plotRange
%     figNo = plotNo - min(plotRange) + 1;
%     figure(figNo);
%     clf
% %     lin = plot(1:size(origYtest, 1), origYtest(:, plotNo), '-');
%       lin = plot(1:size(YtrueTest, 1), YtrueTest(:, plotNo), '-');
%     hold on
% %     lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Y(bestIndNn, plotNo)], ':')];
% %     lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Ypred(:, plotNo)], '--')];
% 
% %     lin = [lin; plot(1:size(YtestGplvmOrig, 1), [YtestGplvmOrig(1:startInd-1, plotNo); Y(bestIndNn, plotNo)], ':')];
%     lin = [lin; plot(1:size(YtrueTest, 1), [YtrueTest(1:startInd-1, plotNo); Ypred(:, plotNo)], '--')];
%     legend('Ytrue','Ypre');
%     ax = gca;
%     xlabel('Time');
%     ylabel('Output');
%     set(ax, 'fontname', 'arial');
%     set(ax, 'fontsize', 20);
%     set(lin, 'lineWidth', 2);
% %     fileName = ['dem' capName 'Reconstruct' num2str(expNo) '_' num2str(plotNo)];
%     %  print('-depsc', ['../tex/diagrams/' fileName])
%     %  print('-deps', ['../tex/diagrams/' fileName 'NoColour'])
%     %plot2svg(['../html/' fileName '.svg'])
%     
%     % make smaller for PNG plot.
%     pos = get(gcf, 'paperposition');
%     origpos = pos;
%     pos(3) = pos(3)/2;
%     pos(4) = pos(4)/2;
%     set(gcf, 'paperposition', pos);
%     fontsize = get(gca, 'fontsize');
%     set(gca, 'fontsize', fontsize/2);
%     lineWidth = get(gca, 'lineWidth');
%     set(gca, 'lineWidth', lineWidth*2);
%     %print('-dpng', ['../html/' fileName])
%     %plot2svg(['../html/' fileName '.svg'])
%     set(gcf, 'paperposition', origpos);
% end

%%%
%prunedModelUpdated = vargplvmPruneModel(modelUpdated,1);
% fileToSave = 'demCmu35Vargplvm';
% if strcmp(name(end-3:end),'Body')
%     n = 'PredBody.mat';
% else
%     n='PredLegs.mat';
% end
% 
% if isfield(model, 'dynamics') & ~isempty(model.dynamics)
% 	prunedModel = vargplvmPruneModel(model);
% 	prunedModelUpdated = vargplvmPruneModel(modelUpdated);
% 	save([fileToSave num2str(expNo) n], 'prunedModel','prunedModelUpdated','Ypred','errStruct');
% else
% 	save([fileToSave num2str(expNo) n], 'model','Ypred','errStruct');
% end
% fprintf(1,'# Saved %s\n',[filename ]);
%%%
