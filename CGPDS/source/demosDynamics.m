% train the CGPDS on missa and dog dataset

close all


%%
% ############################ MISSA ###################################
clear; load opts
fprintf(1,'\n\n#-----  MISSA DEMO ----#\n');

experimentNo = 30;
dataSetName = 'missa';
indPoints = -1; latentDim=12;
initialBetaIters = 0;
reconstrIters = 10000;
itNo=[2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000]; %
dynamicKern={'matern32','white','bias'};
vardistCovarsMult=1.6;
dataSetSplit = 'blocks';
blockSize = 2;
whiteVar = 0.1;
msk = [48 63 78 86 96 111 118];
demHighDimCgpds





%%
% ############################ DOG #######################################
% Train
clear; load opts
fprintf(1,'\n\n#-----  DOG DEMO: Reconstruction ----#\n');
dataSetName = 'dog';
experimentNo=68;
indPoints=-1; latentDim=6;
reconstrIters = 8000;
itNo=[1000 1000 1000 1000 1000 1000 500 500 500 500 1000 1000 1000 1000 1000 1000 1000 500 500]; %16000
periodicPeriod = 2.8840;
dynamicKern={'rbfperiodic','whitefixed','bias','rbf'};
vardistCovarsMult=0.8;
whiteVar = 1e-6;
dataSetSplit = 'custom';
indTr = 1:54;
indTs = 55:61;
learnSecondVariance = 0;
demHighDimCgpds

%%

% % Test
% clear; load opts
% dataSetName = 'dog';
% experimentNo=65;
% dataSetSplit = 'custom';
% indTr = 1:54; indTs = 55:61;
% predWithMs = 1; % Do reconstruction
% reconstrIters = 18000; 
% doSampling = 0;
% if rePredict
%     demHighDimVargplvmTrained
% else
%     load demDogVargplvm65Pred
%     demHighDimVargplvmLoadPred
% end
% % Produce plots (these go to the supplementary)
% fr=reshape(Varmu(5,:),height,width); imagesc(fr); colormap('gray'); 
% print -depsc ../diagrams/supplDogPredGpds5.eps; system('epstopdf ../diagrams/supplDogPredGpds5.eps');
% fr=reshape(Yts(5,:),height,width); imagesc(fr); colormap('gray'); 
% print -depsc ../diagrams/supplDogPredYts5.eps; system('epstopdf ../diagrams/supplDogPredYts5.eps');
% fr=reshape(Varmu(6,:),height,width); imagesc(fr); colormap('gray'); 
% print -depsc ../diagrams/supplDogPredGpds6.eps; system('epstopdf ../diagrams/supplDogPredGpds6.eps');
% fr=reshape(Yts(6,:),height,width); imagesc(fr); colormap('gray'); 
% print -depsc ../diagrams/supplDogPredYts6.eps; system('epstopdf ../diagrams/supplDogPredYts6.eps');

%%

% ############################ CMU ###################################
%rbf
% clear; load opts
% fprintf(1,'\n\n#-----  CMU DEMO: Rbf ----#\n');
% experimentNo=34; 
% itNo = [300 300 400 200 200 300 400 400];
% dynamicKern = {'rbf', 'white', 'bias'};
% vardistCovarsMult = 0.152;
% if retrainModels 
%     if rePredict
%         doReconstr = 1;
%     else
%         doReconstr=0;
%     end
%     demCmu35gplvmVargplvm3;
% elseif rePredict
%     % 'demCmu35gplvmVargplvm34.mat' must be included in the path
%     demCmu35vargplvmReconstructTaylor
% end
% predictPart = 'Legs';  plotRange = [];
% demCmu35VargplvmPlotsScaled
% % demCmu35VargplvmAnimate
% fprintf(1,'# VGPDS RBF error on Legs reconstr:');
% errStruct
% 
% predictPart = 'Body';  plotRange = [];
% demCmu35VargplvmPlotsScaled
% % demCmu35VargplvmAnimate
% fprintf(1,'# VGPDS RBF error on Body reconstr:');
% errStruct
% 
% bar(model.kern.comp{1}.inputScales);
% print -depsc ../diagrams/supplMocapScalesRbf.eps; system('epstopdf ../diagrams/supplMocapScalesRbf.eps');
% 
% 
% % matern32 for legs
% clear; load opts
% fprintf(1,'\n\n#-----  CMU DEMO: Matern32 ----#\n');
% experimentNo=33; 
% itNo = [300 300 400 200 200 300 400 400];
% dynamicKern = {'matern32', 'white', 'bias'};
% vardistCovarsMult = 0.24;
% if retrainModels 
%     if rePredict
%         doReconstr = 1;
%     else
%         doReconstr=0;
%     end
%     demCmu35gplvmVargplvm3;
% elseif rePredict
%     % 'demCmu35gplvmVargplvm33.mat' must be included in the path
%     demCmu35vargplvmReconstructTaylor
% end
% predictPart = 'Legs'; plotRange = 10;
% demCmu35VargplvmPlotsScaled
% print -depsc ../diagrams/supplMocapLeg5GpdsMatern.eps; system('epstopdf ../diagrams/supplMocapLeg5GpdsMatern.eps');
% % demCmu35VargplvmAnimate
% fprintf(1,'# VGPDS Matern error on Legs reconstr:');
% errStruct
% 
% predictPart = 'Body'; plotRange = 28;
% demCmu35VargplvmPlotsScaled
% print -depsc ../diagrams/supplMocapBody28GpdsMatern.eps; system('epstopdf ../diagrams/supplMocapBody28GpdsMatern.eps');
% % demCmu35VargplvmAnimate
% fprintf(1,'# VGPDS Matern error on Body reconstr:');
% errStruct
% close all
% bar(model.kern.comp{1}.inputScales);
% print -depsc ../diagrams/supplMocapScalesMatern.eps; system('epstopdf ../diagrams/supplMocapScalesMatern.eps');


%% ---------
% fprintf(1,'\n\n#---- FINISHED reproducing plots and results!! \n');
% delete opts.mat
% 
% 
% % Record version of MATLAB/Octave
% a = ver('octave');
% if length(a) == 0
%   a = ver('matlab');
% end
% fid = fopen('vers.tex', 'w');
% fprintf(fid, [a.Name ' version ' a.Version]);
% fclose(fid);
% 
% % Record computer architecture.
% fid = fopen('computer.tex', 'w');
% fprintf(fid, ['\\verb+' computer '+']);
% fclose(fid);
% 
% % Record date of run.
% fid = fopen('date.tex', 'w');
% fprintf(fid, datestr(now, 'dd/mm/yyyy'));
% fclose(fid);

