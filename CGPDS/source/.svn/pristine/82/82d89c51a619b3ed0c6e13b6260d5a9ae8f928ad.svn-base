onlyJog=1;
[Y, lbls] = lvmLoadData('cmu35WalkJog');
seq = cumsum(sum(lbls)) - [1:31];
dataSetName = 'cmu35gplvm';
% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
if exist('onlyJog')&& onlyJog
    Y=Y(seq(16)+1:seq(25),:);
    seq=seq-seq(16);
    seq=seq(17:17+8);
end
save([dataSetName 'Jog'],'Y');
fprintf(1,'\n#----------------------------------------------------\n');
fprintf(1,'# Dataset: %s\n',dataSetName);
fprintf(1,'#----------------------------------------------------\n');
% Fix times:
prevSeq = 1;
timeStampsTraining = [];
dt=0.05;
Y=Y(1:seq(1),:);
seq=seq(1:1);
for i=1:length(seq)
    t = ([0:(seq(i)-prevSeq)].*dt)';
    prevSeq = seq(i)+1;
    timeStampsTraining = [timeStampsTraining ;t];
end;
timeStampsTest = ([0:size(Ytest,1)-1].*dt)';
%% 先整个转换存储方式
for d=1:size(Y,2)
    XTemp{d}=[timeStampsTraining;timeStampsTest];
    yTemp{d}=[Y(:,d);Ytest(:,d)];
end
%% 去除掉miss data
% Indices associated with right leg.
legInd = [8:14];

% Indices associated with upper body.
bodyInd = [21:50];

for d=legInd
    XTemp{d}=timeStampsTraining;
    yTemp{d}=Y(:,d);
end

XTestTemp=timeStampsTest;
for d=1:size(Y,2)
    yTestTemp{d}=Ytest(:,d);
end
save('cmu35cmogpL','XTemp', 'yTemp', 'XTestTemp','yTestTemp');