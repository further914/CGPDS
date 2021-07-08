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
% Y=[Y;Ytest];
seq(end) = [];
segments=[1,seq+1];
save('cmu35gpdmTr','Y', 'segments');
segments = 1;
save('cmu35gpdmTs','Ytest','segments');