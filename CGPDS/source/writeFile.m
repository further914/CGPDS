% fid = fopen('demToyMultigp2DVargplvmresult11_20.txt','w');
% fprintf(fid,'RMSEy_1\tRMSEy_2\tRMSEy_3\tRMSEy_4\n');
% for i = 1:10
%     expNo = 10+i;
%     load(['demToyMultigp2DVargplvm' num2str(expNo) 'result2']);
% %     str = [num2str(RMSE(1)) ',' RMSE(2) ',' RMSE(3) ',' RMSE(4) '\n'];
%     fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',RMSE(1),RMSE(2),RMSE(3),RMSE(4));
% end 
% fprintf(fid,'LPy_1\tLPy_2\tLPy_3\tLPy_4\n');
% for i = 1:10
%     expNo = 10+i;
%     load(['demToyMultigp2DVargplvm' num2str(expNo) 'result2']);
%     fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',p(1),p(2),p(3),p(4));
%     
% end 
% fprintf(fid,'vary_1\tvary_2\tvary_3\tvary_4\n');
% for i = 1:10
%     expNo = 10+i;
%     load(['demToyMultigp2DVargplvm' num2str(expNo) 'result2']);
%     fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',var(1),var(2),var(3),var(4));
%     
% end 
% fclose(fid);
% 
% fid = fopen('demToyMultigp2DVarmgplvmresult11_20.txt','w');
% fprintf(fid,'RMSEy_1\tRMSEy_2\tRMSEy_3\tRMSEy_4\n');
% for i = 1:10
%     expNo = 10+i;
%     load(['demToyMultigp2DVarmgplvm' num2str(expNo) 'result2']);
% %     str = [num2str(RMSE(1)) ',' RMSE(2) ',' RMSE(3) ',' RMSE(4) '\n'];
%     fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',RMSE(1),RMSE(2),RMSE(3),RMSE(4));
% end 
% fprintf(fid,'LPy_1\tLPy_2\tLPy_3\tLPy_4\n');
% for i = 1:10
%     expNo = 10+i;
%     load(['demToyMultigp2DVarmgplvm' num2str(expNo) 'result2']);
%     fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',p(1),p(2),p(3),p(4));
%     
% end 
% for i = 1:10
%     expNo = 10+i;
%     load(['demToyMultigp2DVarmgplvm' num2str(expNo) 'result2']);
%     fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',var(1),var(2),var(3),var(4));
%     
% end 
% fclose(fid);
% 
fid = fopen('demToyMultigp2DVargplvmresult41_50.txt','w');
fprintf(fid,'RMSEy_1\tRMSEy_2\tRMSEy_3\tRMSEy_4\n');
for i = 1:10
    expNo = 40+i;
    load(['demToyMultigp2DVargplvm' num2str(expNo) 'result']);
%     str = [num2str(RMSE(1)) ',' RMSE(2) ',' RMSE(3) ',' RMSE(4) '\n'];
    fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',RMSE(1),RMSE(2),RMSE(3),RMSE(4));
end 
fprintf(fid,'LPy_1\tLPy_2\tLPy_3\tLPy_4\n');
for i = 1:10
    expNo = 40+i;
    load(['demToyMultigp2DVargplvm' num2str(expNo) 'result']);
    fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',p(1),p(2),p(3),p(4));
    
end 
fclose(fid);

fid = fopen('demToyMultigp2DVarmgplvmresult41_50.txt','w');
fprintf(fid,'RMSEy_1\tRMSEy_2\tRMSEy_3\tRMSEy_4\n');
for i = 1:10
    expNo = 40+i;
    load(['demToyMultigp2DVarmgplvm' num2str(expNo) 'result']);
%     str = [num2str(RMSE(1)) ',' RMSE(2) ',' RMSE(3) ',' RMSE(4) '\n'];
    fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',RMSE(1),RMSE(2),RMSE(3),RMSE(4));
end 
fprintf(fid,'LPy_1\tLPy_2\tLPy_3\tLPy_4\n');
for i = 1:10
    expNo = 40+i;
    load(['demToyMultigp2DVarmgplvm' num2str(expNo) 'result']);
    fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',p(1),p(2),p(3),p(4));
    
end 
fclose(fid);