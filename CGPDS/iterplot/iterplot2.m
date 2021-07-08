axisN=10*[10 30 50];
svl=[ 6.33 4.19 4.01];
cov1=[ 1.33 0.83 0.71];
sv2=[5.52 4.37 3.50];
cov2=[0.79 0.89 0.32];
sv3=[5.11 3.78 3.33];
cov6=[0.44 0.14 0.08];
svm_2k=[4.89  3.67 3.47];
cov4=[0.36 0.18 0.17];
mvmed=[4.55 3.52 3.21];
cov3=[0.34 0.09 0.11];

clf
errorbar(axisN,svl,cov1,'m-','LineWidth',1.5);
hold on
errorbar(axisN,sv2,cov2,'b-v','LineWidth',1.5);
hold on
errorbar(axisN,sv3,cov6, 'g-<','LineWidth',1.5);
hold on
errorbar(axisN,svm_2k,cov4,'k-*','LineWidth',1.5) 
hold on
errorbar(axisN,mvmed,cov3,'r-o','LineWidth',1.5);

% hold on
% errorbar(axisN,medmk,cov5,'b-s');
 

h=gca;
set(h,'FontSize',16);
xlabel('Training Set Size');
ylabel('RMSE');
legend('sGPR','MTGP','COGP','VGPDS','VDM-GPDS');
print -depsc 2.eps
%%
