axisN=10*[10 30 50];

svl=[1.04 0.55 0.47];
cov1=[0.22 0.13 0.08];
sv2=[0.72 0.47 0.39];
cov2=[0.07 0.04 0.02];
sv3=[0.68 0.52 0.45];
cov6=[0.05 0.02 0.02];
svm_2k=[0.65 0.48 0.44];
cov4=[0.03 0.06 0.01];
mvmed=[0.66 0.45 0.39];
cov3=[0.04 0.01 0.01];


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
print -depsc 7.eps
%%

