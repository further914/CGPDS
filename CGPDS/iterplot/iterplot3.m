axisN=10*[10 30 50];

svl=[ 4.48 2.61 2.08];
cov1=[1.04 0.51 0.45 ];
sv2=[3.20 2.53 1.97];
cov2=[0.38 0.56 0.27];
sv3=[3.18 2.27 1.89];
cov6=[0.23 0.09 0.06];
svm_2k=[2.96 2.11 1.98];
cov4=[0.26 0.12 0.12];
mvmed=[2.68 2.01 1.77];
cov3=[0.22 0.06 0.08];


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
print -depsc 3.eps
%%

