axisN=10*[10 30 50];

svl=[5.01 2.21 2.08];
cov1=[1.06 0.37 0.38];
sv2=[3.27 2.01 1.66];
cov2=[0.35 0.20 0.17];
sv3=[3.26 2.25 1.77];
cov6=[0.25 0.18 0.06];
svm_2k=[3.19 2.31 1.93];
cov4=[0.20 0.38 0.14];
mvmed=[3.19 1.95 1.64];
cov3=[0.26 0.09 0.06];


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
print -depsc 4.eps
%%

