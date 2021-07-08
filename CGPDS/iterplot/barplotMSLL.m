close all;
clear;
model_series(:,:,1) = [3.85 3.88 4.17 6.20 2.95 ;
2.85 3.36 5.08 2.62 2.59;
2.78 3.42 5.24 2.57 2.52];
model_error(:,:,1) =[0.86 0.19 0.69 0.87 0.14;
 0.32 0.35 0.79 0.04 0.04;
0.27 0.27 0.49 0.03 0.04];

model_series(:,:,2) = [4.06 3.10 3.39 5.73 2.34;
2.44 2.77 3.49 2.06 2.00;
2.05 2.59 3.54 2.02 1.93];
model_error(:,:,2) = [1.18 0.27 0.38 0.89 0.15;
0.37 0.19 0.43 0.04 0.03;
0.32 0.24 0.34 0.03 0.04];
model_series(:,:,3) = [4.33 3.44 2.68 2.36 2.44;
2.39 2.40 2.65 2.12 1.92;
2.28 2.28 2.46 2.04 1.83];
model_error(:,:,3) = [1.54 0.49 0.24 0.08 0.12;
0.59 0.12 0.25 0.15 0.04;
0.51 0.19 0.10 0.08 0.02];

model_series(:,:,4) = [2.76 2.11 1.30 0.86 0.95;
0.92 1.27 1.71 0.59 0.53;
0.76 1.27 2.20 0.55 0.44];
model_error(:,:,4) = [1.53 0.55 0.21 0.66 0.10;
0.51 0.19 0.21 0.10 0.04;
0.45 0.15 0.22 0.04 0.03];


color=[0 0 0.75;0 1 0;1 0.5 0];
for n=1:4
figure;
h = bar(model_series(:,:,n));
ch = get(h,'children');
set(ch{5},'FaceColor','r')
% set(ch{2},'FaceVertexCData',[1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4])
% set(ch{3},'FaceVertexCData',[1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4])
% set(ch{4},'FaceVertexCData',[1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4])
% set(ch{5},'FaceVertexCData',[1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4])
set(gca,'FontSize',13);
set(h,'BarWidth',1);    % The bars will now touch each other
% set(gca,'YGrid','on')
% set(gca,'GridLineStyle','-')
set(gca,'XTicklabel','100|300|500')
set(get(gca,'YLabel'),'String','MSLL','FontSize',15)
set(get(gca,'XLabel'),'String','Training Set Size','FontSize',15)
lh = legend([ch{1} ch{2} ch{3} ch{4} ch{5}],'sGPR','MTGP','COGP','VGPDS','VDMM-GPRM');
% set(lh,'Location','BestOutside','Orientation','horizontal')
hold on;
numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
      errorbar(x, model_series(:,i,n), model_error(:,i,n), 'k', 'linestyle', 'none','LineWidth',2);
end

print('-depsc',[num2str(n) '.eps'])
end
