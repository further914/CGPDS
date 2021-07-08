close all;
clear;
model_series(:,:,1) = [6.33	5.52	5.11	4.89	4.55;
4.19	4.37	3.78	3.67	3.52;
4.01	3.5	3.33	3.47	3.21];
model_error(:,:,1) =[1.33	0.79	0.44	0.36	0.34;
0.83	0.89	0.14	0.18	0.09;
0.71	0.32	0.08	0.17	0.11];

model_series(:,:,2) = [4.48	3.2	3.18	2.96	2.68;
2.61	2.53	2.27	2.11	2.01;
2.08	1.97	1.89	1.98	1.77];
model_error(:,:,2) = [1.04	0.38	0.23	0.26	0.22;
0.51	0.56	0.09	0.12	0.06;
0.45	0.27	0.06	0.12	0.08];
model_series(:,:,3) = [5.01	3.27	3.26	3.19	3.19;
2.21	2.01	2.25	2.31	1.95;
2.08	1.66	1.77	1.93	1.64];
model_error(:,:,3) = [2.06	0.35	0.25	0.2	0.26;
0.37	0.2	0.18	0.38	0.09;
0.38	0.17	0.06	0.14	0.06];

model_series(:,:,4) = [1.04	0.72	0.68	0.65	0.66;
0.55	0.47	0.52	0.48	0.45;
0.47	0.39	0.45	0.44	0.39];
model_error(:,:,4) = [0.22	0.07	0.05	0.03	0.04;
0.13	0.04	0.02	0.06	0.01;
0.08	0.02	0.02	0.01	0.01];


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
set(get(gca,'YLabel'),'String','RMSE','FontSize',15)
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
