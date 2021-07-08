function [X Xts] = sarcos_all_pca(Ytr,Yts,dims)
   
       [vi,ui] = pca(Ytr);
       X =  Ytr*ui(:, 1:dims)*diag(1./sqrt(vi(1:1:dims)));
       Xts = Yts*ui(:, 1:dims)*diag(1./sqrt(vi(1:1:dims)));
   end
       

