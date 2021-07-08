function [X Xts] = sarcos_pca(Ytr,Yts,dims)
   jointN = size(Ytr,2)/3;
   X=[];
   Xts=[];
   for i=1:jointN
       chosend = ((i-1)*3+1):i*3;
       [vi,ui] = pca(Ytr(:,chosend));
       Xi =  Ytr(:,chosend)*ui(:, 1:dims)*diag(1./sqrt(vi(1:1:dims)));
       X = [X,Xi];
       if length(Yts)>0
       Xtsi = Yts(:,chosend)*ui(:, 1:dims)*diag(1./sqrt(vi(1:1:dims))); 
       Xts = [Xts,Xtsi];
       end
       
      
   end
       
end

