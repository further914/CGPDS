function [ A ] = matrixPositiveDefinite( A )
% ensure a matrix to be positive definite
if det(A) < 1e-20
    A=(A+A')/2;%对称阵
    n=size(A,1);%A
    [v,l]= eig(A);%返回矩阵的特征向量v（每一列对应于某个特征值的特征向量），l是特征值，是一个对角阵，对角线元素为特征值。
    l(find(l<0))=1e-100;%将小于0的特征值赋值为0
    A=v*l*v';%得到半正定矩阵
    eps = 10^-6;
    A=A+eye(n)*eps;% 正定矩阵
end
end

