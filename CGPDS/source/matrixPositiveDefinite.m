function [ A ] = matrixPositiveDefinite( A )
% ensure a matrix to be positive definite
if det(A) < 1e-20
    A=(A+A')/2;%�Գ���
    n=size(A,1);%A
    [v,l]= eig(A);%���ؾ������������v��ÿһ�ж�Ӧ��ĳ������ֵ��������������l������ֵ����һ���Խ��󣬶Խ���Ԫ��Ϊ����ֵ��
    l(find(l<0))=1e-100;%��С��0������ֵ��ֵΪ0
    A=v*l*v';%�õ�����������
    eps = 10^-6;
    A=A+eye(n)*eps;% ��������
end
end

