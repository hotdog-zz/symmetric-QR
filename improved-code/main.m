clc;
clear all;

% tolΪ����ж�,nΪ�������,xΪ����ֵ,QΪ��������,AΪ�Գƾ���
% eigenvalueΪ���еõ���ԭ���������ֵ
tol=1e-14;
n=2000;
x=unifrnd(-100,100,n,1);
x=sort(x);
D=diag(x);
Q=rand(n);
[Q,R]=qr(Q);
A=Q*D*Q';
[eigenvalue,U]=symmetricQR(A,tol);
plot(sort(eigenvalue)-sort(x))


%{
% ����ͼΪ���һ���������ļ���
figure(1)
plot(sort(x),((sort(eigenvalue)-sort(x))/norm(A))./sort(x));
xlabel("ԭ����ֵ")
ylabel("����ֵ֮��")
figure(2)
imagesc((U*diag((eigenvalue))*U'-A)/norm(A));
colorbar
figure(3)
imagesc(U*U'-eye(n,n));
colorbar
%}