clc;
clear all;

% tol为误差判定,n为矩阵阶数,x为特征值,Q为正交矩阵,A为对称矩阵
% eigenvalue为运行得到的原矩阵的特征值
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
% 三幅图为误差一、二、三的计算
figure(1)
plot(sort(x),((sort(eigenvalue)-sort(x))/norm(A))./sort(x));
xlabel("原特征值")
ylabel("特征值之差")
figure(2)
imagesc((U*diag((eigenvalue))*U'-A)/norm(A));
colorbar
figure(3)
imagesc(U*U'-eye(n,n));
colorbar
%}