function [Q,A]=tridiagonal(A)
    % 矩阵三对角化的实现
    % Q为变换矩阵，A的主对角线和次对角线为三对角化后的矩阵元素
    n=length(A);
    Q=eye(n,n);
    for k=1:n-2
        [v,beta]=house(A(k+1:n,k));
        v=v';
        u=beta*A(k+1:n,k+1:n)*v;
        w=u-(beta*u'*v/2)*v;
        A(k+1,k)=norm(A(k+1:n,k),2);
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-v*w'-w*v';
        
        w=beta*Q(1:n,k+1:n)*v;
        Q(1:n,k+1:n)=Q(1:n,k+1:n)-w*v';
    end
end