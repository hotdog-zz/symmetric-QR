function [alpha,Q]=symmetricQR(A,tol)
    % 隐式对称QR算法的实现
    % tol为将次对角元变为0时的u;输出alpha为特征值，Q为变换矩阵
    
    % 三对角化
    n=max(size(A));
    [Q,A]=tridiagonal(A);
    % 存储三对角化后对角元alpha和次对角元beta
    alpha=diag(A);
    beta=diag(A,-1);
    clear A;
    % QR迭代，直到次对角元均为0停止
    q=0;
    while q<n
        % 满足条件的次对角元元素置零
        for i=1:n-1
            if abs(beta(i))<=tol*(abs(alpha(i))+abs(alpha(i+1)))
                beta(i)=0;
            end  
        end
        % 寻找最大的不可约三对角矩阵
        [p,q]=Findpq(alpha,beta);
        % QR迭代
        if q<n
            [alpha(p+1:n-q),beta(p+1:n-q-1),G]=qrshifts(alpha(p+1:n-q),beta(p+1:n-q-1));
            Q(1:n,p+1:n-q)=Q(1:n,p+1:n-q)*G;
        end
    end
end