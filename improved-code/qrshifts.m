function [alpha,beta,Q]=qrshifts(alpha,beta)
    % QR迭代的实现
    % 输出Q为变换矩阵，alpha,beta为迭代一次后的对角元和次对角元
    n=length(alpha);
    Q=eye(n,n);
    d=(alpha(n-1)-alpha(n))/2;
    miu=alpha(n)-beta(n-1)^2/(d+sign(d)*sqrt(d^2+beta(n-1)^2));
    x=alpha(1)-miu;
    z=beta(1);
    [c,s]=givens(x,z);
    G=[c,s;-s,c];
    Q(1:n,1:2)=Q(1:n,1:2)*G';
    if n==2
        T=[alpha(1),beta(1);beta(1),alpha(2)];
        T=G*T*G';
        alpha(1)=T(1,1);
        beta(1)=T(2,1);
        alpha(2)=T(2,2);
    else
        T=[alpha(1),beta(1);beta(1),alpha(2)];
        T=G*T*G';
        % gamma为每次变换后次次对角元元素
        gamma=s*beta(2);
        beta(2)=c*beta(2);
    end

    for k=2:n-1
        x=T(2,1);
        z=gamma;
        [c,s]=givens(x,z);
        G=[c,s;-s,c];
        Q(1:n,k:k+1)=Q(1:n,k:k+1)*G';
        if k<n-1
            alpha(k-1)=T(1,1);
            beta(k-1)=c*T(2,1)+s*gamma;
            gamma=s*beta(k+1);
            beta(k+1)=c*beta(k+1);
            T=[T(2,2),beta(k);beta(k),alpha(k+1)];
            T=G*T*G';
        else
            alpha(k-1)=T(1,1);
            beta(k-1)=c*T(2,1)+s*gamma;
            T=[T(2,2),beta(k);beta(k),alpha(k+1)];
            T=G*T*G';
            alpha(k)=T(1,1);
            alpha(k+1)=T(2,2);
            beta(k)=T(2,1);
        end
    end
end