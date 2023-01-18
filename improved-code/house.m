function [v,beta]=house(x)
    % household算法的实现
    % 输入为需要变换的向量x;变换矩阵为H=I-beta*v*v'
    n=length(x);
    x=x/norm(x,inf);
    sigma=x(2:n)'*x(2:n);
    v(2:n)=x(2:n);
    if sigma==0
        beta=0;
    else
        alpha=sqrt(x(1)^2+sigma);
        if x(1)<=0
            v(1)=x(1)-alpha;
        else
            v(1)=-sigma/(x(1)+alpha);
        end
        beta=2*v(1)^2/(sigma+v(1)^2);
        v=v/v(1);
    end
end