function [c,s]=givens(a,b)
    % givens变换的实现
    % a,b为输入向量，需要将b变为0;[c,s;-s,c]即为givens变换矩阵
    if b==0
        c=1;
        s=0;
    else
        if abs(b)>abs(a)
            tau=a/b;
            s=1/sqrt(1+tau^2);
            c=s*tau;
        else
            tau=b/a;
            c=1/sqrt(1+tau^2);
            s=c*tau;
        end
    end
end