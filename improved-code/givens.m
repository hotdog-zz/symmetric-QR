function [c,s]=givens(a,b)
    % givens�任��ʵ��
    % a,bΪ������������Ҫ��b��Ϊ0;[c,s;-s,c]��Ϊgivens�任����
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