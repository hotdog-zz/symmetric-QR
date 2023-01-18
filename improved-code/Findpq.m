function [p,q]=Findpq(alpha,beta)
    % 寻找当前矩阵中最大的不可约三对角阵
    % 输出p,q后,p+1行、列至n-q行、列为所要找的矩阵
    n=length(alpha);
    q=0;
    n_p_q=0;
    flag=0;
    for i=n-1:-1:1
        if flag==0 && beta(i)~=0
            flag=1;
            n_p_q=1;
        end
        if flag==1 && beta(i)==0
            break;
        end
        if flag==0
            q=q+1;
        end
        if flag==1
            n_p_q=n_p_q+1;
        end
    end
    if q==n-1 && n_p_q==0
        q=n;
    end
    p=n-n_p_q-q;
end
        
      