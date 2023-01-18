function [alpha,Q]=symmetricQR(A,tol)
    % ��ʽ�Գ�QR�㷨��ʵ��
    % tolΪ���ζԽ�Ԫ��Ϊ0ʱ��u;���alphaΪ����ֵ��QΪ�任����
    
    % ���Խǻ�
    n=max(size(A));
    [Q,A]=tridiagonal(A);
    % �洢���Խǻ���Խ�Ԫalpha�ʹζԽ�Ԫbeta
    alpha=diag(A);
    beta=diag(A,-1);
    clear A;
    % QR������ֱ���ζԽ�Ԫ��Ϊ0ֹͣ
    q=0;
    while q<n
        % ���������ĴζԽ�ԪԪ������
        for i=1:n-1
            if abs(beta(i))<=tol*(abs(alpha(i))+abs(alpha(i+1)))
                beta(i)=0;
            end  
        end
        % Ѱ�����Ĳ���Լ���ԽǾ���
        [p,q]=Findpq(alpha,beta);
        % QR����
        if q<n
            [alpha(p+1:n-q),beta(p+1:n-q-1),G]=qrshifts(alpha(p+1:n-q),beta(p+1:n-q-1));
            Q(1:n,p+1:n-q)=Q(1:n,p+1:n-q)*G;
        end
    end
end