function [Q,S]=ArnoldiMethod(A,n)
    m=length(A);
    b=rand(m,1);
    S=zeros(n+1,n);
    Q(:,1)=b/norm(b);
    for i=2:n+1
        Q(:,i)=A*Q(:,i-1);
        for j=1:i-1
            S(j,i-1)=Q(:,i)'*Q(:,j);
        end
        Q(:,i)=Q(:,i)-Q(:,1:i-1)*S(1:i-1,i-1);
        S(i,i-1)=norm(Q(:,i));
        Q(:,i)=Q(:,i)/S(i,i-1);
    end
    S=S(1:n,:);
    Q=Q(:,1:n);
end