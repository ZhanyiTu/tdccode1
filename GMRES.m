function x=GMRES(A,b,n)
    S=zeros(n+1,n);
    t=zeros(n,1);
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
    Q=Q(:,1:n);
    bb=norm(b)
    b=[norm(b) zeros(1,n)]';
    t=inv(S'*S)*S'*b;
    x=Q*t;
end