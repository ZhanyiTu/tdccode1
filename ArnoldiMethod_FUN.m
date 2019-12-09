A=rand(100);
sigma1=max(abs(eig(A)));
e=zeros(10,1);
for n=1:10
    [Q,S]=ArnoldiMethod(A,n);
    sigma2=max(abs(eig(S)));
    e(n)=abs(sigma1-sigma2);
end
figure
plot(abs(eig(A)),'.');
title('The distribution of absolute values of the eigenvalues')
figure
plot(1:10,log(e));
xlabel('n = # of iterations')
ylabel('log(error)')
title('log(errors) with respect to # of iterations')