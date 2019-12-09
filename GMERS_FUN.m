format long
v=zeros(100,1);
e=zeros(4,1);
N=[4 5 10];
e_=zeros(20,5);
for i=1:100
    r=randperm(3);
    v(i)=N(r(1));
end
v_=v;
v=rand(100,1)+v;
A=diag(v);
P=rand(100);
A=P*A*P^-1;
b=rand(100,1);
x0=A\b;
temp=[5 10 20 90];
for i=1:4
    n=temp(i);
    x=GMRES(A,b,n);
    e(i)=norm(x-x0);
end
for j=1:5
    v=v_+rand(100,1)*10^(j-5);
    A=diag(v);
    A=P*A*P^-1;
    x0=A\b;
    for n=1:22
    x=GMRES(A,b,n);
    e_(n,j)=norm(x-x0);
    end
end
figure
for i=1:5
    plot(1:22,log(e_(:,i)))
    hold on
end
% e(1) denotes the error of solution for n=5
% e(2) denotes the error of solution for n=10
% e(3) denotes the error of solution for n=20
% e(4) denotes the error of solution for n=90
xlabel('# of iterations')
ylabel('log(error)')
legend('¦Å=10^-4','¦Å=10^-3','¦Å=10^-2','¦Å=10^-1','¦Å=10^0','location','northeast')