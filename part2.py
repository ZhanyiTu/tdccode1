
import numpy as np
from numpy.linalg import norm
from numpy.linalg import inv
from numpy.linalg import solve
def GMRES(A,b,n):
    Q = np.zeros((100,n+1))
    S = np.zeros((n+1, n))
    t = np.zeros((n,1))
    kk = b/norm(b) #相当于Q
    kk = kk.flatten()
    Q[:,0] = kk
    for i in np.arange(1, n+1, 1):
        Q[:,i]=np.dot(A,Q[:,i-1])
        for j in np.arange(0, i, 1):
            S[j,i-1]=np.dot(Q[:,i].T,Q[:,j])
        Q[:,i]=Q[:,i]-np.dot(Q[:,0:(i-1)],S[0:(i-1),i-1])
        S[i,i-1]=norm(Q[:,i])
    Q = Q[:,0:n]
    nm = np.array([[norm([b])]])
    zs = np.zeros((1,n))
    b = np.hstack((nm,zs))
    b = b.T
    t = np.dot(np.dot(inv(np.dot(S.T,S)),S.T),b)
    x = np.dot(Q,t)
    x = x.flatten();
    return(x)


v = np.zeros((100,1))#100个长度为1 内容为0的数组
e = np.zeros((4,1))#4个长度为1 内容为0的数组
N = np.array([[4,5,10]])#n为数组4,5,10
e_ = np.zeros((20,5)) #20个长度为5 内容为5个0的数组
for i in range(100):
    r=np.random.permutation(range(3))#随机排列序列 range(3)=[0, 1, 2]
    v[i]=N[:,r[0]]#v每项随机在N里选
v_ = v
v = np.random.rand(100, 1)+v#rand: 生成100个大小为1的数组每个值小于1大于等于0
A = np.diag(v.flatten())#flatten 返回成1维的数组 diag:返回对角线数组 就是100*100的矩阵, 只有中间有数
P = np.random.rand(100,100)#生成100个大小为100 每个值在[0,1)的数组
A = P*A*inv(P)#inv计算逆
b = np.random.rand(100, 1)
b=b.flatten()
x0 = np.dot(np.dot(inv(np.dot(A.T,A)),A.T),b)#两个数组的点积
temp = np.array([[5,10,20,90]])
e=[]
for i in range(4):
    n = temp[:,i]
    n = n[0]
    x = GMRES(A,b,n)
    temp2 = norm(x-x0)
    e.append(temp2)
for j in range(5)
    v =
'''
function x=GMRES(A,b,n)//定义一个叫做gmres的函数
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
    b=[norm(b) zeros(1,n)]';
    t=inv(S'*S)*S'*b;
    x=Q*t;
end


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
legend('ε=10^-4','ε=10^-3','ε=10^-2','ε=10^-1','ε=10^0','location','northeast')

>> u=[10^-4 10^-3 10^-2 10^-1 10^0];
>> v=[8 9 11 15 18];
>> plot(log(u),v);
>> xlabel('log(ε)');
>> ylabel('# of iterations for the convergence');
>> title('choice of ε affects the speed of convergence')
'''
