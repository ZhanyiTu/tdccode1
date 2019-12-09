from numpy import *;
import numpy as np;
from numpy.linalg import norm
from numpy.linalg import inv

# v = mat(zeros())
# e = mat(random.randint(0, size(4, 1)))
# N = mat([4, 5, 10])
# e_ = mat(random.randint(0, size(20, 5)))
# for i in range(100):
#     r = random.permutation(range(3))
#     v[i] = N[:, r[0]]
def GMRES(A,b,n):
    Q = mat(np.zeros((100,n+1)))
    S = mat(np.zeros((n+1, n)))
    t = mat(np.zeros((n,1)))
    Q[:,0] = b/norm(b)
    for i in np.arange(1, n + 1, 1):
        Q[:, i] = A * Q[:, i - 1]
        for j in np.arange(0, i , 1):
            S[j, i - 1] = Q[:, i].T * Q[:, j]
        Q[:, i] = Q[:, i] - Q[:, 0:(i - 1)] * S[0:(i - 1), i - 1]
        S[i, i - 1] = norm(Q[:, i])
        Q[:, i] = Q[:, i] / S[i, i - 1]
    Q = Q[:,0:n ]
    bb = mat(np.zeros((n + 1,1))) #第一个数组为norm(b) 其余n - 1个数组为0的数组
    bb[0, 0] = norm(b)
    t = inv(S.T * S) * S.T * bb#S最中间1列不应该是0 Q应该是正常的 S现在正常了 Q一开始少了一行 现在也正常了 bb少了一行
    x = Q * t
    return(x)


v = mat(np.zeros((100,1)))#100个长度为1 内容为0的数组
e = mat(np.zeros((4,1)))#4个长度为1 内容为0的数组
N = mat(np.array([[4,5,10]]))#n为数组4,5,10
e_ = mat(np.zeros((22,5))) #20个长度为5 内容为5个0的数组
for i in range(100):
    r=np.random.permutation(range(3))#随机排列序列 range(3)=[0, 1, 2]
    v[i]=N[:,r[0]]#v每项随机在N里选
v_ = v # v_是100个5 4 10 排列组合的1维数组

v = mat(np.random.rand(100, 1))+v#rand: 生成100个大小为1的数组每个值小于1大于等于0 v是v_加上0到1的随机数 形式不变
A = mat(np.diag(v.A.flatten()))#flatten 返回成1维的数组或者用1维数组生成对角阵 diag:返回对角线数组 就是100*100的矩阵, 只有中间有数 A是正确的 其余为0中间有数
P = mat(np.random.rand(100,100))#生成100个大小为100 每个值在[0,1)的数组 P是每个值都在0 1之间的随机数 是100100数组 正确
A = (P*A*P).I#inv计算逆 这样的乘法计算是错误的 p是正常的 原来是先乘再求逆 服了.
b = mat(np.random.rand(100, 1))
x0 = A.I * b#A 出现了非常大的数 说明前面有错误 A正常 b正常 A\b A左除b A的逆乘b
temp = [5,10,20,90]
for i in range(4):
    n = temp[i]
    x = GMRES(A, b, n)
    e[i] = norm(x - x0)
for j in range(5):
    v = v_ + np.random.rand(100, 1) * pow(10, j-4)
    A = np.diag(v.A.flatten())
    A = P * A * P.I
    x0 = A.I * b
    for n in arange(1, 23, 1):
        x = GMRES(A, b, n)
        e_[n - 1, j] = norm(x - x0)#x0的计算有问题
print(e_)