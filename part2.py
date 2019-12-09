
import numpy as np
from numpy.linalg import norm
from numpy.linalg import inv
from numpy.linalg import solve
def GMRES(A,b,n):
    Q = np.zeros((100,n+1))
    S = np.zeros((n+1, n))
    t = np.zeros((n,1))
    kk = b/norm(b)
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
    return(x)


v = np.zeros((100,1))
e = np.zeros((4,1))
N = np.array([[4,5,10]])
e_ = np.zeros((20,5))
for i in range(100):
    r=np.random.permutation(range(3))
    v[i]=N[:,r[0]]
v_ = v
v = np.random.rand(100, 1)+v
A = np.diag(v.flatten())
P = np.random.rand(100,100)
A = P*A*inv(P)
b = np.random.rand(100, 1)
b=b.flatten()
x0 = np.dot(np.dot(inv(np.dot(A.T,A)),A.T),b)
temp = np.array([[5,10,20,90]])
e=[]
for i in range(1):
    n = temp[:,i]
    n = n[0]
    x = GMRES(A,b,n)
    e.append(norm(x-x0))
e
