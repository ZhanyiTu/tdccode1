from numpy import *
from numpy.linalg import *

def ArnoldiMethod(A,n):
    m = len(A)
    b = random.rand(m)
    b = transpose([b])
    S = zeros([n+1,n])
    Q = zeros([m,n+1])
    Q[:,0] = list(b/norm(b))
    for i in range(2,n+2):
        Q[:,i-1] = list(transpose([dot(A, Q[:,i-2])]))
        for j in range(1,i):
            S[j-1,i-2] = dot(Q[:,i-1],transpose([Q[:,j-1]]))
            print("j=",j,S)
        Q[:,i-1] = list(transpose([Q[:,i-1]]) - dot(Q[:,0:i-1],transpose([S[0:i-1,i-2]])))
        S[i-1,i-2] = norm(Q[:,i-1])
        Q[:,i-1] = Q[:,i-1]/S[i-1,i-2]
    S = S[0:n,:]
    Q = Q[:,0:n]
    return S,Q

A = np.random.random((100, 100))
e1, v1 = np.linalg.eig(A)
sigma1=np.linalg.norm(e1,ord=np.Inf,axis=0)
error=zeros([10,1])
for n in range(0,10):
    Q,S=ArnoldiMethod(A,n)
    en, vn = np.linalg.eig(A)
    sigma2=np.linalg.norm(en,ord=np.Inf,axis=0)
    error[n]=sigma1-sigma2