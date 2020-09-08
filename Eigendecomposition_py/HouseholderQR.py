import numpy as np

A=np.array([[-1,-4,3],[-2,-5,-3],[2,2,3]])
R=np.copy(A)
n=np.shape(A)[1]
Q=np.identity(n)

for i in range(0,n-1):
    u=np.array([np.copy(R[i:,i])]).T
    u[0]=u[0]-np.linalg.norm(u)
    u=u/np.linalg.norm(u)
    
    H=np.identity(n)
    H[i:,i:]=np.identity(n-i)-2*u@u.T
    Q=H@Q
    
    R=H@R