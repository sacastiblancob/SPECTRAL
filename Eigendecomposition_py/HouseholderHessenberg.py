import numpy as np

A=np.array([[12,3.1,-3.1,1],[-2,-8.4,3,2],[1,2,17.5,2],[1,2.1,3.8,8.3]])
#A=np.array([[12,-2,1,1],[-2,-8.4,3,2.1],[1,3,17.5,3.8],[1,2.1,3.8,8.3]])

n=np.shape(A)[1]
P_tot=np.identity(n)
H=np.copy(A)
for i in range(1,n-1):
    u=np.array([np.copy(H[i:,i-1])]).T
    u[0]=u[0]-np.linalg.norm(u)
    u=u/np.linalg.norm(u)

    P=np.identity(n)
    P[i:,i:]=np.identity(n-i)-2*u@u.T
    P_tot=P@P_tot
    H=P@H@P

Lambdas_A,evectors_A=np.linalg.eig(A)
Lambdas_H,evectors_H=np.linalg.eig(H)
X = P_tot.T@evectors_H
