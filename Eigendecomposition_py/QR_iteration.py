import numpy as np
from numpy import linalg as ln
from scipy import linalg as lns
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import


A=np.array([[-1,4,3,1],[2,-5,3,2],[2,2,-3,4],[0,1,3,-4]])
Lambdas,evectors=np.linalg.eig(A)
A_orig=A
n=100
H = np.identity(A.shape[0])
J = np.identity(A.shape[0])
for i in range(0,n):
    Q,R=ln.qr(A)
#    Q = Q@np.array([[-1,0,0],[0,-1,0],[0,0,1]])
#    R = Q.T@A
    H = Q.T@H
    A=R@Q
    
    error=(np.array(sorted(np.diag(A),key=abs))-
           np.array(sorted(Lambdas,key=abs)))
H = H.T    
E = H@(R@Q)@H.T

origin = [0], [0], [0]
# plt.clf()
fig = plt.figure()
ax = fig.gca(projection=r'3d')
ax.quiver(*origin, A[0,:], A[1,:], A[2,:], color = ['r'], arrow_length_ratio=0.3)
# ax.axis('equal')  #<-- set the axes to the same scale
ax.set_xlim([-6,6])
ax.set_ylim([-6,6])
ax.set_zlim([-6,6])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.grid(b=True, which='major') #<-- plot grid lines
fig.show()

