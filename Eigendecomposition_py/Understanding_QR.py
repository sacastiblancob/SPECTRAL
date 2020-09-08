# -*- coding: utf-8 -*-
"""
UNDERSTANDING QR TRANSFORMATION AND EIGENDECOMPOSITION

SERGIO CASTIBLANCO - 04/2020
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as ln

#DEFINING SOME MATRIX A
A = np.array([[2,1],[1,2]])

#DEFINING THE UNITARY MATRIX
U = np.identity(2)

#DEFINING THE ORIGIN
origin = [0], [0]

#PLOTING THE TRANSFORMATION A
plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['b','g'], angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, U[0,:], U[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#QR DECOMPOSITION
Q,R = ln.qr(A)
Q1 = Q
R1 = R
Q2 = Q@np.array([[-1,0],[0,1]])
R2 = Q2.T@A

#lets understanD what a QR decomposition means, ploting U,Q,R
plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q1[0,:], Q1[1,:], color = ['b','b'], angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R1[0,:], R1[1,:], color = ['r','r'], angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['b','b'], angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R2[0,:], R2[1,:], color = ['r','r'], angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines


#second version of Q
plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#recall from where comes Q
#All the QR decomposition

plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

A1 = np.array([[Q2[0,0],A[0,1]],[Q2[1,0],A[1,1]]])
plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A1[0,:], A1[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

A1 = np.array([[Q1[0,0],A[0,1]],[Q1[1,0],A[1,1]]])
plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A1[0,:], A1[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

A1 = np.array([[Q2[0,0],A[0,1]],[Q2[1,0],A[1,1]]])
plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A1[0,:], A1[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

A1 = np.array([[Q2[0,0],A[0,1]],[Q2[1,0],A[1,1]]])
b = np.dot(Q2[:,0],A[:,1])*Q2[:,0]
A2 = np.array([[Q2[0,0],b[0]],[Q2[1,0],b[1]]])

plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A2[0,:], A2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

A2 = np.array([[Q2[0,0],b[0]],[Q2[1,0],b[1]]])
b = A[:,1]-np.dot(Q2[:,0],A[:,1])*Q2[:,0]
A3 = np.array([[Q2[0,0],b[0]],[Q2[1,0],b[1]]])


plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A3[0,:], A3[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines


plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

##ones all QR decomposition
#lets aply it, first R then Q

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R2[0,:], R2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

R21 = np.array([[R2[0,0],U[0,1]],[R2[1,0],U[1,1]]])

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R21[0,:], R21[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

R22 = np.array([[R2[0,0],U[0,1]],[R2[1,0],U[1,1]]])
b[0] = R2[0,1]
b[1] = U[1,1]
R22 = np.array([[R22[0,0],b[0]],[R22[1,0],b[1]]])


plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R22[0,:], R22[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R2[0,:], R2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#THEN APPLY Q

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines


plt.clf()
plt.quiver(*origin, R2[0,:], R2[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A[0,:], A[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#WHAT IS RQ, FIRST Q AND THEN R
#BECAUSE THEN A = RQ

#first the rotation Q
plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#then the squishing R

Q3 = np.array([[R2[0,0],U[0,1]],[R2[1,0],U[1,1]]])
Q3 = Q3@Q2

plt.clf()
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q3[0,:], Q3[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

Q4 = np.array([[U[0,0],R2[0,1]],[U[1,0],U[1,1]]])
Q4 = Q4@Q3

plt.clf()
plt.quiver(*origin, Q3[0,:], Q3[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q4[0,:], Q4[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

Q5 = np.array([[U[0,0],U[0,1]],[U[1,0],R2[1,1]]])
Q5 = Q5@Q4

plt.clf()
plt.quiver(*origin, Q4[0,:], Q4[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q5[0,:], Q5[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines


A10 = R2@Q2
plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A10[0,:], A10[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#recollecting everything

#HOW TO GET A FROM QR
plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R2[0,:], R2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, R2[0,:], R2[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A[0,:], A[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#HOW TO GET A2 FROM RQ
plt.clf()
plt.quiver(*origin, A10[0,:], A10[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, U[0,:], U[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

plt.clf()
plt.quiver(*origin, R2[0,:], R2[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A10[0,:], A10[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#SO WHAT IS THE QR OF A2
plt.clf()
plt.quiver(*origin, A10[0,:], A10[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#DRAW QR OF A2

#MAKING ALL THE ITERATIONS
plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A10[0,:], A10[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#QR decomposition of A2
Q21,R21 = ln.qr(A10)
Q21 = Q21@np.array([[-1,0],[0,1]])
R21 = Q21.T@A10

plt.clf()
plt.quiver(*origin, A10[0,:], A10[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q21[0,:], Q21[1,:], color = ['b','b'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R21[0,:], R21[1,:], color = ['g','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#computing A3
A20 = R21@Q21
plt.clf()
plt.quiver(*origin, A10[0,:], A10[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A20[0,:], A20[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#QR decomposition of A3
Q31,R31 = ln.qr(A20)
Q31 = Q31@np.array([[-1,0],[0,1]])
R31 = Q31.T@A20

plt.clf()
plt.quiver(*origin, A20[0,:], A20[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q31[0,:], Q31[1,:], color = ['b','b'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R31[0,:], R31[1,:], color = ['g','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#computing A4
A30 = R31@Q31
plt.clf()
plt.quiver(*origin, A20[0,:], A20[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A30[0,:], A30[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#QR decomposition of A4
Q41,R41 = ln.qr(A30)
Q41 = Q41@np.array([[-1,0],[0,1]])
R41 = Q41.T@A30

plt.clf()
plt.quiver(*origin, A30[0,:], A30[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q41[0,:], Q41[1,:], color = ['b','b'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R41[0,:], R41[1,:], color = ['g','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#computing A5
A40 = R41@Q41
plt.clf()
plt.quiver(*origin, A30[0,:], A30[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A40[0,:], A40[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#QR decomposition of A5
Q51,R51 = ln.qr(A40)
Q51 = Q51@np.array([[-1,0],[0,1]])
R51 = Q51.T@A40


#plotting every A's
plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A10[0,:], A10[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A20[0,:], A20[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A30[0,:], A30[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, A40[0,:], A40[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#plotting every Q's
plt.clf()
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q21[0,:], Q21[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q31[0,:], Q31[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q41[0,:], Q41[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q51[0,:], Q51[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#plotting every R's
plt.clf()
plt.quiver(*origin, R2[0,:], R2[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R21[0,:], R21[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R31[0,:], R31[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R41[0,:], R41[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, R51[0,:], R51[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#plotting A and autovectors of A
X = np.array([[1,-1],[1,1]])
X = X*(1/np.sqrt(2))
plt.clf()
plt.quiver(*origin, A[0,:], A[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, X[0,:], X[1,:], color = ['b','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines


E1,V1 = ln.eig(A)

plt.clf()
plt.quiver(*origin, Q2[0,:], Q2[1,:], color = ['b','b'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q21[0,:], Q21[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q31[0,:], Q31[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q41[0,:], Q41[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, Q51[0,:], Q51[1,:], color = ['gray','gray'],angles='xy', scale_units='xy', scale=1)
plt.quiver(*origin, V1[0,:], V1[1,:], color = ['g','g'],angles='xy', scale_units='xy', scale=1)
plt.axis('equal')  #<-- set the axes to the same scale
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.grid(b=True, which='major') #<-- plot grid lines

#AUTOVECTORES
H = Q51@Q41@Q31@Q21@Q2
#H = Q21@Q2
#H = Q31@H






