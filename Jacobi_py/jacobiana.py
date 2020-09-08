# -*- coding: utf-8 -*-
"""
SERGIO CASTIBLANCO
01-04-2020

THIS CODE IS MADE FOR ANALYSE THE PROPERTIES OF THE VANDERMONDE MATRICES OF 
JACOBI POLYNOMIALS BASIS

"""
import os
os.chdir("F:\SPECTRAL\Codes\Jacobi_py")

import numpy as np
import funspe as fs
#import matplotlib.pyplot as plt
#from matplotlib import style
from numpy import linalg as ln

#INITIAL INPUTS alpha, beta and N
#   N = polynomial order
#   if alpha=beta=0.0 --> Legendre Polynomials
#   if alpha=beta=1.0 --> Chebyshev Polynomials
#
# and remind that if you have an polynomial order N you have N+1 nodes

alpha = 0.0
beta = 0.0
N = 4
M = N+1     #number of nodes

#COMPUTING THE APPROPIATE X AND WEIGHTS (GAUSSIAN QUADRATURE) FOR THE JACOBI BASIS

#SELECT IF YOU WANT THE DOMAIN OPEN (-1,1) OR CLOSED [-1,1]. 1 FOR OPEN, 2 FOR CLOSE
dom = 1
if (dom==1):
    x,w = fs.jacobigq(alpha,beta,N)
else:
    x = fs.jacobigl(alpha,beta,N)
    xdummy,w = fs.jacobigq(alpha+1.0,beta+1.0,N-2)
    w = np.concatenate(([0.1],w,[0.1]),axis=None)

#COMPUTING THE VANDERMONDE MATRIX P AND HIS RESPECTIVE INVERSE
P = fs.jacobip(x,alpha,beta,N).T
PI = ln.inv(P)

#codigo de burgers2D, w para N=4 es el siguiente
w2 = np.array([0.10000000000000001,0.54444444444444429,0.71111111111111114,0.54444444444444429,0.10000000000000001])
#and in matrix form is
W2 = np.diag(w2)
#in the same code of Burgers2D, theres a matrix B that is the Legendre polynomials
#evaluated in x, is the same as P orthogonal
B = np.array([[0.50000000000000000,0.50000000000000000,0.50000000000000000,0.50000000000000000,0.50000000000000000],[-1.5000000000000000,-0.98198050606196563,0.0000000000000000,0.98198050606196563,1.5000000000000000],[2.5000000000000000,0.35714285714285698,-1.2500000000000000,0.35714285714285698,2.5000000000000000],[-3.5000000000000000,0.98198050606196596,0.0000000000000000,-0.98198050606196596,3.5000000000000000],[2.5000000000000000,-1.0714285714285716,0.93750000000000000,-1.0714285714285716,2.5000000000000000]])
#then appears a matrix C that is
c = np.array([0.5,1.5,2.5,3.5,2.5])
C = np.diag(c)
#those C are made in the next way
#where n is equal to the number of points that we are ussing, in this case n=5
#   C = 0.0
# 	W = 0.0
# 	L = 0.0
# 	do i = 1, n - 1
# 	 C(i, i) = (i - 1.) + 0.5
# 	 W(i, i) = wg(i)
# 	enddo
# 	C(n,n) = n / 2.
#then a matrix M is computed in the next way:
#   M = matmul(C, transpose(B))         --> This, I belive, is some kind of normalization
#                                           Im gonna call this matrix B2
#   M = matmul(M, W)                    --> This is the matrix B*W
B2 = C@(B.T)
M = B2@W2
U3 = B2.T@W2@B2
U4 = B@W2@B.T
#some of those U3 or U4 must be the Identity

# proof of P.T*W*P = I
W = np.diag(w)
U = (P.T)@W@P
U5 = (P.T)@W2@P
RW = np.sqrt(W)
U2 = RW@P
I = (U2.T)@U2

#the W that must be
W3 = ln.inv(P.T)@np.identity(5)@ln.inv(P)
W5 = ln.inv(B)@np.identity(5)@ln.inv(B.T)
W6 = ln.inv(B2.T)@np.identity(5)@ln.inv(B2)

Z = ln.inv(P)@P

#PROPERTIES OF THE VANDERMONDE MATRIX P AND HIS INVERSE

#determinant
Pdet = ln.det(P)

#condition number
Pcond = ln.cond(P)

#rank
Prank = ln.matrix_rank(P)

#QR decomposition
Q,R = ln.qr(P)

#eigendecomposition
L,X = ln.eig(P)

#svd
S,V,Dh = ln.svd(P)
D = Dh.T












