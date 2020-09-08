# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 16:36:31 2020

@author: SERGIO
"""


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
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt
#from matplotlib import style
from numpy import linalg as ln
import time

#INITIAL INPUTS alpha, beta and N
#   N = polynomial order
#   if alpha=beta=0.0 --> Legendre Polynomials
#   if alpha=beta=1.0 --> Chebyshev Polynomials
#
# and remind that if you have an polynomial order N you have N+1 nodes

alpha = 0.0
beta = 0.0
n = 4
m = n+1     #number of nodes

x = fs.jacobigl(alpha,beta,n)

#COMPUTING THE VANDERMONDE MATRIX P
B = fs.jacobip(x,alpha,beta,n).T
B3 = fs.legendre(x,n).T

#w = np.array([0.10000000000000001,0.54444444444444429,0.71111111111111114,0.54444444444444429,0.10000000000000001])
#w = np.array([1/10,49/90,32/45,49/90,1/10])
w = 2/(m*(m-1)*(B3[:,n]**2))
W = np.diag(w)

#computing u
un = x**2
um = ln.solve(B,un)

left = B.T@W@B@um
rig = B.T@W@un
#indeed left is equal to rigth

left3 = B3.T@W@B3@um
rig3 = B3.T@W@un
#neither in this case left is not equal to rigth

#seems that only if the matrix B comes from the orthonormal version of Legendre, that equality is true
#maybe the orthogonal basis are influenced by their norms

#computing the integral between -1 and 1
INT = np.dot(un,w)

#unitary
U = B.T@W@B
C = ln.inv(U)
U3 = B3.T@W@B3
C3 = ln.inv(U3)

##############################################################################
# localFil.f from Burgers Code (J. Escobar)
##############################################################################

# Matrix B, same as the orthogonal version in this code
B2 = B3

#then appears a matrix C that is
#c2 = np.array([0.5,1.5,2.5,3.5,2.5]) --> c2 for N=5 (number of points), extracted from Burgers2D codes
c2 = np.arange(1,m+1)
c2 = (c2-1)+0.5
c2[n] = m/2
C2 = np.diag(c2)
#those C are made in the next way
#where n is equal to the number of points that we are ussing, in this case n=5
#   C = 0.0
# 	W = 0.0
# 	L = 0.0
# 	do i = 1, n - 1
# 	 C(i, i) = (i - 1.) + 0.5
# 	 W(i, i) = wg(i)
# 	enddo
# 	C(n,n) = n / 2.      #this seems an error, becaus n must be the grade of the polynomial not number of points
#   W(n,n) = wg(n)
#	M = matmul(C, transpose(B))

#The last 2.5 should be 2.0 becuase N/2 = 4/2 --> filter.pdf
#The matrix W is diag(GLL weights)
# M = inv(B.T@W@B)@B.T@W
# M = C@B.T@W   ;   C = inv(B.T@W@B)   --> filter.pdf

M = C@B.T@W   #--> computing M with B from orthonormal basis
M2 = C2@B2.T@W              #--> computing M with B from orthogonal basis and C of Burgers
M3 = C3@B2.T@W               #--> computing M with B from orthonormal and C well computed from filter.pdf 

#then the matrix L is computed, with the filter definition
pf = 50
alp = -np.log(1.0e-16)       # Machine precision according to Diamessis 2008

l = np.arange(1,m+1)
l = np.exp(-alp * ((l - 1.0) / (m - 1.0)) ** pf)

L = np.diag(l)

#finally F = B@L@M

F = B@L@M                   #--> computing M with B from orthonormal basis
F2 = B2@L@M2                #--> computing M with B from orthogonal basis and C of Burgers
F3 = B2@L@M3                #--> computing M with B from orthonormal and C well computed from filter.pdf
#all F are the same

##############################################################################
# Sumedh Thesis
##############################################################################

U4 = np.sqrt(W)@B
I1 = ln.inv(U4)
I2 = U4.T
#I1 and I2 are almost the same, but we have the problem that the numerical integral is not perfect
#but, let's continue with the Sumedh proposition

##############################################################################
# PROPERTIES ANALYISIS
##############################################################################

############### not so good numerical integration ##########################

# plot of the las term in the numerical integration
dim1 = 100
alpha = 0
beta = 0
H=np.zeros((4,dim1-3))
for i in range(4,dim1+1):
    
    xi = fs.jacobigl(alpha,beta,i)
    P1 = fs.jacobip(xi,alpha,beta,i).T
    P2 = fs.legendre(xi,i).T
    
    j = i+1
    w1 = 2/(j*(j-1)*(P2[:,i]**2))
    W1 = np.diag(w1)
    RW = np.diag(np.sqrt(w1))
    
    I1 = P1.T@W1@P1
    I2 = ln.inv(I1)
    
    H[0,i-4] = i
    H[1,i-4] = I1[i,i]
    H[2,i-4] = I2[i,i]
#    H[3,i-4] = 1/(2+(1/i))
    H[3,i-4] = i/(2*i+1)

    
plt.plot(H[0,:], H[1,:],label=("Last of P.T@W@P"))
plt.plot(H[0,:], H[2,:],label=("Last of inv(P.T@W@P)"))
plt.ylabel(r'Condition Number')
plt.xlabel(r'Matrix Dimension')
plt.title('Condition numbers')
plt.legend(loc='upper left')

###############condition numbers##########################

#plot of condition number of vandermonde matrices
dim1 = 100
alpha = 0
beta = 0
H=np.zeros((3,dim1-3))
for i in range(4,dim1+1):
    
    xi = fs.jacobigl(alpha,beta,i)
    P1 = fs.jacobip(xi,alpha,beta,i).T
    P2 = fs.legendre(xi,i).T
    
    H[0,i-4] = i
    H[1,i-4] = ln.cond(P1)
    H[2,i-4] = ln.cond(P2)

    
plt.semilogy(H[0,:], H[1,:],label=("Orthonormal P"))
plt.semilogy(H[0,:], H[2,:],label=("Orthogonal P"))
plt.ylabel(r'Condition Number')
plt.xlabel(r'Matrix Dimension')
plt.title('Condition numbers')
plt.legend(loc='upper left')

#plot of condition number of gauss weigths times vandermonde matrices


dim1 = 100
alpha = 0
beta = 0
H=np.zeros((4,dim1-3))
for i in range(4,dim1+1):
    
    xi = fs.jacobigl(alpha,beta,i)
    P1 = fs.jacobip(xi,alpha,beta,i).T
    P2 = fs.legendre(xi,i).T
    
    j = i+1
    w1 = 2/(j*(j-1)*(P2[:,i]**2))
    W1 = np.diag(w1)
    RW = np.diag(np.sqrt(w1))
    
    # r1 = (2/(j*(j-1)))*(P1[0,i]*P1[0,i]+P1[i,i]*P1[i,i])
    # r = (-j*((j-1)**3)*(2**(2*j-1))*((np.math.factorial(j-2))**4))/((2*j-1)*((np.math.factorial(2*j-2))**3))
    # r = r*P1[:,i]*P1[:,i]

    
    WP1 = W1@P1
    WP2 = W1@P2
    RWP = RW@P1
    
    H[0,i-4] = i
    H[1,i-4] = ln.cond(WP1)
    H[2,i-4] = ln.cond(WP2)
    H[3,i-4] = ln.cond(RWP)
    
plt.semilogy(H[0,:], H[1,:],label=("W*Orthonormal P"))
plt.semilogy(H[0,:], H[2,:],label=("W*Orthogonal P"))
plt.semilogy(H[0,:], H[3,:],label=("sqrt(W)*Orthonormal P"))
plt.ylabel(r'Condition Number')
plt.xlabel(r'Matrix Dimension')
plt.title('Condition numbers')
plt.legend(loc='upper left')

#plot of condition number of matrices P.T@W@P, filter.pdf equation 11
dim1 = 100
alpha = 0
beta = 0
H=np.zeros((4,dim1-3))
for i in range(4,dim1+1):
    
    xi = fs.jacobigl(alpha,beta,i)
    P1 = fs.jacobip(xi,alpha,beta,i).T
    P2 = fs.legendre(xi,i).T
    
    j = i+1
    w1 = 2/(j*(j-1)*(P2[:,i]**2))
    W1 = np.diag(w1)
    RW = np.diag(np.sqrt(w1))
    
    MP1 = P1.T@W1@P1
    MP2 = P2.T@W1@P2
    
    H[0,i-4] = i
    H[1,i-4] = ln.cond(MP1)
    H[2,i-4] = ln.cond(MP2)
    
plt.semilogy(H[0,:], H[1,:],label=("W*Orthonormal P.T"))
plt.semilogy(H[0,:], H[2,:],label=("W*Orthogonal P.T"))
plt.ylabel(r'Condition Number')
plt.xlabel(r'Matrix Dimension')
plt.title('Condition numbers')
plt.legend(loc='upper left')

#plot of condition number of matrices M = inv(P.T@W@P)@P.T@W, filter.pdf equation 11
dim1 = 100
alpha = 0
beta = 0
H=np.zeros((4,dim1-3))
for i in range(4,dim1+1):
    
    xi = fs.jacobigl(alpha,beta,i)
    P1 = fs.jacobip(xi,alpha,beta,i).T
    P2 = fs.legendre(xi,i).T
    
    
    c1 = np.ones(i+1)
    c1[i] = i/(2*i+1)
    C1 = np.diag(c1)
    
    c2 =  np.arange(0,i+1) + 1/2
    c2[i] = i/2
    C2 = np.diag(c2)
    
    j = i+1
    w1 = 2/(j*(j-1)*(P2[:,i]**2))
    W1 = np.diag(w1)
    RW = np.diag(np.sqrt(w1))
    
    MP1 = C1@P1.T@W1
    MP2 = C2@P2.T@W1
    MP3 = P1.T@W1
    
    H[0,i-4] = i
    H[1,i-4] = ln.cond(MP1)
    H[2,i-4] = ln.cond(MP2)
    H[3,i-4] = ln.cond(MP3)
    
plt.semilogy(H[0,:], H[1,:],label=("M from orthonormal"))
plt.semilogy(H[0,:], H[2,:],label=("M from orthogonal"))
plt.semilogy(H[0,:], H[3,:],label=("M from orthonormal if integration exact"))
plt.ylabel(r'Condition Number')
plt.xlabel(r'Matrix Dimension')
plt.title('Condition numbers')
plt.legend(loc='upper left')


#plot of condition number of matrices F = B@L@M, filter.pdf equation 11
pf = 5
alp = -np.log(1.0e-16)       # Machine precision according to Diamesis 2008

dim1 = 100
alpha = 0
beta = 0
H=np.zeros((4,dim1-3))
for i in range(4,dim1+1):
    
    xi = fs.jacobigl(alpha,beta,i)
    P1 = fs.jacobip(xi,alpha,beta,i).T
    P2 = fs.legendre(xi,i).T
    
    
    c1 = np.ones(i+1)
    c1[i] = i/(2*i+1)
    C1 = np.diag(c1)
    
    c2 =  np.arange(0,i+1) + 1/2
    c2[i] = i/2
    C2 = np.diag(c2)

    j = i+1    
    l = np.arange(1,j+1)
    l = np.exp(-alp * ((l - 1.0) / (i - 1.0)) ** pf)
    L = np.diag(l)
    
    w1 = 2/(j*(j-1)*(P2[:,i]**2))
    W1 = np.diag(w1)
    RW = np.diag(np.sqrt(w1))
    
    FP1 = P1@L@C1@P1.T@W1
    FP2 = P2@L@C2@P2.T@W1
    FP3 = P1@L@P1.T@W1
    
    H[0,i-4] = i
    H[1,i-4] = ln.cond(FP1)
    H[2,i-4] = ln.cond(FP2)
    H[3,i-4] = ln.cond(FP3)
    
plt.semilogy(H[0,:], H[1,:],label=("F from orthonormal"))
plt.semilogy(H[0,:], H[2,:],label=("F from orthogonal"))
plt.semilogy(H[0,:], H[3,:],label=("F from orthonormal if integration exact"))
plt.ylabel(r'Condition Number')
plt.xlabel(r'Matrix Dimension')
plt.title('Condition numbers')
plt.legend(loc='upper left')

#FILTER KILLS THE SPACE IN FACT
#TAKE THE RANK OF THE MATRIX AND GET LOWER, THE FILTER MATRIX HAVE AN ACTUAL NULLSPACE
    
############### descomposition an things ##########################
i = 10
j = i+1

#vandermonde matrix
xi = fs.jacobigl(alpha,beta,i)
P1 = fs.jacobip(xi,alpha,beta,i).T
P2 = fs.legendre(xi,i).T

#c for inv(P.T@W@P)
c1 = np.ones(i+1)
c1[i] = i/(2*i+1)
C1 = np.diag(c1)

c2 =  np.arange(0,i+1) + 1/2
c2[i] = i/2
C2 = np.diag(c2)

#filter function matrix
pf = 50
alp = -np.log(1.0e-16)       # Machine precision according to Diamesis 2008

j = i+1    
l = np.arange(1,j+1)
l = np.exp(-alp * ((l - 1.0) / (i - 1.0)) ** pf)
L = np.diag(l)

#gll weights
w1 = 2/(j*(j-1)*(P2[:,i]**2))
W1 = np.diag(w1)
RW = np.diag(np.sqrt(w1))

WP1 = W1@P1
WP2 = W1@P2
RWP = RW@P1

#Matrix M = C@P.T@W
MP1 = C1@P1.T@W1
MP2 = C2@P2.T@W1
MP3 = P1.T@W1

#Filter matrix P@L@M
FP1 = P1@L@C1@P1.T@W1
FP2 = P2@L@C2@P2.T@W1
FP3 = P1@L@P1.T@W1



#A for analysis porpouse
A1 = P1
A2 = P2

#trace
trace1 = np.trace(A1)
trace2 = np.trace(A2)

#rank
rank1 = ln.matrix_rank(A1)
rank2 = ln.matrix_rank(A2)

#determinant
det1 = ln.det(A1)
det2 = ln.det(A2)

# there's no nullspace, neither image nor nullity

#eigendecomposition
D1, V1 = ln.eig(A1)
D2, V2 = ln.eig(A2)

#singular value decomposition
U1, S1, Vh1 = ln.svd(A1)
V1 = Vh1.T
U2, S2, Vh2 = ln.svd(A2)
V2 = Vh2.T

#plotting the singular values
plt.semilogx(S1,np.zeros(j),'bs',label=("Singular values for A1"))
plt.semilogx(S2,np.ones(j),'rs',label=("Singular values for A2"))
plt.legend(loc='center left')
plt.title('Singular Values')

#plotting the singular spectrum
for i in range(0, i+1):
    plt.subplot(2, 2, 1)
    plt.plot(xi, U1[:,i], 'b')
    plt.ylim([-1.0, 1.0])
    plt.xlim([-1.05, 1.05])
    plt.ylabel(r'Left SVD, A1')
    plt.title('Left SVD, A1')
    
    plt.subplot(2, 2, 2)
    plt.plot(xi, U2[:,i], 'r')
    plt.ylim([-1.0, 1.0])
    plt.xlim([-1.05, 1.05])
    plt.ylabel(r'Left SVD, A2')
    plt.title('Left SVD, A2')
    
    plt.subplot(2, 2, 3)
    plt.plot(xi, V1[:,i], 'b')
    plt.ylim([-1.0, 1.0])
    plt.xlim([-1.05, 1.05])
    plt.ylabel(r'Rigth SVD, A1')
    plt.title('Rigth SVD, A1')
    
    plt.subplot(2, 2, 4)
    plt.plot(xi, V2[:,i], 'r')
    plt.ylim([-1.0, 1.0])
    plt.xlim([-1.05, 1.05])
    plt.ylabel(r'Rigth SVD, A2')
    plt.title('Rigth SVD, A2')
    
    plt.draw()    
    titulo = 'SPECTRUM'
    plt.suptitle(titulo)

    time.sleep(5)
    
    
plt.draw()
titulo = 'SPECTRUM'
plt.suptitle(titulo)























