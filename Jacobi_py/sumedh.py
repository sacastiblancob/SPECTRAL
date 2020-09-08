# -*- coding: utf-8 -*-
"""
SERGIO CASTIBLANCO
01-04-2020

THIS CODE IS MADE FOR ANALYSE THE PROPERTIES OF THE VANDERMONDE MATRICES OF 
JACOBI POLYNOMIALS BASIS

"""
import os
os.chdir("/media/aldair/FILES/SPECTRAL/Codes/Jacobi_py")

import numpy as np
import funspe as fs
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt
#from matplotlib import style
from numpy import linalg as ln
from matplotlib import style
import scipy.special as scs

# alpha = 0.0
# beta = 0.0
# n = 4
# m = n+1     #number of nodes

# x = fs.jacobigl(alpha,beta,n)

# #COMPUTING THE VANDERMONDE MATRIX P
# B = fs.jacobip(x,alpha,beta,n).T        #orthonormal vandermonde version
# B3 = fs.legendre(x,n).T                 #orthogonal vandermonde version

# #GLL w
# w = 2/(m*(m-1)*(B3[:,n]**2))
# W = np.diag(w)

# #this ought be the identity
# I = B.T@W@B

# #if I is identity then this ought be unitary
# U = np.sqrt(W)@B
# I2 = U@U.T
# I3 = U.T@U

# #trough QR
# Q,R = ln.qr(U)

# #guessing about unitarity
# U2 = np.copy(U)
# U2[:,n] = U2[:,n]/np.abs(R[n,n])
# #since R stores the "difference" of the matrix A with an unitary matrix
# #seem like the error in the numerical integration affect the n column of np.sqrt(W)@B

# Q2,R2 = ln.qr(U2)

# I4 = U2@U2.T

# plot of the "numerical integration residual"
# dim1 = 20
# alpha = 0
# beta = 0
# H=np.zeros((5,dim1-3))
# for i in range(4,dim1+1):
    
#     xi = fs.jacobigl(alpha,beta,i)
#     P1 = fs.jacobip(xi,alpha,beta,i).T
#     P2 = fs.legendre(xi,i).T
    
#     j = i+1
#     w1 = 2/(j*(j-1)*(P2[:,i]**2))
#     W1 = np.diag(w1)
#     RW = np.diag(np.sqrt(w1))
    
#     Q,R = ln.qr(RW@P1)
    
#     H[0,i-4] = i
#     H[1,i-4] = np.abs(R[i,i])**2
#     H[2,i-4] = ln.norm(Q)**2
#     H[3,i-4] = ln.norm(RW@P1)**2
#     H[4,i-4] = (2*i+1)/i

    
# plt.plot(H[0,:], H[1,:],label=("R[n,n]"))
# plt.plot(H[0,:], H[2,:],label=("normf(Q)"))
# plt.plot(H[0,:], H[3,:],label=("normf(RW@B)"))
# plt.ylabel(r'Condition Number')
# plt.xlabel(r'Matrix Dimension')
# plt.title('Condition numbers')
# plt.legend(loc='upper left')

#numerical reasons

# i = 7       #polynomial degree
# j = i+1     #number of nodes
# alpha = 0
# beta = 0

# #vandermonde matrix
# xi = fs.jacobigl(alpha,beta,i)
# P1 = fs.jacobip(xi,alpha,beta,i).T
# P2 = fs.legendre(xi,i).T

# #gll weights
# w1 = 2/(j*(j-1)*(P2[:,i]**2))
# W1 = np.diag(w1)
# RW = np.diag(np.sqrt(w1))

# #c for inv(P.T@W@P)
# c1 = np.ones(i+1)
# c1[i] = i/(2*i+1)
# C1 = np.diag(c1)

# c2 =  np.arange(0,i+1) + 1/2
# c2[i] = i/2
# C2 = np.diag(c2)

# #Sumedh alternative
# U1 = RW@P1
# U2 = RW@P1@(np.diag(np.sqrt(c1)))

# #filter function matrix

# ##CONSTRUCTION OF FILTERING MATRIX
# # OP=1 , EXPONENTIAL FILTER         exp(-alp * ((l - 1.0) / (i - 1.0)) ** pf)
# #                                   p=order
# # OP=2 , LANCZOS FILTER             (pi*x)^(-1) * sin(pi*x)  --> vandeven P=1
# # OP=3 , RAISED COSINE              [1 + cos(pi*x)]/2        --> vandeven P=2
# # OP=4 , SHARPENED RAISED COSINE    ([1 + cos(pi*x)]/2)^(4)*(35-84*([1 + cos(pi*x)]/2)^(1) + 70*([1 + cos(pi*x)]/2)^(2) - 20*([1 + cos(pi*x)]/2)^(3))
# #                                                            --> vandeven P=7   
# # OP=5 , ERFC FILTER                0.5*erfc(2*p^(0.5)*x)
# #                                   p=order
# #                                   0.5*erfc(2*p^(0.5)*(abs(x)-0.5)*sqrt((-log(1-4*(x-0.5)^2))/(4*(x-0.5)^2)))
# # OP=6 , LOG-ERFC FILTER            0.5*erfc(2*p^(0.5)*(abs(x)-0.5)*sqrt((-log(1-4*(x-0.5)^2))/(4*(x-0.5)^2)))
# #                                   p=order
# #For OP=1, OP=5 and OP=6 you should introduce order filter pf accordingly

# ####
# OP = 4
# ####

# alp = -np.log(1.0e-16)       # Machine precision according to Diamesis 2008
# pf = 30

# j = i+1    
# xauxf = np.arange(1,j+1)
# xauxf = ((xauxf - 1.0) / (j - 1.0))
# if OP==1:
#     l = np.exp(-alp * xauxf ** pf)
#     L = np.diag(l)
# elif OP==2:
#     l = (np.pi*xauxf)**(-1) * np.sin(np.pi*xauxf)
#     l[0] = 1.0
#     L = np.diag(l)
# elif OP==3:
#     l = (1+np.cos(np.pi*xauxf))/2
#     L = np.diag(l)
# elif OP==4:
#     l = (1+np.cos(np.pi*xauxf))/2
#     l = l**4*(35 - 84*l + 70*l**2 - 20*l**3)
#     L = np.diag(l)
# elif OP==5:
#     #l = 0.5*scs.erfc(2*pf**0.5*(np.abs(xauxf) - 0.5)*np.sqrt((-np.log(1-4*(xauxf-0.5)**2))/(4*(xauxf-0.5)**2)))
#     l = 0.5*scs.erfc(2*pf**0.5*(np.abs(xauxf)-0.5))
#     L = np.diag(l)
# elif OP==6:
#     l = 0.5*scs.erfc(2*pf**0.5*(np.abs(xauxf) - 0.5)*np.sqrt((-np.log(1-4*(xauxf-0.5)**2))/(4*(xauxf-0.5)**2)))
#     if (np.mod(i,2)==0):
#         l[int(i/2)]=0.5
#     L = np.diag(l)
# #what if L is the identity
# #L = np.identity(j)


# #M = inv[B.T@W@B]@B.T@W, lets proof the properties of B@W in all the possibilities
# WP1 = W1@P1
# WP2 = W1@P2
# RWP1 = U1
# RWP2 = U2

# #Matrix M = C@P.T@W
# MP1 = C1@P1.T@W1
# MP2 = C2@P2.T@W1

# #Filter matrix P@L@M
# F1 = P1@L@C1@P1.T@W1
# F2 = P2@L@C2@P2.T@W1
# #by this way, F1=F2
# #SVD of F1
# Q1,S1,Z1 = ln.svd(F1)
# #rank of F1
# r1 = ln.matrix_rank(F1)
# #condition number of F1, condition number has no sense since F is singular
# #con1 = ln.cond(F1)

# #filter matrix from Sumedh U@L@inv(U), in the unitary case inv(U)=U.T
# F3 = U1@L@ln.inv(U1)
# F4 = U2@L@U2.T
# #F5 = U1@L@U1.T, not even compute it, if L is de identity, F should be the identity also
# #that is not the case
# #by this way, F3=F4
# #SVD of F4
# Q2,S2,Z2 = ln.svd(F4)
# #rank of F4
# r2 = ln.matrix_rank(F4)
# #condition number of F4, condition number has no sense since F is singular
# #con2 = ln.cond(F4)



# #lets go trough the SVD
# QS1 = Q1@np.diag(S1)
# QS2 = Q2@np.diag(S2)

# o = 1

# plt.ion()
# plt.figure(1, figsize=(27, 15))
# style.use('ggplot')
# #plotting the singular spectrum
# for k in range(0, i+1):
# #    plt.clf()
    
#     plt.subplot(2, 8, o)
#     plt.plot(xi, QS1[:,k],'r',label=("Q@S"))
#     plt.ylim([-1.0, 1.0])
#     plt.xlim([-1.05, 1.05])
    
#     #plt.subplot(1, 2, 1)
#     plt.plot(xi, Z1.T[:,k],'b',label=("Z"))
#     plt.ylim([-1.0, 1.0])
#     plt.xlim([-1.05, 1.05])
#     plt.title('N, i='+str(k+1) + ' , \u03C3 = ' + str(round(S1[k],2)) )
    
#     if (k!=0 and k!=4):
#         plt.yticks(np.arange(-1, 1, step=0.25),[])
    
#     plt.subplot(2, 8, o+1)
#     plt.plot(xi, QS2[:,k],'r',label=("Q@S"))
#     plt.ylim([-1.0, 1.0])
#     plt.xlim([-1.05, 1.05])
    
# #    plt.subplot(2, 2, 4)
#     plt.plot(xi, Z2.T[:,k],'b',label=("Z"))
#     plt.ylim([-1.0, 1.0])
#     plt.xlim([-1.05, 1.05])
#     plt.title('U, i='+str(k+1) + ' , \u03C3 = ' + str(round(S2[k],2)) )
    
#     plt.yticks(np.arange(-1, 1, step=0.25),[])
    
#     plt.draw()    
#     titulo = 'Trough SVD'
#     plt.suptitle(titulo)
    
#     o = o+2
    
# #    input("Press Enter to continue...")
#     plt.pause(1)
    
# plt.draw()
# titulo = 'SPECTRUM'
# plt.suptitle(titulo)


# =============================================================================
# PROOF OVER A FUCNTION OF FILTER
# =============================================================================

i = 16
j = i+1
alpha = 0
beta = 0

#vandermonde matrix
xi = fs.jacobigl(alpha,beta,i)
P1 = fs.jacobip(xi,alpha,beta,i).T
P2 = fs.legendre(xi,i).T

#gll weights
w1 = 2/(j*(j-1)*(P2[:,i]**2))
W1 = np.diag(w1)
RW = np.diag(np.sqrt(w1))

#c for inv(P.T@W@P)
c1 = np.ones(i+1)
c1[i] = i/(2*i+1)
C1 = np.diag(c1)

c2 =  np.arange(0,i+1) + 1/2
c2[i] = i/2
C2 = np.diag(c2)

#Sumedh alternative
U1 = RW@P1
U2 = RW@P1@(np.diag(np.sqrt(c1)))

#filter function matrix

# #CONSTRUCTION OF FILTERING MATRIX
# OP=1 , EXPONENTIAL FILTER         exp(-alp * ((l - 1.0) / (i - 1.0)) ** pf)
#                                   p=order
# OP=2 , LANCZOS FILTER             (pi*x)^(-1) * sin(pi*x)  --> vandeven P=1
# OP=3 , RAISED COSINE              [1 + cos(pi*x)]/2        --> vandeven P=2
# OP=4 , SHARPENED RAISED COSINE    ([1 + cos(pi*x)]/2)^(4)*(35-84*([1 + cos(pi*x)]/2)^(1) + 70*([1 + cos(pi*x)]/2)^(2) - 20*([1 + cos(pi*x)]/2)^(3))
#                                                             --> vandeven P=7   
# OP=5 , ERFC FILTER                0.5*erfc(2*p^(0.5)*x)
#                                   p=order
#                                   0.5*erfc(2*p^(0.5)*(abs(x)-0.5)*sqrt((-log(1-4*(x-0.5)^2))/(4*(x-0.5)^2)))
# OP=6 , LOG-ERFC FILTER            0.5*erfc(2*p^(0.5)*(abs(x)-0.5)*sqrt((-log(1-4*(x-0.5)^2))/(4*(x-0.5)^2)))
#                                   p=order
# For OP=1, OP=5 and OP=6 you should introduce order filter pf accordingly

####
OP = 1
####

alp = -np.log(1.0e-16)       # Machine precision according to Diamesis 2008
pf = 15

j = i+1    
xauxf = np.arange(1,j+1)
xauxf = ((xauxf - 1.0) / (j - 1.0))
if OP==1:
    l = np.exp(-alp * xauxf ** pf)
    L = np.diag(l)
elif OP==2:
    l = (np.pi*xauxf)**(-1) * np.sin(np.pi*xauxf)
    l[0] = 1.0
    L = np.diag(l)
elif OP==3:
    l = (1+np.cos(np.pi*xauxf))/2
    L = np.diag(l)
elif OP==4:
    l = (1+np.cos(np.pi*xauxf))/2
    l = l**4*(35 - 84*l + 70*l**2 - 20*l**3)
    L = np.diag(l)
elif OP==5:
    #l = 0.5*scs.erfc(2*pf**0.5*(np.abs(xauxf) - 0.5)*np.sqrt((-np.log(1-4*(xauxf-0.5)**2))/(4*(xauxf-0.5)**2)))
    l = 0.5*scs.erfc(2*pf**0.5*(np.abs(xauxf)-0.5))
    L = np.diag(l)
elif OP==6:
    l = 0.5*scs.erfc(2*pf**0.5*(np.abs(xauxf) - 0.5)*np.sqrt((-np.log(1-4*(xauxf-0.5)**2))/(4*(xauxf-0.5)**2)))
    if (np.mod(i,2)==0):
        l[int(i/2)]=0.5
    L = np.diag(l)

#what if L is the identity
# L = np.identity(j)
# L[i,i] = 0
# L[i-1,i-1] = 0

#M = inv[B.T@W@B]@B.T@W, lets proof the properties of B@W in all the possibilities
WP1 = W1@P1
WP2 = W1@P2
RWP1 = U1
RWP2 = U2

#Matrix M = C@P.T@W
MP1 = C1@P1.T@W1
MP2 = C2@P2.T@W1

#Filter matrix P@L@M
F1 = P1@L@C1@P1.T@W1
F4 = U2@L@U2.T

#function

kas = 15
times = 20
norm1 = np.zeros((kas,times+1))
norm1[:,0] = 1
norm2 = np.zeros((kas,times+1))
norm2[:,0] = 1

plt.ion()
plt.figure(1, figsize=(27, 15))
style.use('ggplot')

colors = np.arange(0.05,0.85,step=(0.8/kas))

for k in range(0,kas):
    # ko = 20
    yi = np.sin(2*(k+1)*xi*np.pi)
    yo1 = yi
    yo2 = yi
       
    for j in range(0,times):
        yf1 = F1@yo1
        yf2 = F4@yo2
        yn1 = ln.norm(yf1)/ln.norm(yi)
        yn2 = ln.norm(yf2)/ln.norm(yi)
        yo1 = yf1
        yo2 = yf2
        norm1[k,j+1] = yn1
        norm2[k,j+1] = yn2
    
    plt.subplot(1, 2, 1)
    plt.plot(range(0,times+1), norm1[k,:],c=str(colors[k]))
    plt.title('Non-Unitary Filter', fontsize=14)
    plt.ylim([0.6, 1.4])
    plt.xlim([-1, times])
    plt.xlabel("Number of Applications", fontsize=12)
    plt.ylabel("Norm of the vecotor U", fontsize=12)
    
    plt.subplot(1, 2, 2)
    plt.plot(range(0,times+1), norm2[k,:],c=str(colors[k]))
    plt.title('Unitary Filter', fontsize=14)
    plt.ylim([0.6, 1.4])
    plt.xlim([-1, times])
    plt.xlabel("Number of Applications", fontsize=12)
    plt.ylabel("Norm of the vecotor U", fontsize=12)
    
    plt.draw()
    titulo = '# OF FILTER SUCCESSIVE APPLICATIONS vs NORM OF VECTOR \n u0 = sin(2*pi*k*x)'
    plt.suptitle(titulo, fontsize=18)
    plt.pause(0.01)



    
# times = 5
# yi = np.sin(xi*np.pi) + (1/4)*np.sin(4*xi*np.pi) + (1/8)*np.sin(8*xi*np.pi) #+ (1/16)*np.sin(30*xi*np.pi)
# yi2 = np.sin(xi*np.pi) + (1/4)*np.sin(4*xi*np.pi) + (1/8)*np.sin(8*xi*np.pi)
# #yi = np.sin(60*xi*np.pi)
# yo1 = yi
# yo2 = yi
# norm1 = np.zeros((1,times+1))
# norm1[:,0] = 1
# norm2 = np.zeros((1,times+1))
# norm2[:,0] = 1

# plt.ion()
# plt.figure(1, figsize=(27, 15))
# style.use('ggplot')

# for j in range(0,times):
    
#     plt.clf()
    
#     yf1 = F1@yo1
#     yf2 = F4@yo2
#     yn1 = ln.norm(yf1)/ln.norm(yi)
#     yn2 = ln.norm(yf2)/ln.norm(yi)
#     yo1 = yf1
#     yo2 = yf2
#     norm1[0,j+1] = yn1
#     norm2[0,j+1] = yn2

#     plt.subplot(1, 2, 1)
# #    plt.subplot(1, 1, 1)
# #    plt.semilogy(xi, np.abs(yf1-yi),'r')
# #    plt.semilogy(xi, np.abs(yf2-yi),'b')
#     plt.plot(xi, yf1,'r')
#     plt.plot(xi, yi,'g')
#     plt.title('Non-Unitary Filter')
# #    plt.ylim([-1.5, 1.5])
#     plt.xlim([-1, 1])
    
#     plt.subplot(1, 2, 2)
# #    plt.semilogy(xi, np.abs(yf2-yi2),'b')
#     plt.plot(xi, yf2,'b')
#     plt.plot(xi, yi,'g')
#     plt.title('Unitary Filter')
# #    plt.title('Original')
# #    plt.ylim([-1.5, 1.5])
#     plt.xlim([-1, 1])
    
#     plt.draw()
#     titulo = 'Filters'
#     plt.suptitle(titulo)
#     plt.pause(1)

# =============================================================================
# PROOF ABOUT FILTER ORDER
# =============================================================================


# i = 25
# j = i+1
# alpha = 0
# beta = 0

# #vandermonde matrix
# xi = fs.jacobigl(alpha,beta,i)
# P1 = fs.jacobip(xi,alpha,beta,i).T
# P2 = fs.legendre(xi,i).T

# #gll weights
# w1 = 2/(j*(j-1)*(P2[:,i]**2))
# W1 = np.diag(w1)
# RW = np.diag(np.sqrt(w1))

# #c for inv(P.T@W@P)
# c1 = np.ones(i+1)
# c1[i] = i/(2*i+1)
# C1 = np.diag(c1)

# c2 =  np.arange(0,i+1) + 1/2
# c2[i] = i/2
# C2 = np.diag(c2)

# #Sumedh alternative
# U1 = RW@P1
# U2 = RW@P1@(np.diag(np.sqrt(c1)))


# #M = inv[B.T@W@B]@B.T@W, lets proof the properties of B@W in all the possibilities
# WP1 = W1@P1
# WP2 = W1@P2
# RWP1 = U1
# RWP2 = U2

# #Matrix M = C@P.T@W
# MP1 = C1@P1.T@W1
# MP2 = C2@P2.T@W1


# #filter function matrix

# alp = -np.log(1.0e-16)       # Machine precision according to Diamesis 2008
# pfs = np.arange(1, 25, step=0.25)

# sU = []
# sN = []

# ####
# OP = 5      #only OP 1,5 and 6 are supported
# ####

# alp = -np.log(1.0e-16)       # Machine precision according to Diamesis 2008

# j = i+1    
# xauxf = np.arange(1,j+1)
# xauxf = ((xauxf - 1.0) / (j - 1.0))

# for pf in pfs:
#     #pf = 10
   
#     #j = i+1    
#     #l = np.arange(1,j+1)
#     #l = np.exp(-alp * ((l - 1.0) / (j - 1.0)) ** pf)
#     #L = np.diag(l)
#     #what if L is the identity
#     #L = np.identity(j)
    
#     if OP==1:
#         l = np.exp(-alp * xauxf ** pf)
#         L = np.diag(l)
#     elif OP==5:
#         #l = 0.5*scs.erfc(2*pf**0.5*(np.abs(xauxf) - 0.5)*np.sqrt((-np.log(1-4*(xauxf-0.5)**2))/(4*(xauxf-0.5)**2)))
#         l = 0.5*scs.erfc(2*pf**0.5*(np.abs(xauxf)-0.5))
#         L = np.diag(l)
#     elif OP==6:
#         l = 0.5*scs.erfc(2*pf**0.5*(np.abs(xauxf) - 0.5)*np.sqrt((-np.log(1-4*(xauxf-0.5)**2))/(4*(xauxf-0.5)**2)))
#         if (np.mod(i,2)==0):
#             l[int(i/2)]=0.5
#         L = np.diag(l)
    
#     #Filter matrix P@L@M
#     F1 = P1@L@C1@P1.T@W1
#     F2 = P2@L@C2@P2.T@W1
#     #by this way, F1=F2
#     #SVD of F1
#     Q1,S1,Z1 = ln.svd(F1)
#     sN.append(S1[0])
#     #rank of F1
#     #r1 = ln.matrix_rank(F1)
#     #condition number of F1, condition number has no sense since F is singular
#     #con1 = ln.cond(F1)
    
#     #filter matrix from Sumedh U@L@inv(U), in the unitary case inv(U)=U.T
#     #F3 = U1@L@ln.inv(U1)
#     F4 = U2@L@U2.T
#     #F5 = U1@L@U1.T, not even compute it, if L is de identity, F should be the identity also
#     #that is not the case
#     #by this way, F3=F4
#     #SVD of F4
#     Q2,S2,Z2 = ln.svd(F4)
#     sU.append(S2[0])
#     #rank of F4
#     #r2 = ln.matrix_rank(F4)
#     #condition number of F4, condition number has no sense since F is singular
#     #con2 = ln.cond(F4)

# sU = np.array(sU)
# sN = np.array(sN)

# plt.plot(pfs, sU,'b', label=("Unitary"))
# plt.plot(pfs, sN,'r', label=("Non-Unitary"))
# plt.title('Highest Singular Value')
# plt.ylim([0, 2])
# plt.xlim([0, 25])
# plt.xlabel("Filter Order P")
# plt.ylabel("\u03C31")
# plt.legend(loc=4)

# for i in range(0,26):
#     L[i,i] = 1/(i+1)






























