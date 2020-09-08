# -*- coding: utf-8 -*-
"""
SERGIO CASTIBLANCO
31-03-2020

THIS CODE IS MADE FOR PLOT THE JACOBI POLYNOMIALS

"""
import os
os.chdir("F:\SPECTRAL\Codes\Jacobi_py")

import numpy as np
import funspe as fs
import matplotlib.pyplot as plt
from matplotlib import style
#from numpy import linalg as ln


#COMPUTING THE JACOBI BASIS FUNCTIONS

N = 5

#if alpha=beta=0.0 we are talking about Legendre Polynomials
#if alpha=beta=1.0 we are talking about Chebyshev Polynomials
alpha = 0.0
beta = 0.0
x = np.arange(-1,1.0005,0.001)
Por = fs.legendre(x,N)
Pon = fs.jacobip(x,alpha,beta,N)
zeros, w = fs.jacobigq(alpha,beta,N-1)

#PLOT
#do you want the just the orthogonal or the orthonormal version?
# 1 for orthogonal
# 2 for orthonormal
ortho = 2
Pplot = Pon
if (ortho == 1):
    Pplot = Por

plt.ion()
plt.figure(1, figsize=(11, 8.5))
style.use('ggplot')

# plt.subplot(1, 1, 1)
# plt.plot(x, Pplot[0,:])
# plt.title('Legendre Polynomials')
# plt.xlabel(r'X')
# plt.ylabel(r'L(x)')
# plt.draw()

for i in range(0,N+1):
    
    #plt.clf()
    
    plt.subplot(1, 1, 1)
    #plt.plot(x, Pplot[i,:], 'b')
    plt.plot(x, Pplot[i,:],label=("L"+str(i)+"(x)"))
    plt.xlim([-1.05, 1.05])
    plt.ylim([ np.min(Pplot)-0.1,np.max(Pplot)+0.1])
    plt.ylabel(r'L(X)')
    plt.title('Legendre Polynomials')
    
        
    plt.draw()
    plt.pause(0.01)

plt.subplot(1, 1, 1)
plt.plot(zeros,np.zeros(N),'bs')
plt.draw()


