#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Sergio Aldair Castiblanco Ballesteros
#Universidad Nacional de Colombia
#Ingeniería Civil

#This code computes the first derivative of a function c(x), using the Fourier
#pseudo-spectral method (that is, at collocation points). That function c(x)
#is user defined. The dealiasing process is performed
#
#list of variables:
#    c = userdefined function c(x)
#    dc = analytical solution of the first derivative of c(x)
#    N = number of collocation points in the domain
#    L = length of the domain
#    x = Fourier nodes
#    gap = distance between Fourier nodes
#    w = exp(2*pi*i/N)   ;   where i=sqrt(-1)
#    F = basic Fourier approximation matrix, for 0 <= k <= N-1
#        (nodal to modal):   F*Cm = Cn
#    Fs = Fourier approximation matrix from Steven G. JOhnson notes
#        (nodal to modal):   Fs*Cm = Cn (same as the Peyret's book)
#    DFs = Fourier firs derivative matrix from Steven G. Johnson notes
#        (nodal to modal derivative):    DFs*Cm = dCn
#    DFp = Fourier first derivative matrix from Peyret's book
#        (nodal to nodal derivative):    DFp*Cn = dCn
#    Cn = Nodal values of the function c(x) (Physical space)
#    Cm = Modal values of the function c(x) (Spectral space)
#    dCn = Nodal values of the first derivative of the function c(x)
#    dCm = Modal values of the first derivative of the function c(x)    
#    error = error calculated between analytical dc and numerical dCn
#    *Please note that de Steven G. Johnson approximation is the same of the
#       Pyret's book. In the first one the derivative is computed from the
#       modal values of C, in the second case the derivative is computed
#       directly from the nodal values of C. What's better, or maybe useful,
#       for you?
#
#For the purpous exposed above, the metodology is based on professor Strang
#MIT Fourier Courses, on the notes of FFT-based differentiation, professor 
#Steven G. Johnson, MIT, and on the Spectral Methods for Incompressible
#Viscous Flow book of Roger Peyret.

# =============================================================================
# Import libraries and setting the folder direction
# =============================================================================

import numpy as np
from numpy import linalg
import math
import os
import matplotlib.pyplot as plt
from matplotlib import style
#os.chdir('/home/debian3/Documents/Sergio/Proyecto/Codigos/Fourier_py')

# =============================================================================
# Entry values
# =============================================================================
N = 16              #Number of points in the domain
                    #Number of points in the domain of fourier N-1
L = 2*math.pi
gap = L/N
x = np.linspace(0,L-gap,N)

# =============================================================================
# Evaluating the user defined function in the x array
# =============================================================================
#c = x*(x-math.pi)*(x-2*math.pi)
#dc = 3*x**2 - 6*math.pi*x + 2*math.pi**2

#c = np.sin(x)
#dc = np.cos(x)

c = 4*np.exp(-8*(x-np.pi)**2)
dc = -8*2*x*4*np.exp(-8*(x-np.pi)**2)

#c = (30/np.sqrt(4*math.pi))*np.exp(-((x-(L/2))**2)/(4))
#dc = (30/np.sqrt(4*math.pi))*(-(x-(L/2))/2)*np.exp(-((x-(L/2))**2)/(4))

#c = np.tanh(x-(L/2))
#dc = (1/(np.cosh(x-(L/2)))**2)

#c = np.tanh(4*(x-math.pi))*x*(x-2*math.pi)*(1/10)
#dc = (1/10)*(2*(x-math.pi)*np.tanh(4*(x-math.pi)) + 4*x*(x-2*math.pi)*(1/(np.cosh(4*math.pi-4*x))**2))

#c = np.sin(8*x)
#dc = 8*np.cos(8*x)

#c = (x-math.pi)**3
#dc = 3*(x-math.pi)**2

# =============================================================================
# Fourier approximation Matrix NxN from (0,N-1) (Strang)
# =============================================================================
F = np.zeros((N,N),dtype=np.complex)
w = np.exp(complex(0,2*math.pi/N))
for i in range(0,N):
    for j in range(0,N):
        F[i,j]=w**(i*j)

# =============================================================================
# Fourier waveless approximation matrix NxN from (-N/2 to N/2)
#   (Jhonson and Peyret)
# =============================================================================
Fs = np.ones((N,N),dtype=np.complex)
if N%2==0:
    Fs[:,int(N/2)]=np.cos(math.pi*N*x/L) 
for k in range(1,int(np.ceil(N/2))):
    for n in range(0,N):
        Fs[n,k] = w**(n*k)
        Fs[n,N-k] = w**(-n*k)

# =============================================================================
# Fourier waveless approximation first derivative matrix NxN from (-N/2 to N/2)
#   (Jhonson)        
# =============================================================================
DFs = np.zeros((N,N),dtype=np.complex)
for k in range(1,int(np.ceil(N/2))):
    for n in range(0,N):
        DFs[n,k] = (2*math.pi/L)*1j*k*w**(n*k)
        DFs[n,N-k] = -(2*math.pi/L)*1j*k*w**(-n*k)

# =============================================================================
# Fourier waveless approximation first derivative matrix NxN from (-N/2 to N/2)
#   (Peyret) 
# =============================================================================
DFp = np.zeros((N,N))
if N%2!=0:
    for k in range(0,N):
        for n in range(0,N):
            if n!=k:
                DFp[k,n]=((-1)**(k+n))/(2*np.sin((x[k]-x[n])/2))
                
else:
    for k in range(0,N):
        for n in range(0,N):
            if n!=k:
                DFp[k,n]=0.5*((-1)**(k+n))/(np.tan((x[k]-x[n])/2))

# =============================================================================
# Computing the solution of the function and firs derivative approximations
# =============================================================================
Cm = linalg.solve(F,c)          #Direct Fourier approximation (Strang)
Cms = linalg.solve(Fs,c)        #Waveless Fourier approximation (Jhonson)
Cn = np.dot(F,Cm)               #Direct Fourier approximation (Strang)
Cns = np.dot(Fs,Cms)            #Waveless Fourier approximation (Jhonson)
dCn = np.dot(DFs,Cm)           #Waveless Fourier der. approximation (Jhonson)
error = np.abs(dCn-dc)          #Error of the derivative approximation
dCp = np.dot(DFp,c)             #Pyret's direct way for the derivative
errorp = np.abs(dCp-dc)         #Error of Pyret's derivative method
norm = linalg.norm(error[1:N-1])










