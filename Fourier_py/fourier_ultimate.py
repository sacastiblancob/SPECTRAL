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

import os
os.chdir("F:\SPECTRAL\Codes\Fourier_py")

import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from matplotlib import style
#os.chdir('/home/debian3/Documents/Sergio/Proyecto/Codigos/Fourier_py')

# =============================================================================
# Initial values
# =============================================================================
N = 16             #Number of points in the domain
                    #Number of points in the domain of fourier N-1
L = 2*np.pi
gap = L/N
x = np.linspace(-np.pi,np.pi-gap,N)

# =============================================================================
# Evaluating the user defined function in the x array
# =============================================================================
#c = x*(x-math.pi)*(x-2*math.pi)
#dc = 3*x**2 - 6*math.pi*x + 2*math.pi**2

#c = np.sin(x)
#dc = np.cos(x)

c = (30/np.sqrt(4*np.pi))*np.exp(-(x**2)/4)
dc = (30/np.sqrt(4*np.pi))*(-(x-(L/2))/2)*np.exp(-(x**2)/4)

#c = np.tanh(x)
#dc = (1/(np.cosh(x))**2)

#c = np.sin(8*x)
#dc = 8*np.cos(8*x)

#c = x**3
#dc = 3*x**2

# =============================================================================
# Fourier approximation Matrix NxN from (-N/2,N/2)
# =============================================================================
F = np.zeros((N,N),dtype=np.complex)
if N%2==0:    
    k = np.array(range(int(-N/2),int((N/2))))
else:
    k = np.array(range(int(-N/2),int((N/2)+1)))

for i in k:
    for j in range(0,N):
        F[i,j]=np.exp(1j*i*x[j])
        
# =============================================================================
# Transform
# =============================================================================
cm = linalg.solve(F,c)
abscm = np.abs(cm)
cn = np.dot(F,cm)
error = np.abs(c-cn)
meannorm = np.sqrt(np.dot(error,error))
infnorm = np.max(error)

# =============================================================================
# Postprocessing
# =============================================================================
a=N*10
xaux = np.linspace(-np.pi,np.pi,a)
yaux = (30/np.sqrt(4*np.pi))*np.exp(-(xaux**2)/4)
freq = np.zeros((N,a),dtype=np.complex)

for i in range(0,N):
    freq[i,:]=cm[i]*np.exp(xaux*1j*k[i])

fcos = freq.real
fsin = freq.imag
aprox = np.zeros(a)
for i in range(0,a):
    aprox[i] = -(np.sum(fcos[:,i])+np.sum(fsin[:,i]))+2*np.mean(c)
    
## =============================================================================
## Plotting
## =============================================================================
#
plt.ion()
style.use('ggplot')
plt.figure(1, figsize=(11.5, 11.5))

#Mass plot

fig1 = plt.subplot(1,1,1)
plot1 = fig1.bar(k,abscm)
fig1.set_xlabel(r'K')
fig1.set_ylabel(r'C modal')
fig1.set_title('Energy Spectrum of F (-N/2 to N/2)')

plt.ion()
style.use('ggplot')
plt.figure(1, figsize=(11.5, 11.5))

fig1 = plt.subplot(1, 2, 1)
plot1 = fig1.plot(xaux,yaux)
plot2 = fig1.plot(xaux,aprox)
plot3 = fig1.plot(x,c,'o')
fig1.set_xlim(x[0], x[N-1]+gap)
fig1.set_ylim([np.floor(np.min(c))-1,np.ceil(np.max(c))+1])
fig1.set_xlabel(r'x')
fig1.set_ylabel(r'C')
#    fig1.tick_params(axis='both', which='major', labelsize=6)
fig1.set_title('Analytical solution and F aprox')

fig2 = plt.subplot(1, 2, 2)
fig2.set_xlim(x[0], x[N-1]+gap)
fig2.set_ylim([-np.ceil(np.max(c))-1,np.ceil(np.max(c))+1])
fig2.set_xlabel(r'x')
fig2.set_ylabel(r'C')
fig2.set_title('Fourier Decomposition')


for i in range(0,N):
  
   plt.clf()
   plot4 = fig2.plot(xaux,fcos[i,:])
   plot5 = fig2.plot(xaux,fsin[i,:])
   #    fig2.tick_params(axis='both', which='major', labelsize=6)

   plt.draw()
   plt.pause(0.1)