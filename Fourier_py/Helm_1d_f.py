#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Sergio Aldair Castiblanco Ballesteros
#Universidad Nacional de Colombia
#Ingeniería Civil

#This code compute the solution of the Helmholtz equation with the
#Fourier approximation  

#The helmholtz equation is:
#    eps*d2C/dx2 - C = eps*f(x)
#    where: f(x) --> user defined function
    

# =============================================================================
# Import libraries and setting the folder direction
# =============================================================================

import numpy as np
from numpy import linalg
import math
import os
import matplotlib.pyplot as plt
from matplotlib import style
os.chdir('/home/debian3/Documents/Sergio/Proyecto/Codigos/Fourier_py')

# =============================================================================
# Entry values
# =============================================================================
eps = 0.5
N1 = 5001             #Number of points in the domain
N = N1-1            #Number of points in the domain of fourier N-1
L = 2*math.pi
gap = L/N
x = np.linspace(0,L-gap,N)
D = 0.02            #Diffusion coefficient
t = 3               #time (to in steady condition)

# =============================================================================
# Evaluating the user defined function in the x array
# =============================================================================
y = (10/np.sqrt(4*math.pi*eps))*np.exp(-((x-math.pi)**2)/(4*eps))

# =============================================================================
# Computing the Fourier Matrix NXN
# =============================================================================
F = np.zeros((N,N),dtype=np.complex)
w = np.exp(complex(0,2*math.pi/N))
e = np.zeros((N,N),dtype=np.complex)
for i in range(0,N):
    for j in range(0,N):
        F[i,j]=-w**(i*j)*(eps*(j**2)+1)

# =============================================================================
# Computing the complex Fourier Frequencies of the solution C
# =============================================================================
f = linalg.solve(F,eps*y)

# =============================================================================
# Solution with the Fourier Equation F=sum(Ci*)
# =============================================================================
sol = np.dot(F,f)
sol2 = sol.real
dec = np.zeros((N,N),dtype=np.complex)
for i in range(0,N):
    dec[i,:]=f*F[i,:]
decr = dec.real
deci = dec.imag

# =============================================================================
# Plot
# =============================================================================
plt.ion()
style.use('ggplot')
plt.figure(1, figsize=(11.5, 11.5))

fig1 = plt.subplot(1, 2, 1)
plot1 = fig1.plot(x,y)
fig1.set_xlim(x[0], L)
fig1.set_ylim([np.floor(np.min(y))-1,np.ceil(np.max(y))+1])
fig1.set_xlabel(r'x')
fig1.set_ylabel(r'C')
#    fig1.tick_params(axis='both', which='major', labelsize=6)
fig1.set_title('Analytical solution')

fig2 = plt.subplot(1, 2, 2)
fig2.set_xlim(x[0], L)
fig2.set_ylim([np.floor(np.min((decr,deci)))-1,np.ceil(np.max((decr,deci)))+1])
fig2.set_xlabel(r'x')
fig2.set_ylabel(r'C')
fig2.set_title('Fourier Decomposition')

for i in range(0,N):
    
#    plt.clf()
    plot2 = fig2.plot(x,decr[:,i])
    plot3 = fig2.plot(x,deci[:,i])
    #    fig2.tick_params(axis='both', which='major', labelsize=6)

    plt.draw()
    plt.pause(0.1)









