#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Sergio Aldair Castiblanco Ballesteros
#Universidad Nacional de Colombia
#Ingeniería Civil

#This code compute the Fourier coefficients for an array with the values
#of any function defined periodic between 0 and 2pi
#Based on Strang MIT Fourier Courses 
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
N = 16                      #Number of points in the domain
L = 2*math.pi
gap = L/N
x = np.linspace(0,L-gap,N)
K = np.array(range(0,N))

# =============================================================================
# Evaluating the user defined function in the x array
# =============================================================================
#y = 3/(5-4*np.cos(x))
#y = np.tanh(x-(L/2))
y = (30/np.sqrt(4*math.pi))*np.exp(-((x-(L/2))**2)/(4))
#y = np.sin(x-(L/2))

# =============================================================================
# Computing the Fourier Matrix NXN
# =============================================================================
F = np.zeros((N,N),dtype=np.complex)
w = np.exp(complex(0,2*math.pi/N))
for i in range(0,N):
    for j in range(0,N):
        F[i,j]=w**(i*j)

# =============================================================================
# Another point of view (FFT-based differentiation of Steven G. Johnson, MIT)
# =============================================================================
FS = np.ones((N,N),dtype=np.complex)
if N%2==0:
    FS[:,int(N/2)]=np.cos(math.pi*N*x/L) 
for k in range(1,int(np.ceil(N/2))):
    for n in range(0,N):
        FS[n,k] = w**(n*k)
        FS[n,N-k] = w**(-n*k)

# =============================================================================
# Computing the complex Fourier Frequencies of the solution C
# =============================================================================
f = linalg.solve(F,y)
absf = np.abs(f)
fs = linalg.solve(FS,y)

# =============================================================================
# Solution with the Fourier Equation F=sum(Ci*)
# =============================================================================
yn = np.dot(F,f)
yns = np.dot(FS,fs)
yn2 = yn.real
yns2 = yns.real
error = np.abs(y-yn2)
errors = np.abs(y-yns2)
dec = np.zeros((N,N),dtype=np.complex)
for i in range(0,N):
    dec[i,:]=f*F[i,:]
decr = dec.real
deci = dec.imag

a=N*20
xaux = np.linspace(0,2*math.pi,a)
yaux = (30/np.sqrt(4*math.pi))*np.exp(-((xaux-(L/2))**2)/(4))
freq = np.zeros((N,a),dtype=np.complex)
freqs = np.zeros((N,a),dtype=np.complex)

for i in range(0,N):
    freq[i,:]=f[i]*np.exp(xaux*1j*i)

freqs[0,:] = fs[0]
if N%2==0:
    freqs[int(N/2),:]=fs[int(N/2)]*np.cos(math.pi*N*xaux/L)
for i in range(1,int(np.ceil(N/2))):
    freqs[i,:]=fs[i]*np.exp(2*math.pi*1j*i*xaux/L)
    freqs[N-i,:]=fs[N-i]*np.exp(-2*math.pi*1j*i*xaux/L)
    
fcos = freq.real
fcoss = freqs.real
fsin = freq.imag
fsins = freqs.imag
aprox = np.zeros(a)
aproxs = np.zeros(a)
for i in range(0,a):
    aprox[i] = np.sum(fcos[:,i])+np.sum(fsin[:,i])
    aproxs[i] = np.sum(fcoss[:,i])+np.sum(fsins[:,i])

# =============================================================================
# Plot
# =============================================================================
plt.ion()
style.use('ggplot')
plt.figure(1, figsize=(11.5, 11.5))

#Mass plot
fig1 = plt.subplot(1,1,1)
plot1 = fig1.bar(K,absf)
fig1.set_xlabel(r'K')
fig1.set_ylabel(r'C modal')
fig1.set_title('Energy Spectrum of F (Strang 0-N)')
#
fig1 = plt.subplot(1, 2, 1)
plot1 = fig1.plot(xaux,yaux)
plot4 = fig1.plot(xaux,aprox)
plot6 = fig1.plot(xaux,aproxs)
plot5 = fig1.plot(x,y,'o')
fig1.set_xlim(x[0], L)
fig1.set_ylim([np.floor(np.min(y))-1,np.ceil(np.max(y))+1])
fig1.set_xlabel(r'x')
fig1.set_ylabel(r'C')
#    fig1.tick_params(axis='both', which='major', labelsize=6)
fig1.set_title('Analytical solution and F aprox')

fig2 = plt.subplot(1, 2, 2)
fig2.set_xlim(x[0], L)
fig2.set_ylim([np.floor(np.min((decr,deci)))-1,np.ceil(np.max((decr,deci)))+1])
fig2.set_xlabel(r'x')
fig2.set_ylabel(r'C')
fig2.set_title('Fourier Decomposition')

for i in range(0,N):
    
#    plt.clf()
    plot2 = fig2.plot(xaux,fcos[i,:])
    plot3 = fig2.plot(xaux,fsin[i,:])
    #    fig2.tick_params(axis='both', which='major', labelsize=6)

    plt.draw()
    plt.pause(0.1)


