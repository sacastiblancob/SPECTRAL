#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Sergio Aldair Castiblanco Ballesteros
#Universidad Nacional de Colombia
#Ingeniería Civil

#This code solve the advection equation 1D
#               dC/dt + u dc/dx = f(x)
#
# =============================================================================
# Import libraries and setting the folder direction
# =============================================================================

import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from matplotlib import style
from functions import FTrans
from functions import IFTrans

# =============================================================================
# Entry and grid values
# =============================================================================

#space grid
N = 128         #Only it's functional with evens
L = 2*np.pi
gap = L/N
x = np.linspace(-np.pi,np.pi-gap,N)

#Wave frequencies K
if N%2==0:
    K = np.linspace(-int(N/2),int(N/2)-1,N)
else:
    K = np.linspace(-int(N/2),int(N/2),N)

#velocity
u = 3

#Time grid
to = 0
tf = 2
t = np.linspace(to,tf,50)
Ts = np.size(t)
dt = t[1] - t[0]

# =============================================================================
# Solution
# =============================================================================

#Storage vectors
C = np.zeros((N,Ts))
A = np.zeros((N,Ts))
error = np.zeros(Ts)

#Initial conditions
C[:,0] = np.cos(2*x)
A[:,0] = np.cos(2*x)
error[0] = linalg.norm(C[:,0]-A[:,0])
#C[:,0] = x*(x-np.pi)*(x+np.pi)/10
#C[:,0] = np.exp(-16*x**2)
#C[:,0] = np.tanh(x)
Com = FTrans(C[:,0],N)

#Integrator factor
gamma = -1j*K*u
I = np.exp(1*gamma*dt)

#Time loop
for time in range(1,Ts):
    Cm = I*Com
    C[:,time] = IFTrans(Cm,N)
    A[:,time] = A[:,0] = np.cos(2*(x - u*t[time]))
    error[time] = linalg.norm(C[:,time]-A[:,time]) 
    Com = Cm

# =============================================================================
# Plot
# =============================================================================
plt.ion()
style.use('ggplot')
plt.figure(1, figsize=(11.5, 11.5))

fig1 = plt.subplot(1, 1, 1)
fig1.set_xlim(x[0], L)
fig1.set_ylim([np.floor(np.min(C[:,0]))-1,np.ceil(np.max(C[:,0]))+1])
fig1.set_xlabel(r'x')
fig1.set_ylabel(r'C')
fig1.set_title('Advection')

for i in range(0,Ts):
    
    plt.clf()
    fig1 = plt.subplot(1, 1, 1)
    fig1.set_xlim(x[0], x[N-1])
    fig1.set_ylim([np.floor(np.min(C[:,0]))-1,np.ceil(np.max(C[:,0]))+1])
    fig1.set_xlabel(r'x')
    fig1.set_ylabel(r'C')
    fig1.set_title('Advection')
    plot1 = fig1.plot(x,C[:,i])
    plot2 = fig1.plot(x,A[:,i])
    #    fig2.tick_params(axis='both', which='major', labelsize=6)

    plt.draw()
    plt.pause(0.1)












