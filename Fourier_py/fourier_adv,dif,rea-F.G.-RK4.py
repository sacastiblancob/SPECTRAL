#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Sergio Aldair Castiblanco Ballesteros
#Universidad Nacional de Colombia
#Ingeniería Civil

#This code solve the advection-diffusion-reaction equation 1D
#               dC/dt + u dc/dx = D d2c/dx2 + r*c
#Time discretization: First order (Euler)
#
# =============================================================================
# Import libraries and setting the folder direction
# =============================================================================

import os
os.chdir("F:\SPECTRAL\Codes\Fourier_py")

import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from matplotlib import style
from functions import FT as FTrans
from functions import IFT as IFTrans
from time import time as tm

# =============================================================================
# Entry and grid values
# =============================================================================

#space grid
N = 150         #Only it's functional with evens
L = 2*np.pi
gap = L/N
x = np.linspace(-np.pi,np.pi-gap,N)

#Wave frequencies K
if N%2==0:
    K = np.linspace(-int(N/2),int(N/2)-1,N)
else:
    K = np.linspace(-int(N/2),int(N/2),N)

#constants
u = 3
D = 0.05
r = 0

#Time grid
to = 2
tf = 2.5
t = np.linspace(to,tf,50)
Ts = np.size(t)
dt = t[1] - t[0]

# =============================================================================
# Solution
# =============================================================================

#Storage vectors
C = np.zeros((N,Ts))
A = np.zeros((N,Ts))
error = np.zeros((N,Ts))

#Initial conditions
#a = np.pi/2
a = 0
C[:,0] = 1/np.sqrt(4*np.pi*D*(to))*np.exp(-((x+a)**2)/(4*D*(to)))
A[:,0] = 1/np.sqrt(4*np.pi*D*(to))*np.exp(-((x+a)**2)/(4*D*(to)))
            #Analytical of advection-diffusion equation
#C[:,0] = 4*np.exp(-8*(x)**2)
#A[:,0] = 4*np.exp(-8*(x)**2)
#C[:,0] = np.cos(2*(x-to*u))
#A[:,0] = np.cos(2*(x-to*u))
error[:,0] = np.abs(C[:,0]-A[:,0])
#C[:,0] = np.cos(2*x)
#C[:,0] = x*(x-np.pi)*(x+np.pi)/10
#C[:,0] = 4*np.exp(-8*(x+np.pi/2)**2)
#C[:,0] = np.tanh(x)

Com = FTrans(C[:,0],N)

#Integrator factor
gamma = -1j*K*u
I = np.exp(1*gamma*dt)

#Time loop
tini = tm()
for time in range(1,Ts):
    Cm1 = -(D*K**2 + r)*Com
    Cm2 = -(D*K**2 + r)*(Com+Cm1*(dt/2))
    Cm3 = -(D*K**2 + r)*(Com+Cm2*(dt/2))
    Cm4 = -(D*K**2 + r)*(Com+Cm3*dt)
    Cm = Com*I + np.exp(1*gamma*dt)*(dt/6)*(Cm1 + 2*Cm2 + 2*Cm3 + Cm4)
    C[:,time] = IFTrans(Cm,N)
#    A[:,time] = 4*np.exp(-8*(x-(t[time]-to)*u)**2)
#    A[:,time] = np.cos(2*(x-(t[time])*u))
    A[:,time] = 1/np.sqrt(4*np.pi*D*(t[time]))*np.exp(-((x+a-u*(t[time]-to))**2)/(4*D*(t[time])))
            #Analytical of advection-diffusion equation
    error[:,time] = np.abs(C[:,time]-A[:,time])
    Com = Cm
tfin = tm()
tt = tfin - tini
# =============================================================================
# Plotting
# =============================================================================

plt.ion()
style.use('ggplot')
plt.figure(1, figsize=(11.5, 11.5))

fig1 = plt.subplot(1, 1, 1)
fig1.set_xlim(x[0], L)
fig1.set_ylim([np.floor(np.min(C[:,0]))-1,np.ceil(np.max(C[:,0]))+1])
fig1.set_xlabel(r'x')
fig1.set_ylabel(r'C')
fig1.set_title('Advection-Diffusion-Reaction , Fourier-Galerkin')

for i in range(0,Ts):
    
    plt.clf()
    fig1 = plt.subplot(1, 1, 1)
    fig1.set_xlim(x[0], x[N-1])
    fig1.set_ylim([np.floor(np.min(C[:,0]))-1,np.ceil(np.max(C[:,0]))+1])
    fig1.set_xlabel(r'x')
    fig1.set_ylabel(r'C')
    fig1.set_title('Advection-Diffusion-Reaction , Fourier-Galerkin')
    plot1 = fig1.plot(x,C[:,i])
    plot2 = fig1.plot(x,A[:,i])
    #    fig2.tick_params(axis='both', which='major', labelsize=6)

    plt.draw()
    plt.pause(0.1)