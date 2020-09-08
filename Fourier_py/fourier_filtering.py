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

import numpy as np
#from numpy import linalg
import matplotlib.pyplot as plt
from matplotlib import style
from functions import FT
from functions import IFT
from scipy import special

# =============================================================================
# Entry and grid values
# =============================================================================

#space grid
N = 100         #Only it's functional with evens
L = 2*np.pi
gap = L/N
x = np.linspace(-np.pi,np.pi-gap,N)

#Wave frequencies K
if N%2==0:
    K = np.linspace(-int(N/2),int(N/2)-1,N)
else:
    K = np.linspace(-int(N/2),int(N/2),N)

#constants
u = 6
D = 0
r = 0

#filtering parameters
P = 5                               #Order of the filter
alfa = -np.log(np.spacing(1))       #Filter constant alfa, machine precission
gammafe = np.exp(-alfa*(K/N)**P)


Xbv = np.abs(x) - 0.5
gammabv = 0.5*(1 - special.erf(2*Xbv*((-N*np.log(1-4*Xbv**2)/(4*Xbv**2))**0.5)))

#Time grid
to = 2
tf = 4
t = np.linspace(to,tf,100)
Ts = np.size(t)
dt = t[1] - t[0]

# =============================================================================
# Solution
# =============================================================================

#Storage vectors
C = np.zeros((N,Ts))
Cf = np.zeros((N,Ts))
Cfbv = np.zeros((N,Ts))
A = np.zeros((N,Ts))
CM = np.zeros((N,Ts))
E = np.zeros(Ts,dtype=np.complex)
error = np.zeros((N,Ts))

#Initial conditions
#a = np.pi/2
a = 0
#C[:,0] = 10/np.sqrt(4*np.pi*D*(to))*np.exp(-((x+a)**2)/(4*D*(to)))
#A[:,0] = 10/np.sqrt(4*np.pi*D*(to))*np.exp(-((x+a)**2)/(4*D*(to)))
            #Analytical of advection-diffusion equation

C[10:20,0] = 1
#C[:,0] = 16*np.exp(-8*(x)**2)
#A[:,0] = 16*np.exp(-8*(x)**2)
#C[:,0] = np.cos(2*(x-to*u))
#A[:,0] = np.cos(2*(x-to*u))
#error[:,0] = np.abs(C[:,0]-A[:,0])
#C[:,0] = np.cos(2*x)
#C[:,0] = x*(x-np.pi)*(x+np.pi)/10
#C[:,0] = 4*np.exp(-8*(x+np.pi/2)**2)
#C[:,0] = np.tanh(x)
Com = FT(C[:,0],N)
E[0] = np.sum(np.abs(Com)**2)
CM[:,0] = np.abs(Com)

#Integrator factor
gamma = -1j*K*u - D*K**2 - r
I = np.exp(1*gamma*dt)

#Time loop
for time in range(1,Ts):
    Cm = I*Com
    Cmf = Cm*gammafe
    Cmbv = Cm*gammabv
    CM[:,time] = np.abs(Cm)
    E[time] = np.sum(np.abs(Cm)**2)
    C[:,time] = IFT(Cm,N)
    Cf[:,time] = IFT(Cmf,N)
    Cfbv[:,time] = IFT(Cmbv,N)
#    A[:,time] = 16*np.exp(-8*(x-(t[time]-to)*u)**2)
#    A[:,time] = np.cos(2*(x-(t[time])*u))
#    A[:,time] = 10/np.sqrt(4*np.pi*D*(t[time]))*np.exp(-((x+a-u*(t[time]-to))**2)/(4*D*(t[time])))
            #Analytical of advection-diffusion equation
#    error[:,time] = np.abs(C[:,time]-A[:,time])
    Com = Cm
#print(error[:,Ts-1])
#print(E)
# =============================================================================
# Plotting
# =============================================================================

#plt.ion()
#style.use('ggplot')
#plt.figure(1, figsize=(11.5, 11.5))
#
##Energy plot
#fig1 = plt.subplot(1,1,1)
#plot1 = fig1.plot(t,E)
#fig1.set_xlabel(r'Time')
#fig1.set_ylabel(r'Energy')
#fig1.set_title(r'Time Energy Plot')

##Mass plot
#
#for i in range(0,Ts)
#    plt.clf()
#    fig1 = plt.subplot(1,1,1)
#    plot1 = fig1.bar(K,abscm)
#    fig1.set_xlabel(r'K')
#    fig1.set_ylabel(r'C modal')
#    fig1.set_title('Energy Spectrum of F (-N/2 to N/2)')

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
#    plot2 = fig1.plot(x,A[:,i])
    plot3 = fig1.plot(x,Cf[:,i])
#    plot4 = fig1.plot(x,Cfbv[:,i])
    #    fig2.tick_params(axis='both', which='major', labelsize=6)

    plt.draw()
    plt.pause(0.1)

