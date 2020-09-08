#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Sergio Aldair Castiblanco Ballesteros
#Universidad Nacional de Colombia
#Ingeniería Civil

#This code solve the next differential equation:
#
#               dc/dx = dc
#where C is the solution, and dc is an user defined function and also, as you
#can see, the firs derivative of C
#
#And find the Modal values of C (Fourier Frequencies of the solution C), and
#the Nodal values of C ussing the modal results and the Fourier linear formula:
#               Cn = sum_k(Cm*exp(i*k*x));   Cn=C nodal; Cm = C modal
#Based on Strang MIT Fourier Courses and the notes of FFT-based differentiation
#of Steven G. Johnson, MIT
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
N = 11              #Number of points in the domain
                    #Number of points in the domain of fourier N-1
L = 2*math.pi
gap = L/N
x = np.linspace(0,L-gap,N)

# =============================================================================
# Evaluating the user defined function in the x array
# =============================================================================
c = (x-math.pi)**3
dc = 3*(x-math.pi)**2
#y = np.tanh(x-(L/2))
#y = (30/np.sqrt(4*math.pi))*np.exp(-((x-(L/2))**2)/(4))

# =============================================================================
# Computing the Fourier Matrix NxN from (0,N-1)
# =============================================================================
F = np.zeros((N,N),dtype=np.complex)
w = np.exp(complex(0,2*math.pi/N))
for i in range(0,N):
    for j in range(0,N):
        F[i,j]=w**(i*(j+N))

# =============================================================================
# Computing the first derivative Fourier Matrix NxN from (0,N-1)
# =============================================================================
DF = np.zeros((N,N),dtype=np.complex)
for i in range(0,N):
    for j in range(0,N):
        DF[i,j]=(j-N)*1j*w**(i*(j-N))
# =============================================================================
# Computing the complex Fourier Frequencies of the solution: Cmodal (Cm)
# =============================================================================
Cm = linalg.solve(DF,dc)
Cn = np.dot(F,Cm)










