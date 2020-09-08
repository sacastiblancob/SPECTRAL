#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Sergio Aldair Castiblanco Ballesteros
#Universidad Nacional de Colombia
#Ingeniería Civil

#Fourier transform and related things

#Fourier Transform
def FT(cn,N):
    
    import numpy as np
    
    L = 2*np.pi
    gap = L/N
    x = np.linspace(-np.pi,np.pi-gap,N)
    
    F = np.zeros((N,N),dtype=np.complex)
    if N%2==0:    
        k = np.array(range(int(-N/2),int((N/2))))
    else:
        k = np.array(range(int(-N/2),int((N/2)+1)))

    for i in range(0,N):
        for j in range(0,N):
            F[i,j]=(1/N)*np.exp(-1j*k[i]*x[j])
    
    cm = np.dot(F,cn)
    return cm

#Inverse Fourier transform
def IFT(cm,N):
    
    import numpy as np
    
    L = 2*np.pi
    gap = L/N
    x = np.linspace(-np.pi,np.pi-gap,N)
    
    F = np.zeros((N,N),dtype=np.complex)
    if N%2==0:    
        k = np.array(range(int(-N/2),int((N/2))))
    else:
        k = np.array(range(int(-N/2),int((N/2)+1)))

    for i in range(0,N):
        for j in range(0,N):
            F[i,j]=np.exp(1j*k[i]*x[j])
    
    cn = np.dot(F,cm)
    return cn
