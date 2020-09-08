#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Sergio Aldair Castiblanco Ballesteros
#Universidad Nacional de Colombia
#Ingeniería Civil

#Fourier transform and related things

#Fourier Transform
def FT(cn):
    
    import numpy as np
    
    N = np.size(cn)
    
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
def IFT(cm):
    
    import numpy as np
    
    N = np.size(cm)
    
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

def OEP(a):     #Odd-Even Permutation (OEP)
    
    import numpy as np
    
    N = np.size(a)
    b = np.zeros((N))
    for i in range(0,int(N/2)):
        b[i] = a[2*i]
        b[i+int(N/2)] = a[2*i+1] 
    
    return b
    
def MOEP(cn):     #Multi Odd-Even Permutation (MOEP)
    
    import numpy as np
    from DFT import OEP
    
    N = np.size(cn)
    ca = cn
    if N%2 == 0:
        op = 0
        a = int(N)
        while a%2 == 0:
            cont = np.array(range(0,a))
            for i in range(0,int(N/a)):
                ca[cont+i*a] = OEP(ca[cont+i*a])
            a = int(a/2)
            op = op + 1
    else:
        op = 0
    return ca,op
    
def FFT(cn):
    
    import numpy as np
    from DFT import FT
    from DFT import MOEP
    
    N = np.size(cn)
    ca,op = MOEP(cn)
    
    if op ==0:
        cm = FT(cn)
    
    else:
        ...transformada rapida
    