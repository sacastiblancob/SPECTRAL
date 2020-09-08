from scipy import special
import math
import numpy as np

N1 = 17             #Number of points in the domain
N = N1-1            #Number of points in the domain of fourier N-1
L = 2*math.pi
gap = L/N
x = np.linspace(0,L-gap,N)
t = 3

b = np.zeros(N)
a = np.zeros(N,dtype=np.complex)  
for i in range(0,N):
    b[i] = special.jv(i+1,math.pi)
    a[i] = np.sin(i*math.pi/2)*b[i]*np.exp(1j*i*t)
