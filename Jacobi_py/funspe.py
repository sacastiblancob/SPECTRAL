# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:51:55 2020
USEFULL FUNCTIONS FOR COMPUTE LEGENDRE AND JACOBI POLINOMIALS OF ORDER N AND
VALUES OF alpha AND beta.

Note that the result of the jacobip(x,alpha,beta,N) function are the Jacobi's
polynomials normalized, that means they are orthonormal. The next expresion
is valid for every jacobi polynomial:
    integral(P(x)*P(x)dx) between(-1 and 1) = 1.0

In contrast the result of the legendre(x) function are the Legendre polynomials
not normalized, that means that they are not orthonormal. For legendre
normalized result, you should multiply every Li(x) by 1/sqrt(2/(2*i+1)) where
i is the order of the Legendre polynomial.

In the same way if you want to go from the normalized to the non-normalized
version of Legendre polynomials, you should multiply the jacobip result by
sqrt(2/(2*i+1)), where again, i is the order of the Legendre polynomial.

@author: SERGIO
"""




def legendre(x,N):
    import numpy as np
    #computing the size of x
    M = np.size(x)
    P = np.zeros((N+1,M))
    P[0,] = 1
    P[1,] = x
    for k in range(1,N):
        P[k+1,] = ((2*k+1)/(k+1))*x*P[k,] - (k/(k+1))*P[k-1,]
#for orthonormal version uncoment next 2 lines
#    for k in range(0,N+1):
#        P[k,] = P[k,]*np.sqrt((2*k+1)/2)
    return P

def legendrenorm(x,N):
    import numpy as np
    #computing the size of x
    M = np.size(x)
    P = np.zeros((N+1,M))
    P[0,] = 1/np.sqrt(2)
    P[1,] = np.sqrt(3/2)*x
    for k in range(1,N):
        an = np.sqrt(k**2/((2*k+1)*(2*k-1)))
        an1 = np.sqrt((k+1)**2/((2*(k+1)+1)*(2*(k+1)-1)))
        P[k+1,] = (x*P[k,] - an*P[k-1,])/an1
    
    return P

def jacobigq(alpha,beta,N):
    # Jan. S. Hetshaven adaptation from Matlab
    # Purpose: Compute the N’th order Gauss quadrature points, x,
    # and weights, w, associated with the Jacobi
    # polynomial, of type (alpha,beta) > -1 ( <> -0.5).
    
    import math
    import numpy as np
    from numpy import linalg as ln
    
    if (N==0):
        x=(alpha-beta)/(alpha+beta+2)
        w=2 
        return x,w 
    h1 = 2*np.arange(0,N+1)+alpha+beta
    J = (np.diag(-1/2*(alpha**2-beta**2)/(h1+2)/h1)+ np.diag(2/(h1[0:N]+2)*np.sqrt(np.arange(1,N+1)*(np.arange(1,N+1)+alpha+beta)*(np.arange(1,N+1)+alpha)*(np.arange(1,N+1)+beta)/(h1[0:N]+1)/(h1[0:N]+3)),1))
    if (alpha+beta<10*np.finfo(float).eps):
        J[0,0]=0.0
    J = J + J.T
    D,V = ln.eig(J)
    x = np.sort(D)
    w = (V[0,:])**2*2**(alpha+beta+1)/(alpha+beta+1)*math.gamma(alpha+1)*math.gamma(beta+1)/math.gamma(alpha+beta+1)
    w = w[np.argsort(D)]
    
    return x,w
    
def jacobip(x,alpha,beta,N):
    # Jan. S. Hetshaven adaptation from Matlab
    # function [P] = JacobiP(x,alpha,beta,N)
    # Purpose: Evaluate Jacobi Polynomials of type (alpha,beta) > -1
    # (alpha+beta <> -1) at points x for all orders to N and returns
    # P[N:length(xp))]
    # Note : They are normalized to be orthonormal.
    # Turn points into row if needed.
    
    import math
    import numpy as np
    
    xp = x
    dims = np.size(xp)
    PL = np.zeros((N+1,dims))
    gamma0 = 2**(alpha+beta+1)/(alpha+beta+1)*math.gamma(alpha+1)*math.gamma(beta+1)/math.gamma(alpha+beta+1);
    PL[0,] = 1.0/np.sqrt(gamma0)
    if (N==0):
        return PL
    gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0
    PL[1,] = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/np.sqrt(gamma1);
    if (N==1):
        return PL
    # Repeat value in recurrence.
    aold = 2/(2+alpha+beta)*np.sqrt((alpha+1)*(beta+1)/(alpha+beta+3))
    # Forward recurrence using the symmetry of the recurrence.
    for i in range(1,N):
        h1 = 2*i+alpha+beta
        anew = 2/(h1+2)*np.sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3))
        bnew = - (alpha**2-beta**2)/h1/(h1+2)
        PL[i+1,:] = 1/anew*( -aold*PL[i-1,:] + (xp-bnew)*PL[i,:])
        aold =anew
    return PL

def jacobigl(alpha,beta,N):
    # Jan. S. Hetshaven adaptation from Matlab
    # function [x] = JacobiGL(alpha,beta,N)
    # Purpose: Compute the N’th order Gauss Lobatto quadrature
    # points, x, associated with the Jacobi polynomial,
    # of type (alpha,beta) > -1 ( <> -0.5).
    
    import numpy as np
    import funspe as fs
    
    x = np.zeros(N+1)
    if(N==1):
        x[0] = -1.0
        x[1] = 1.0
        return x
    xint, w = fs.jacobigq(alpha+1,beta+1,N-2)
    x = np.concatenate(([-1.0],xint,[1.0]),axis=None)
    
    return x
    







    