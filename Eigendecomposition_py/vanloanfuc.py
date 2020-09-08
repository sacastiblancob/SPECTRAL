# -*- coding: utf-8 -*-
"""
Created on Wed May 13 13:49:25 2020

@author: SERGIO
"""


def house(x):
    import numpy as np
    n = len(x)
    sig = np.dot(x[1:n],x[1:n])
    v = np.append(1,x[1:n])
    
    if (sig==0):
        b = 0
    else:
        mu = np.sqrt(x[0]+sig)
        if (x[0]<=0):
            v[0] = x[0] - mu
        else:
            v[0] = -sig/(x[0] + mu)
        b = (2*v[0]**2)/(sig + v[0]**2)
        v = v/v[0]
    return v,b