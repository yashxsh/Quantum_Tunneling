import numpy as np 
from basic_functions import c,En,phi,h,xspace


def time_to_reach(x):
    global c,En,phi,h,xspace
    t = 0  
    dt = 0.035
    nr = 0 
    dr = 0 
    idx = np.argmin(np.abs(xspace - x))
    for i in range(time_scale:=800):
        time_evolution = c*np.exp(-1j*En*t/h)
        psi = np.dot(time_evolution,phi)
        nr+= t*abs(psi[idx])**2
        dr+= abs(psi[idx])**2 
        t+=dt 
    
    return nr/dr



detector = 50 
print(time_to_reach(detector))

        



