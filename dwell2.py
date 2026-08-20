from basic_functions import h,m,c,En,phi,dx,xspace,psi0,p,V0,xmax,xmin,div,modpsisq
import numpy as np 
import matplotlib.pyplot as plt
import scipy  

def dwell(L):
    V = np.where(np.abs(xspace) <= (L/2), V0, 0) 
    scalar = (-h**2) / (2 * m * dx**2)
    main_diag = (V[1:-1] * 2 * m * dx**2) / (-h**2) - 2
    d = main_diag * scalar 
    off_diag = np.ones(div - 3) * scalar

    En, phi = scipy.linalg.eigh_tridiagonal(d, off_diag)
    phi = phi.T 
    c = np.dot(np.conj(phi), psi0) 
    dwell_time = 0
    t = 0 
    dt = 0.035
    for i in range(time_period:=500):
      time_evolution = c*np.exp(-1j*En*t/h)
      psi = np.dot(time_evolution,phi)
      t+=dt
      dwell_time+= modpsisq(psi,-L/2,L/2)*dt
    return dwell_time

L = np.arange(0,20,0.5)
d_time = []
for i in L:
   d_time.append(dwell(i))

plt.plot(L,d_time)
plt.xlabel("Length of barrier")
plt.ylabel("dwell time for V0 = 18.5")
plt.show()
    



