import numpy as np 
import scipy 
import matplotlib.pyplot as plt 
m =1 
h = 1
           
def peak_trace(V0,L):

    t_value = []
    x_value = []
    
    p = 6 
    x0 = -30 
    sigma = 8 
    div=1001
    space=200
    xmin = -space/2
    xmax = space/2
    
    xspace = np.linspace(xmin,xmax,div) #building x-space
    dx = xspace[1]-xspace[0]
    V = np.where(np.abs(xspace) < (L/2), V0, 0)

    # creating hamiltonian 
    scalar = (-h**2) / (2 * m * dx**2)
    main_diag = (V[1:-1] * 2 * m * dx**2) / (-h**2) - 2
    d = main_diag * scalar 
    off_diag = np.ones(div - 3) * scalar


    #finding eignevalues and eignevectors 

    En, phi = scipy.linalg.eigh_tridiagonal(d, off_diag)
    phi = phi.T

    psi0 = ((np.exp(-((xspace[1:-1]-x0)/sigma)**2))*np.exp(1j*p*xspace[1:-1]))
    A = np.sum((abs(psi0)**2)*dx)
    psi0 = psi0 / np.sqrt(A)

    #finding coefficients 
    c = np.dot(np.conj(phi), psi0)
    t = 0 
    dt = 0.035
    for i in range(time_scale:=800):
        time_evolution = c*np.exp(-1j*En*t/h)
        psi = np.dot(time_evolution,phi)
        t_value.append(t)
        if V0==0:
           x_value.append(xspace[1:-1][np.argmax(np.abs(psi)**2)])
        else:
            check = np.sum(np.abs(psi[xspace[1:-1]>L/2])**2)
            if check >= 0.000001:    
                mask = xspace[1:-1] > L/2
                psi_trans = psi[mask]
                x_trans = xspace[1:-1][mask]
                peak_idx = np.argmax(np.abs(psi_trans)**2)
                x_value.append(x_trans[peak_idx])
            else:
                x_value.append(np.nan)
        t+=dt
    return t_value,x_value

import matplotlib.pyplot as plt 

free_particle = peak_trace(0,3)
barrier = peak_trace(18.5,3)
plt.xlabel("time")
plt.ylabel("position of the peak")
plt.plot(free_particle[0],free_particle[1],label="free particle")
plt.plot(barrier[0],barrier[1],label="barrier")
plt.legend()
plt.show()
    