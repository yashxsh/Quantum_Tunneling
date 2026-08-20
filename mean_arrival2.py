import numpy as np 
import scipy 
import matplotlib.pyplot as plt 
m =1 
h = 1
           
def arrival_time(x,V0,L):
    
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


free_particle = arrival_time(50,0,0)
plt.plot([0,20],[free_particle,free_particle], label = "Free particle")
L_list = np.arange(0.5,20,0.5)
time = [] 
for i in L_list:
    time.append(arrival_time(50,V0:=19,i))

plt.plot(L_list,time, label = "With barrier")
plt.xlabel("Length of barrier")
plt.ylabel("arrival time")
plt.legend()
plt.show()