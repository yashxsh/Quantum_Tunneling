import numpy as np 
import scipy 
m =1 
h = 1
             
V0 = 18.5
L = 0.5
p = 6 
x0 = -40 
sigma = 2.5 
div=1001
space=300
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

    

def modpsisq(psi,lower,higher):
    mask = (xspace[1:-1]<=higher) & (xspace[1:-1]>=lower)
    return np.sum(abs(psi[mask])**2) * dx 






