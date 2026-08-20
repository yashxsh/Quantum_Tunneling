import numpy as np 
import scipy 
import matplotlib.pyplot as plt 

m =1  #mass of particle
h = 1 #reduced planck constant 
             
V0 = 18.5  # Potential of the potential barrier 
L = 3   # length of potential barrier, the barrier will be from x = -L/2 to x = L/2
p = 6      
x0 = -40 
sigma = 0.6
div=1001   #number of discrete divisions(units) in the space
space=300
xmin = -space/2
xmax = space/2
    
xspace = np.linspace(xmin,xmax,div)   #building x-space
dx = xspace[1]-xspace[0]


#buikding initial Gaussian wave-function
psi0 = ((np.exp(-((xspace[1:-1]-x0)/sigma)**2))*np.exp(1j*p*xspace[1:-1]))
A = np.sum((abs(psi0)**2)*dx)  #norm of psi0
psi0 = psi0 / np.sqrt(A)      #normalization of psi0


V = np.where(np.abs(xspace) < (L/2), V0, 0)  #assigning potential to each point in x-space.

# creating hamiltonian 

scalar = (-h**2) / (2 * m * dx**2)
main_diag = (V[1:-1] * 2 * m * dx**2) / (-h**2) - 2
d = main_diag * scalar 
off_diag = np.ones(div - 3) * scalar


#finding eignevalues and eignevectors 

En, phi = scipy.linalg.eigh_tridiagonal(d, off_diag)   
phi = phi.T    #we transpose because linalg makes the eigenvectors as columns


#finding coefficients 
c = np.dot(np.conj(phi), psi0)    

#stevens
def frontx(psi):
    global xspace
    i = 0
    sum = 0 
    density = np.abs(psi)**2
    while sum<0.999 and i < len(density) - 1:
        sum+=density[i]*dx
        i+=1
    return xspace[1:-1][i]


t= 0 
dt = 0.035

t_value = []
x_value = []
for i in range(time_range:=800):
    time_evolution = c*np.exp(-1j*En*t/h)
    psi = np.dot(time_evolution,phi)
    t_value.append(t)
    x_value.append(frontx(psi))
    t+=dt

k = p/h 

steven_vel = h*k/m 

plt.plot([min(t_value),max(t_value)],[steven_vel*min(t_value) + x_value[0],steven_vel*max(t_value)+x_value[0]],label="steven")
plt.plot(t_value,x_value,label="simulation")

plt.xlabel("time")
plt.ylabel("position of 99.9th percentile ")
plt.ylim(-L/2,L/2)
plt.legend()
plt.show()