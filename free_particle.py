import numpy as np 
import scipy 
m =1  #mass of particle
h = 1 #reduced planck constant 
             
V0 = 0  # Potential of the potential barrier 
L = 0.5    # length of potential barrier, the barrier will be from x = -L/2 to x = L/2
p = 6      
x0 = -40 
sigma = 2.5 
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

def modpsisq(psi,lower,higher):
    mask = (xspace[1:-1]<=higher) & (xspace[1:-1]>=lower)
    return np.sum(abs(psi[mask])**2) * dx 

import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation 
import numpy as np 


def animation(time_scale=800,dt = 0.035):
    global c,En,h,phi
    final_soln =[]
    t = 0
    for i in range(time_scale):
        time_evolution = c*np.exp(-1j*En*t/h)
        psi = np.dot(time_evolution,phi)
        final_soln.append(psi)
        t+=dt

    return final_soln 

final_soln = animation()

fig,axis = plt.subplots()
axis.set_xlim(-100,100)
axis.set_ylim(0,0.3)
axis.set_xlabel("x")
axis.set_ylabel(r'$|\psi|^2$')
probab, = axis.plot([],[])


probability_list = []

def update_frames(frame):
    global c,En,h,phi,final_soln,xspace
    probab.set_data(xspace[1:-1],abs(final_soln[frame])**2)
    probability_list.append(modpsisq(final_soln[frame],xmin,xmax))
    return probab,


anim = FuncAnimation(fig = fig , func=update_frames,frames=len(final_soln),interval=10,repeat=True)
plt.plot(xspace,V)
anim.save("barrier.gif",writer="Pillow")
plt.show()
print("probability=",probability_list[0])
print("maximum relative error=",(max(abs(np.array(probability_list)))-probability_list[0])/probability_list[0])