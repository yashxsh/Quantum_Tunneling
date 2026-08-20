from basic_functions import h,m,c,En,phi,dx,xspace,V,xmin,xmax,modpsisq
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









