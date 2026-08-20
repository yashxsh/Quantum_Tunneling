from basic_functions import h,m,c,En,phi,dx,xspace,V,L,V0,p,xmax,xmin,modpsisq
import numpy as np 

dwell_time = 0
t = 0 
dt = 0.035
for i in range(time_period:=2000):
  time_evolution = c*np.exp(-1j*En*t/h)
  psi = np.dot(time_evolution,phi)
  t+=dt
  dwell_time+= modpsisq(psi,-L/2,L/2)*dt

print(dwell_time)