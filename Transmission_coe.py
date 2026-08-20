from basic_functions import h,m,c,En,phi,L,V0,p,xmax,modpsisq
import numpy as np 

t = 0 
dt = 0.035
for i in range(time_period:=500):
  time_evolution = c*np.exp(-1j*En*t/h)
  psi = np.dot(time_evolution,phi)
  t+=dt

T_numerical = modpsisq(psi,L/2,xmax)
print(T_numerical)