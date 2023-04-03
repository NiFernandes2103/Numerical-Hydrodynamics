import numpy as np
from KTalgorithm import *

def Euler(q , C , dt):

    '''
    Foward Euler timestep
    '''
    E = q + dt * C

    return E


def modified_RungeKutta(q, C, dt, s=2):

    bl = 0.5

    A = q.shape[0]
    try:
      B = q.shape[1]
      y = np.empty([s,A,B])
    except:
      y = np.empty([s,A])

    y[0:] = Euler(q,C,dt) # foward euler step

    for i in range(s):
        y[i] = bl * q + (1 - bl)*Euler(y[0],C,dt)

    return y[s-1]

def Heuns(q,f,dt,t):

  k1 = dt*f
  k2 = dt*f(t + 1,q + k1)

  return q + 1/2 * (k1 + k2)

def RK4(y0,f,h,t):
  
  k1 = h * (f(t, y0))
  k2 = h * (f((t+h/2), (y0+k1/2)))
  k3 = h * (f((t+h/2), (y0+k2/2)))
  k4 = h * (f((t+h), (y0+k3)))
  k = (k1+2*k2+2*k3+k4)/6
  yn = y0 + k

  return yn


