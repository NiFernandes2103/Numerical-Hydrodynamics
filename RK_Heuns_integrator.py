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


def heun_integrator(t,y):
    pass