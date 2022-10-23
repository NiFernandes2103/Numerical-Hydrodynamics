import numpy as np
from KTalgorithm import *

def E(q , C , dt):
    '''
    Foward Euler timestep

    '''
    E = q + dt * C

    return E


def modified_RungeKutta(q, C, dt, s=2):

    bl = 0.5
    y = np.array(s)
    E = E(q,C,dt)   # foward euler step

    y[0] = E

    for i in range(s):
        y[i] = bl * q + (1 - bl)*E(y[0],C,dt)

    return y[s]


def heun_integrator(t,y):
    pass