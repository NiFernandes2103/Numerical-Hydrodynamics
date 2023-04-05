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

  k1 = dt*f(t,q)
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

  
def integrator(scheme, time, y0, dtmax, BC, method = "Heuns", args=None):

  '''
  This is an integrator that evolves a

  scheme     is the method to get dy/dt e.g. KTscheme
  time       is the current time
  y0         is the current state
  dtmax      is the upperbound of dt set by the user
  BC         is a function that enforces the boundary conditions
  method     is the method used in the integrator
  args       are additional arguments for scheme
  '''

  if args is not None:
        # Wrap the user's scheme in lambdas to hide the
        # additional parameters.  Pass in the original fun as a keyword
        # argument to keep it in the scope of the lambda.
        try:
            _ = [*(args)]
        except TypeError as exp:
            suggestion_tuple = (
                "Supplied 'args' cannot be unpacked. Please supply `args`"
                f" as a tuple (e.g. `args=({args},)`)")
            raise TypeError(suggestion_tuple) from exp

        scheme = lambda t, x, scheme = scheme: scheme(t, x, *args)
        
  t, tEnd = time

  Y = [y0]
  y = y0
  N = args[1].shape[0]
  
  while t < tEnd: 

    C = scheme

    if args is not None:
      dtlocal = np.min(args[0]/local_propagation_speed(y[0:N],y[N:2*N]/y[0:N],y[2*N:], args[2], args[3]/args[4]))
    dt  = np.minimum(dtmax,dtlocal) 

    if method == "Heuns":
      y = Heuns(y,C,dt,t)
    if method == "RK4":
      y = RK4(y,C,dt,t)
    
    #Apply Boundary conditions

    BC(y)

    Y.append(y)

    t = t+dt
    
  return Y


