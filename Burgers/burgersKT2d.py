import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import integrate
from matplotlib.ticker import LinearLocator
from mpl_toolkits import mplot3d
import time
start_time = time.time()

'''
These are auxiliary functions in for the gradient 
'''

def minmod2(x,y):
  return (np.sign(x) + np.sign(y))*np.minimum(np.abs(x), np.abs(y))/2

def minmod3(x,y,z):
  return minmod2(x,minmod2(y,z))


def getGradient(f, dx, axis=0, theta=1):
    """
    Calculate the gradients of a field
    f        is a matrix of the field
    dx       is the cell size in x direction
    axis     is the axis of x-direction
    f_dx     is a matrix of derivative of f in the x-direction
    theta    is the flux limiter 1 <= theta <= 2
    
    """


    df_dx = np.zeros(f.shape)
    n = f.shape[axis]   
    K = np.arange(0, n)
    Kp1 = np.roll(K, -1)
    Km1 = np.roll(K, 1)

    if axis == 0:
      df_dx =  minmod3( theta * ( f - f[Km1][:, K] )/dx, (f[Kp1][:, K] - f[Km1][:, K] ) / (2*dx),theta * ( f[Kp1][:, K] - f ) / dx)
    elif axis == 1:
      df_dx = minmod3( theta * ( f - f[K][:, Km1] )/dx, (f[K][:, Kp1] - f[K][:, Km1] ) / (2*dx),theta * ( f[K][:, Kp1] - f ) / dx)

    return df_dx

def local_propagation_speed(u): 
  
   '''
    Get the local propagation speeds using the eigenvalues 
    of the flux matrix of the non relativistic IS equations

    u          is a matrix of density
   
    '''
    
   return u


def extrapolateInSpaceToFace(q, q_dx, dx, axis=0):
    """
    Calculate the gradients of a field
    q        is a matrix of the field
    q_dx     is a matrix of the field x-derivatives
    dx       is the cell size
    q_XL     is a matrix of spatial-extrapolated values on `left' face along x-axis 
    q_XR     is a matrix of spatial-extrapolated values on `right' face along x-axis 
    """

    n,_ = q.shape

    K = np.arange(0, n)
    Kp1 = np.roll(K, -1)
    Km1 = np.roll(K, 1)


    qP_XL = np.zeros_like(q)
    qP_XR = np.zeros_like(q)
    qM_XR = np.zeros_like(q)
    qM_XL = np.zeros_like(q)
       
    if axis == 0:   
      qP_XL = q - q_dx * dx/2
      qP_XR = q[Kp1][:, K] - q_dx[Kp1][:, K] * dx/2
      qM_XR = q + q_dx * dx/2
      qM_XL = q[Km1][:, K] + q_dx[Km1][:, K] * dx/2

    elif axis == 1:
      qP_XL = q - q_dx * dx/2
      qP_XR = q[K][:, Kp1] - q_dx[K][:, Kp1] * dx/2
      qM_XR = q + q_dx * dx/2
      qM_XL = q[K][:, Km1] + q_dx[K][:, Km1] * dx/2

    
    return qM_XL, qP_XL, qM_XR, qP_XR

def getFlux(u_P, u_M):

  """

  Calculate fluxes between 2 states with local Kurganov Tadmor rule 
  rho_P        is a matrix of left-state  density
  rho_M        is a matrix of right-state density

  """
  
  u = 0.5*(u_P + u_M)
  
  # compute fluxes (local Kurganov-Tadmor)

  flux_Mass   = 0.25*(u_P**2 + u_M**2)
  
  # find wavespeeds

  C_P = local_propagation_speed(u_P) # max propagation speed from the left

  C_M = local_propagation_speed(u_M) # max propagation speed from the right

  C = np.maximum(C_M, C_P)

  # add stabilizing diffusive term
  flux_Mass -= C * 0.5 * (u_P - u_M)

  return flux_Mass

def applyFluxes(flux_H1_X, flux_H2_X, flux_H1_Y, flux_H2_Y, dx, dy, J = 0):
    """
    Apply fluxes to conserved variables
    H         is a matrix of the conserved variable field
    flux_H1_X is a matrix of the x-dir fluxes from the right 
    flux_H2_X is a matrix of the x-dir fluxes from the left
    flux_H1_Y is a matrix of the y-dir fluxes from the right 
    flux_H2_Y is a matrix of the y-dir fluxes from the left    
    dx        is the cell size in the x direction
    dy        is the cell size in the y direction
    """
    C = 0

    # update solution
    C -= (flux_H1_X - flux_H2_X ) / dx
    C -= (flux_H1_Y - flux_H2_Y ) / dy
    C += J
    
    return C

def Burgers(t_init,q,dx,dy,X,Y,explicit_in_time=0):
  """ Burgers equation solve for each time step """

  #-----------------------------------------------------------------------------------------------------------------------------------#

  # Generate Initial Conditions  

  u = q

  '''
  u[0] = u[1]
  u[-1] = u[-2]
  '''

#-----------------------------------------------------------------------------------------------------------------------------------#

  # calculate gradients
  # getGradient(f, dx, axis, theta=1)

  u_dx = getGradient(u,  dx, 0)
  u_dy = getGradient(u,  dy, 1)

#-----------------------------------------------------------------------------------------------------------------------------------#

  # extrapolate in space to face centers

  uM_XL, uP_XL, uM_XR, uP_XR = extrapolateInSpaceToFace(u, u_dx, dx,0)
  uM_YL, uP_YL, uM_YR, uP_YR = extrapolateInSpaceToFace(u, u_dy, dy,1)

#-----------------------------------------------------------------------------------------------------------------------------------# 

  # compute fluxes (local Kurganov-Tadmor)

  flux_u_XR = getFlux(uP_XR, uM_XR)
  flux_u_XL = getFlux(uP_XL, uM_XL)
  flux_u_YR = getFlux(uP_YR, uM_YR)
  flux_u_YL = getFlux(uP_YL, uM_YL)

#-----------------------------------------------------------------------------------------------------------------------------------#

  # get time derivative

  dudt   = applyFluxes( flux_u_XR, flux_u_XL, flux_u_YR, flux_u_YL, dx, dy)

#-----------------------------------------------------------------------------------------------------------------------------------#

  return dudt

def Euler(q , C , dt):
    '''
    Foward Euler timestep

    q   is the field we are updating
    C   is the dq/dt
    dt  is the timestep


    '''
    E = q + dt * C

    return E


def modified_RungeKutta(q, C, dt, s=2):

    '''
    Modified Runge-Kutta integrator

    q   is the field we are updating
    C   is the dq/dt
    dt  is the timestep
    s   is the optional parameter for the order of the integrator

    This function updates the q field by one timestep using the numerical
    derivative of q over time (Heun's Method if s = 2)
    I did not implement higher order RK yet
    '''

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

def Heuns(q,f,dt,t,ExplicitTimeDependence):

  k1 = dt*f(t,q)
  if ExplicitTimeDependence:
    k2 = dt*f(t + dt,q + dt*k1)
  else:
    k2 = dt*(Euler(q + dt*k1,f(t,q + dt*k1),dt))

  return q + 1/2 * (k1 + k2)

def RK4(y0,f,h,t):
  
  k1 = h * (f(t, y0))
  k2 = h * (f((t+h/2), (y0+k1/2)))
  k3 = h * (f((t+h/2), (y0+k2/2)))
  k4 = h * (f((t+h), (y0+k3)))
  k = (k1+2*k2+2*k3+k4)/6
  yn = y0 + k

  return yn

def integrator(scheme, time, q0, dtmax, BC, method = "Heuns", args=None):

  '''
  This is an integrator that evolves a

  scheme     is the method to get dy/dt e.g. KTscheme
  time       is the current time
  q0         is the current state
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

  Q = [q0]
  q = q0
  N = int(args[2].shape[0])  

  while t < tEnd: 

    C = scheme

    # condition to ensure that the time steps are small enough so that
    # waves do not interfere with each other 
    courant_number = np.min(np.divide( dx, local_propagation_speed(q), out=np.zeros_like(dx*np.ones(q.shape)), where=local_propagation_speed(q)!=0))

    dt  =  np.minimum(dtmax, 0.5*courant_number) 

    if method == "Heuns":
      q = Heuns(q,C,dt,t,args[-1])
    if method == "RK4":
      q = RK4(q,C,dt,t)
    if method ==  "modified_RungeKutta":
      q = modified_RungeKutta(q,C(t,q),dt)

    #Apply Boundary conditions

    # BC(q)

    Q.append(q)

    t = t+dt
    
  return Q


# First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


# Simulation parameters
N                      = 100 # resolution
boxsize                = 4.  # in some unit system l
t                      = 0   # s 
tEnd                   = 0.1   # time at the end

# Define Mesh
dx = boxsize / N   # box size
dy = dx
xlin = np.linspace(-0.5*(boxsize-0.5*dx), 0.5*(boxsize-0.5*dx), N) # simulation limits
X,Y = np.meshgrid(xlin,xlin)    

R = np.sqrt(X**2 + Y**2)

u = 1*(R <= boxsize*0.25) + 0.125*(R > boxsize*0.5)

solution = integrator(Burgers, (t,tEnd), u, 0.01, None, args=(dx,dy,X,Y,True))

# Customize the z axis.
ax.set_zlim(0, 3.01)
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)


print(len(solution))

plt.imshow(solution[200])

surf = ax.plot_surface(X,Y,solution[0], cmap="plasma")
ax.view_init()

plt.show()
fig.canvas.draw()

print("--- %s seconds ---" % (time.time() - start_time))