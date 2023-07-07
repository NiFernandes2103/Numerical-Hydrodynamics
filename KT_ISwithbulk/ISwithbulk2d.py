import numpy as np
import matplotlib.pyplot as plt
from Ktmethods2d import *

def KTschemeNonRelativisticIS(t, IC, dx, dy, xlin, gamma, zeta, tau_nu, theta=1):

    """ Finite Volume simulation """

    #-----------------------------------------------------------------------------------------------------------------------------------#

    # Generate Initial Conditions  

    ''' Initial conditions for rho '''
    rho = IC[0:N]

    ''' Initial conditions for v'''
    vx = np.divide(IC[N:2*N] , rho, out=np.zeros_like(IC[N:2*N]), where=rho!=0)

    vy = np.divide(IC[2*N:3*N] , rho, out=np.zeros_like(IC[2*N:3*N]), where=rho!=0)

    ''' Pressure due to equation of state '''
    P = (np.abs(rho))**gamma

    ''' B Constant '''
    B = zeta/tau_nu

    ''' Pi initial condition '''
    Pi = IC[3*N:]
   

    #-----------------------------------------------------------------------------------------------------------------------------------#

    # get Conserved variables
    vol = dx*dx
    Mass, Momx, Momy = getConserved( rho, vx, vy, vol)


    # get Primitive variables
    rho, vx, vy, P = getPrimitive( Mass, Momx, Momy, gamma, vol )


    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # calculate gradients
    # getGradient(f, dx, theta=1)
    
    rho_dx = getGradient(rho, dx, 0, theta)
    vx_dx  = getGradient(vx,  dx, 0, theta)
    vy_dx  = getGradient(vy,  dx, 0, theta)
    P_dx   = getGradient(P,   dx, 0, theta)
    Pi_dx  = getGradient(Pi,  dx, 0, theta)

    rho_dy = getGradient(rho, dy, 1, theta)
    vx_dy  = getGradient(vx,  dy, 1, theta)
    vy_dy  = getGradient(vy,  dy, 1, theta)
    P_dy   = getGradient(P,   dy, 1, theta)
    Pi_dy  = getGradient(Pi,  dy, 1, theta)

    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # extrapolate in space to face centers
    # input extrapolateInSpaceToFace(q, q_dx, dx)
    # output qM_XL, qP_XL, qM_XR, qP_XR

    rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR = extrapolateInSpaceToFace(rho, rho_dx, dx, 0)
    vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR  = extrapolateInSpaceToFace(vx,  vx_dx,   dx, 0)
    vyM_XL,  vyP_XL,  vyM_XR,  vyP_XR  = extrapolateInSpaceToFace(vy,  vy_dx,   dx, 0)
    PM_XL,   PP_XL,   PM_XR,   PP_XR   = extrapolateInSpaceToFace(P,   P_dx,   dx, 0)
    PiM_XL,  PiP_XL,  PiM_XR,  PiP_XR  = extrapolateInSpaceToFace(Pi,  Pi_dx,  dx, 0)


    rhoM_YL, rhoP_YL, rhoM_YR, rhoP_YR = extrapolateInSpaceToFace(rho, rho_dy, dy, 1)
    vxM_YL,  vxP_YL,  vxM_YR,  vxP_YR  = extrapolateInSpaceToFace(vx,  vx_dy,  dy, 1)
    vyM_YL,  vyP_YL,  vyM_YR,  vyP_YR  = extrapolateInSpaceToFace(vy,  vy_dy,  dy, 1)
    PM_YL,   PP_YL,   PM_YR,   PP_YR   = extrapolateInSpaceToFace(P,   P_dy,   dy, 1)
    PiM_YL,  PiP_YL,  PiM_YR,  PiP_YR  = extrapolateInSpaceToFace(Pi,  Pi_dy,  dy, 1)

    #-----------------------------------------------------------------------------------------------------------------------------------# 
    
    # compute fluxes (local Kurganov-Tadmor)
    # input getFlux(rho_P, rho_M, vx_P, vx_M, Pi_P, Pi_M, P_P, P_M, gamma, B)
    # output flux_Mass, flux_Momx, flux_Pi_v

    flux_Mass_XR, flux_Momx_XR, flux_Momy_XR, flux_Pi_vxR = getXFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR,
                                                vyP_XR, vyM_XR, PiP_XR, PiM_XR, PP_XR, PM_XR, gamma, B)

    flux_Mass_XL, flux_Momx_XL, flux_Momy_XL, flux_Pi_vxL = getXFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL,
                                            vyP_XL, vyM_XL, PiP_XL, PiM_XL, PP_XL, PM_XL, gamma, B)

    flux_Mass_YR, flux_Momx_YR, flux_Momy_YR, flux_Pi_vyR = getYFlux(rhoP_YR, rhoM_YR, vxP_YR, vxM_YR,
                                            vyP_YR, vyM_YR, PiP_YR, PiM_YR, PP_YR, PM_YR, gamma, B)

    flux_Mass_YL, flux_Momx_YL, flux_Momy_YL, flux_Pi_vyL = getYFlux(rhoP_YL, rhoM_YL, vxP_YL, vxM_YL,
                                            vyP_YL, vyM_YL, PiP_YL, PiM_YL, PP_YL, PM_YL, gamma, B)

    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # get time derivative 

    J = - Pi/tau_nu
    timederivative_rho = applyFluxes(  flux_Mass_XR,   flux_Mass_XL, flux_Mass_YR,   flux_Mass_YL,  dx,  dy)
    timederivative_Momx = applyFluxes( flux_Momx_XR,   flux_Momx_XL, flux_Momx_YR,   flux_Momx_YL,  dx,  dy)
    timederivative_Momy = applyFluxes( flux_Momy_XR,   flux_Momy_XL, flux_Momy_YR,   flux_Momy_YL,  dx,  dy) 
    timederivative_Pi  = applyFluxes(  flux_Pi_vxR,    flux_Pi_vxL,  flux_Pi_vyR,    flux_Pi_vyL,   dx,  dy, J)

    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    return np.vstack((timederivative_rho, timederivative_Momx, timederivative_Momy, timederivative_Pi))

def Heuns(q,f,dt,t):

  k1 = dt*f(t,q)
  k2 = dt*f(t + dt,q + k1)

  return q + 1/2 * (k1 + k2)

def HeunswithFowardEuler(q,f,dt,t):

  k1 = dt*f(t,q)
  k2 = dt*(f(t,q) +  dt*f(t,q + k1))

  return q + 0.5 * (k1 + k2)

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
  N = int(args[2].shape[0])  
  while t < tEnd: 

    C = scheme

    if args is not None:

      rho = y[0:N]
      vx  = y[N:2*N]/y[0:N]
      vy  = y[2*N:3*N]/y[0:N]
      Pi  = y[3*N:]
      gamma = args[2]

      B =  args[4]/args[5]

      mlpsx = np.max(local_propagation_speed(rho, vx, Pi, gamma, B))
      mlpsy = np.max(local_propagation_speed(rho, vy, Pi, gamma, B))

      mlps = max(mlpsx,mlpsy)

      if mlps != 0:

        dtlocal = 0.3*args[0] / mlps

      else:

        dtlocal = 0

    if dtlocal > 0:
      dt  = np.minimum(dtmax,dtlocal) 
    else:
      dt = dtmax

    dt = max(dt,10**6)


    if method == "Heuns":
      y = Heuns(y,C,dt,t)
    if method == "RK4":
      y = RK4(y,C,dt,t)
    if method ==  "modified_RungeKutta":
      y = modified_RungeKutta(y,C(t,y),dt)

    #Apply Boundary conditions

    BC(y)

    Y.append(y)
    print('t=',t)
    t = t + dt
    
  return Y




def applyBC(y):

  rho = y[0:N] 
  vx = y[N:2*N]/rho 
  Pi = y[2*N:] 

  #Absorbing boundary conditions
  
  
  rho[0]    = rho[1]
  rho[-1]   = rho[-2]


  vx[0]     = vx[1]
  vx[1]     = vx[-2]      

  Pi[0]     = 0   
  Pi[-1]    = 0




t                      = 0   # s 
tEnd                   = 0.1   # time at the end
tOut                   = 0.01 # time of each output

N                      = 300 # resolution
boxsize                = 1.  # in some unit system l
gamma                  = 1 # adiabatic index
zeta                   = 1 # bulk viscosity coefficient
eta                    = 1
tau_nu                 = 1
theta                  = 1


# Define Mesh
dx = boxsize / N   # box size
dy = dx
vol = dx**2        # volume of each box
xlin = np.linspace(-0.5*(boxsize-0.5*dx), 0.5*(boxsize-0.5*dx), N)# simulation limits
Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
s = X.shape
R = np.sqrt(X**2 + Y**2)

# initial condition of density

rho = (1*(R <= 1*0.25) + 0.25*(R > 1*0.25))
plt.imshow(rho)
plt.show()
#rho = 1*(X < boxsize*0.5) + 0.125*(X >= boxsize*0.5)

# initial condition of velocity
vx = np.zeros(s)
#vx = 0.5*np.ones(xlin.shape)
#vx = 1*(X < boxsize*0.5) + 0*(X >= boxsize*0.5)
#vx = np.abs((xlin - (boxsize-0.5*dx)*0.5)/16)

vy = np.zeros(s)
#vy = 0.5*np.ones(xlin.shape)

# initial condition of Pi scalar
Pi = np.zeros(s)


IC = np.vstack((rho, rho*vx, rho*vy, Pi)) # here the initial conditions are stacked in a vector 
                            # rho is IC[0:N] for example

# dx, dy, xlin, gamma, zeta, tau_nu, BC, theta=1

# solution = integrator(KTschemeNonRelativisticIS, (t, tEnd), IC, 0.01, method="RK4", args=(dx, dy, xlin, gamma, zeta, tau_nu, eta, theta))


solution = integrator(KTschemeNonRelativisticIS, (t,tEnd), IC, tOut, applyBC, args=(dx, dy, xlin, gamma, zeta, tau_nu, theta))



i=0
while i < len(solution):
  rho = solution[i][:N][int(N/2)]
  ux = solution[i][N:2*N][int(N/2)]
  uy = solution[i][2*N:3*N][int(N/2)]
  Pi = solution[i][3*N:][int(N/2)]
  ur = np.sqrt(ux**2 + uy**2)
  uphy = (X[:][int(N/2)]*uy - Y[:][int(N/2)]*ux)/(R[:][int(N/2)]**2)
  plt.plot(rho)
  plt.show()
  plt.plot(uphy)
  plt.show()
  i += 10