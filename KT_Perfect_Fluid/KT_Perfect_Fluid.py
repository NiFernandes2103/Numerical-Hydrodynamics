import numpy as np
from KT_relativistic_methods import *

def KTschemeRelativisticFluid(t,IC, dx, N, gamma, P0, theta=1, Eos = polytrope):

    """ Finite Volume simulation """
  
    #-----------------------------------------------------------------------------------------------------------------------------------#

    # Generate Initial Conditions  

    """ Initialize conservative variables """
    D  = IC[0:N]
    Sx = IC[N:2*N]
    tau  = IC[2*N:3*N]

    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # calculate gradients
    # getGradient(f, dx, axis=0, theta=1)

    D_dx = getGradient(D,  dx, 0, theta)
    Sx_dx  = getGradient(Sx,   dx, 0, theta)
    tau_dx = getGradient(tau,  dx, 0, theta)
    
    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # extrapolate in space to face centers
    # input extrapolateInSpaceToFace(q, q_dx, dx, axis=0)
    # output qM_XL, qP_XL, qM_XR, qP_XR

    DM_XL, DP_XL, DM_XR, DP_XR = extrapolateInSpaceToFace(D, D_dx, dx, 0)
    SxM_XL,  SxP_XL,  SxM_XR,  SxP_XR  = extrapolateInSpaceToFace(Sx,  Sx_dx,   dx, 0)
    tauM_XL,  tauP_XL,  tauM_XR,  tauP_XR  = extrapolateInSpaceToFace(tau,  tau_dx,   dx, 0)
    

    #-----------------------------------------------------------------------------------------------------------------------------------# 
    
    # compute fluxes (local Kurganov-Tadmor)
    # def getXFlux(rho_P, rho_M, vx_P, vx_M, vy_P, vy_M, eps_P, eps_M, P_P, P_M, gamma, P0, Eos=polytrope)
    # output flux_D, flux_Sx, flux_tau
 
    flux_D_XR, flux_Sx_XR, flux_tau_XR = getXFlux(DP_XR, DM_XR, SxP_XR, SxM_XR,
                                            tauP_XR, tauM_XR, gamma, P0, Eos)

    flux_D_XL, flux_Sx_XL, flux_tau_XL = getXFlux(DP_XL, DM_XL, SxP_XL, SxM_XL,
                                            tauP_XL, tauM_XL, gamma, P0, Eos)

    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # update solution

    timederivative_D = applyFluxes( flux_D_XR,   flux_D_XL, dx)
    timederivative_Sx = applyFluxes( flux_Sx_XR,   flux_Sx_XL, dx) 
    timederivative_tau = applyFluxes( flux_tau_XR,   flux_tau_XL, dx) 

    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    return np.hstack((timederivative_D,timederivative_Sx,timederivative_tau))
  

def integrator(scheme, time, q0, dtmax, BC, Eos=polytrope, method = "Heuns", args=None):

  """  
  This is an integrator that evolves a

  scheme     is the method to get dy/dt e.g. KTscheme that returns return 
    np.hstack((timederivative_rho,timederivative_Momx,timederivative_Momy,
    timederivative_Pixx,timederivative_Pixy,timederivative_Piyx,timederivative_Piyy))
  time       is the current time 
  q0         is the current state (np.vstack((rho,Momx,Momy,Pixx,Pixy,Piyx,Piyy))
  dtmax      is the upperbound of dt set by the user
  BC         is a function that enforces the boundary conditions
  method     is the method used in the integrator
  args       are additional arguments for scheme
  
  Later implement more Equations of State (Eos)
  """

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
        
  # write the parameters passed to the function
  t, tEnd = time
  Q = [q0]
  q = np.array(q0)
  dx = float(args[0])
  N = int(args[1])  
  outputCount = 1
  gamma   = args[2]
  P0      = args[3]


  while t < tEnd: 

    C = scheme
    courant_number = dx

    if (np.finfo(float).eps > courant_number):
      print("slow update")

    dt  =  np.minimum(dtmax, 0.4*courant_number) 
    

    # get conserved variables
  

    # choose the scheme to integrate(evolve over time) the system 
    if method == "Heuns":
      q = Heuns(q,C,dt,t)
    if method == "RK4":
      q = RK4(q,C,dt,t)
    if method == "Modified_RK":
      q = explicit_modified_RK(q,C,dt,t)
    
    
    # recover primitive variables
    rho,vx,eps = getPrimitive(q[0:N],q[N:2*N],q[2*N:3*N],gamma,P0,0.5*np.ones(N),100,Eos)
    P = Eos(rho,gamma,P0)
    p = np.vstack((rho,vx,P))

    #Apply Boundary conditions
    BC(q)

    t = t+dt
    
    if t >= dtmax*outputCount:
      Q.append(q)
      #M.append(np.sum(rho*args[0]*args[0]))
      print('{:.2f}/{:.2f}'.format(t,tEnd))
      outputCount += 1


  return Q


def BC(q):
  D = q[:N]
  Sx = q[N:2*N]
  tau = q[2*N:3*N]

  D[1]    = D[2]
  D[0]    = D[1]
  
  D[-2]   = D[-3]
  D[-1]   = D[-2]


  #Sx[0]    = -Sx[1]
  #Sx[-1]   = -Sx[-2]


  tau[1]  = tau[2]
  tau[0]  = tau[1]
  tau[-2] = tau[-3]
  tau[-1] = tau[-2]
  


  

t                      = 0    # s 
tEnd                   = 1    # time at the end
tOut                   = 0.01 # time of each output
N                      = 1000 # resolution
boxsize                = 4.   # in some unit system l
gamma                  = 1.4  # adiabatic index
P0                     = 1    # pressure constant
theta                  = 1    # flux limiter parameter


# Define Mesh
dx = boxsize / N   # box size
vol = dx**2        # volume of each box
a = (0.5*dx  - boxsize)*0.5
b = (boxsize - 0.5*dx)*0.5
xlin = np.linspace(a, b, N)# simulation limits

parameters = [t,tEnd,tOut,N,boxsize,gamma,theta,a,b]


""" initial condition of density""" # max = 1, min = 0

#rho = (1.5*(R1 <= 1) + 1*(R1 > 1))
rho = 0.5*(((1 - ((xlin)**2) )**4 )*(np.abs(xlin) < 1) + 0.1*np.ones(xlin.shape))# Mauricio`s funtion advice 
#rho = (1/(R))*(R>0)*(R<1) + 0.1*np.ones(s)
#rho = 0.1*(xlin >= 0) + 0.5*( xlin < 0)
#rho = np.ones(N) 

""" initial condition of velocity""" # max magnitude of 1
vx = np.zeros(N)
#vx = -1*np.sin(Theta)*(R < 1)
#vx = 0.9*(xlin < 0)

#vx = -0.99*(xlin < 0) + 0.99*(xlin >= 0)
#vx = np.abs((xlin - (boxsize-0.5*dx)*0.5)/16)

P = polytrope(rho,gamma,P0)

W = Lorentz_factor(vx)
print(np.max(W))
eps = P/3

D,Sx,tau = getConserved(rho,vx,eps,gamma,P0,polytrope)

IC = np.hstack((D,Sx,tau)) # here the initial conditions are stacked in a vector 
                            # rho is IC[0:N] for example

# input (dx, dy, xlin, gamma, BC, theta=1)
# output solution list of arrays that are 3N x N in the order (rho,vx,P)
solution = integrator(KTschemeRelativisticFluid, (t, tEnd), IC, 0.01, BC, polytrope, method="Modified_RK", args=(dx, N, gamma, P0, theta))

np.savetxt('PerfectFluidRelativistic_parameters',parameters)
np.save('PerfectFluidRelativistic',solution)


print(solution[-1])
