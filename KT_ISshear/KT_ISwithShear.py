import numpy as np
from KTmethods2d import *

def KTschemeNonRelativisticIS(t,IC, dx, dy, N, gamma, zeta, tau_nu, eta, P0, theta=1):

    """ Finite Volume simulation """
  
    #-----------------------------------------------------------------------------------------------------------------------------------#

    # Generate Initial Conditions  

    ''' Initial conditions for rho '''
    rho = IC[0:N]


    ''' Initial conditions for v'''

    Momx = IC[N:2*N]
    Momy = IC[2*N:3*N]
    vx = IC[N:2*N] / rho
    vy = IC[2*N:3*N] / rho


    ''' Pressure due to equation of state '''
    P = P0*rho**gamma
    

    ''' Pi initial condition '''
    Pixx = IC[3*N:4*N]
    Pixy = IC[4*N:5*N]
    Piyx = IC[5*N:6*N]
    Piyy = IC[6*N:]


    # get Speed of sound
    cs = getSpeedOfSound(rho, gamma, P0)

    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # calculate gradients
    # getGradient(f, dx, axis=0, theta=1)

    rho_dx = getGradient(rho,  dx, 0, theta)
    vx_dx  = getGradient(vx,   dx, 0, theta)
    vy_dx  = getGradient(vy,   dx, 0, theta)
    P_dx   = getGradient(P,    dx, 0, theta)
    Pixx_dx  = getGradient(Pixx,   dx, 0, theta)
    Pixy_dx  = getGradient(Pixy,   dx, 0, theta)
    Piyx_dx  = getGradient(Piyx,   dx, 0, theta)
    Piyy_dx  = getGradient(Piyy,   dx, 0, theta)

    rho_dy = getGradient(rho,  dy, 1, theta)
    vx_dy  = getGradient(vx,   dy, 1, theta)
    vy_dy  = getGradient(vy,   dy, 1, theta)
    P_dy   = getGradient(P,    dy, 1, theta)
    Pixx_dy  = getGradient(Pixx,   dy, 1, theta)
    Pixy_dy  = getGradient(Pixy,   dx, 1, theta)
    Piyx_dy  = getGradient(Piyx,   dx, 1, theta)
    Piyy_dy  = getGradient(Piyy,   dx, 1, theta)

    
    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # extrapolate in space to face centers
    # input extrapolateInSpaceToFace(q, q_dx, dx, axis=0)
    # output qM_XL, qP_XL, qM_XR, qP_XR

    rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR = extrapolateInSpaceToFace(rho, rho_dx, dx, 0)
    vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR  = extrapolateInSpaceToFace(vx,  vx_dx,   dx, 0)
    vyM_XL,  vyP_XL,  vyM_XR,  vyP_XR  = extrapolateInSpaceToFace(vy,  vy_dx,   dx, 0)
    PixxM_XL,  PixxP_XL,  PixxM_XR,  PixxP_XR  = extrapolateInSpaceToFace(Pixx,  Pixx_dx,  dx, 0)
    PixyM_XL,  PixyP_XL,  PixyM_XR,  PixyP_XR  = extrapolateInSpaceToFace(Pixy,  Pixy_dx,  dx, 0)
    PiyxM_XL,  PiyxP_XL,  PiyxM_XR,  PiyxP_XR  = extrapolateInSpaceToFace(Piyx,  Piyx_dx,  dx, 0)
    PiyyM_XL,  PiyyP_XL,  PiyyM_XR,  PiyyP_XR  = extrapolateInSpaceToFace(Piyy,  Piyy_dx,  dx, 0)

    rhoM_YL, rhoP_YL, rhoM_YR, rhoP_YR = extrapolateInSpaceToFace(rho, rho_dy, dy, 1)
    vxM_YL,  vxP_YL,  vxM_YR,  vxP_YR  = extrapolateInSpaceToFace(vx,  vx_dy,  dy, 1)
    vyM_YL,  vyP_YL,  vyM_YR,  vyP_YR  = extrapolateInSpaceToFace(vy,  vy_dy,  dy, 1)
    PixxM_YL,  PixxP_YL,  PixxM_YR,  PixxP_YR  = extrapolateInSpaceToFace(Pixx,  Pixx_dy,  dy, 1)
    PixyM_YL,  PixyP_YL,  PixyM_YR,  PixyP_YR  = extrapolateInSpaceToFace(Pixy,  Pixy_dy,  dy, 1)
    PiyxM_YL,  PiyxP_YL,  PiyxM_YR,  PiyxP_YR  = extrapolateInSpaceToFace(Piyx,  Piyx_dy,  dy, 1)
    PiyyM_YL,  PiyyP_YL,  PiyyM_YR,  PiyyP_YR  = extrapolateInSpaceToFace(Piyy,  Piyy_dy,  dy, 1)

    #-----------------------------------------------------------------------------------------------------------------------------------# 
    
    # compute fluxes (local Kurganov-Tadmor)
    # getXFlux(rho_P, rho_M, vx_P, vx_M, vy_P, vy_M, Pixx_P, Pixx_M, Pixy_P,
    #          Pixy_M, Piyx_P,Piyx_M, Piyy_P, Piyy_M, P_P, P_M, gamma, eta,
    #          zeta, tau_nu):
    # getYFlux(rho_P, rho_M, vx_P, vx_M, vy_P, vy_M, Pixx_P, Pixx_M, Pixy_P,
    #          Pixy_M, Piyx_P,Piyx_M, Piyy_P, Piyy_M, P_P, P_M, gamma, eta,
    #          zeta, tau_nu):
    # output flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vx, flux_Pixy_vx, flux_Piyx_vx, flux_Piyy_vx
    # output flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vy, flux_Pixy_vy, flux_Piyx_vy, flux_Piyy_vy
 
    

    flux_Mass_XR, flux_Momx_XR, flux_Momy_XR, flux_Pixx_vxR, flux_Pixy_vxR, flux_Piyx_vxR, flux_Piyy_vxR = getXFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR,
                                            vyP_XR, vyM_XR, PixxP_XR, PixxM_XR,
                                            PixyP_XR, PixyM_XR, PiyxP_XR, PiyxM_XR,
                                            PiyyP_XR, PiyyM_XR, gamma, eta, zeta, tau_nu, P0)

    flux_Mass_XL, flux_Momx_XL, flux_Momy_XL, flux_Pixx_vxL, flux_Pixy_vxL, flux_Piyx_vxL, flux_Piyy_vxL = getXFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL,
                                            vyP_XL, vyM_XL, PixxP_XL, PixxM_XL,
                                            PixyP_XL, PixyM_XL, PiyxP_XL, PiyxM_XL,
                                            PiyyP_XL, PiyyM_XL, gamma, eta, zeta, tau_nu, P0)

    flux_Mass_YR, flux_Momx_YR, flux_Momy_YR, flux_Pixx_vyR, flux_Pixy_vyR, flux_Piyx_vyR, flux_Piyy_vyR = getYFlux(rhoP_YR, rhoM_YR, vxP_YR, vxM_YR,
                                            vyP_YR, vyM_YR, PixxP_YR, PixxM_YR,
                                            PixyP_YR, PixyM_YR, PiyxP_YR, PiyxM_YR,
                                            PiyyP_YR, PiyyM_YR, gamma, eta, zeta, tau_nu, P0)

    flux_Mass_YL, flux_Momx_YL, flux_Momy_YL, flux_Pixx_vyL, flux_Pixy_vyL, flux_Piyx_vyL, flux_Piyy_vyL = getYFlux(rhoP_YL, rhoM_YL, vxP_YL, vxM_YL,
                                            vyP_YL, vyM_YL, PixxP_YL, PixxM_YL,
                                            PixyP_YL, PixyM_YL, PiyxP_YL, PiyxM_YL,
                                            PiyyP_YL, PiyyM_YL, gamma, eta, zeta, tau_nu, P0)
    
    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # update solution

    Jxx = - Pixx/tau_nu
    Jxy = - Pixy/tau_nu
    Jyx = - Piyx/tau_nu
    Jyy = - Piyy/tau_nu


    timederivative_rho = applyFluxes( flux_Mass_XR,   flux_Mass_XL, flux_Mass_YR,   flux_Mass_YL,  dx, dy)
    timederivative_Momx = applyFluxes( flux_Momx_XR,   flux_Momx_XL, flux_Momx_YR,   flux_Momx_YL, dx, dy) 
    timederivative_Momy = applyFluxes( flux_Momy_XR,   flux_Momy_XL, flux_Momy_YR,   flux_Momy_YL, dx, dy) 
    timederivative_Pixx  = applyFluxes( flux_Pixx_vxR,    flux_Pixx_vxL, flux_Pixx_vyR,    flux_Pixx_vyL, dx, dy, Jxx)
    timederivative_Pixy  = applyFluxes( flux_Pixy_vxR,    flux_Pixy_vxL, flux_Pixy_vyR,    flux_Pixy_vyL, dx, dy, Jxy)
    timederivative_Piyx  = applyFluxes( flux_Piyx_vxR,    flux_Piyx_vxL, flux_Piyx_vyR,    flux_Piyx_vyL, dx, dy, Jyx)
    timederivative_Piyy  = applyFluxes( flux_Piyy_vxR,    flux_Piyy_vxL, flux_Piyy_vyR,    flux_Piyy_vyL, dx, dy, Jyy)


    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    return np.vstack((timederivative_rho,timederivative_Momx,timederivative_Momy,timederivative_Pixx,timederivative_Pixy,timederivative_Piyx,timederivative_Piyy))
  

def integrator(scheme, time, q0, dtmax, BC, method = "Heuns", args=None):

  '''
  This is an integrator that evolves a

  scheme     is the method to get dy/dt e.g. KTscheme that returns return 
    np.vstack((timederivative_rho,timederivative_Momx,timederivative_Momy,
    timederivative_Pixx,timederivative_Pixy,timederivative_Piyx,timederivative_Piyy))
  time       is the current time 
  q0         is the current state (np.vstack((rho,Momx,Momy,Pixx,Pixy,Piyx,Piyy))
  dtmax      is the upperbound of dt set by the user
  BC         is a function that enforces the boundary conditions
  method     is the method used in the integrator
  args       are additional arguments for scheme
  
  Later implement more Equations of State (Eos)
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
        
  # write the parameters passed to the function
  t, tEnd = time
  Q = [q0]
  q = q0
  N = int(args[2])  
  outputCount = 1
  gamma   = args[3]
  zeta    = args[4]
  tau_nu  = args[5]
  eta     = args[6]
  P0      = args[7]

  t0 = t
  while t < tEnd: 

    C = scheme

    # compute density to get necessary stability conditions met
    rho = q[0:N]

    #get speed of sound
    cs = getSpeedOfSound(rho,gamma,P0)

    # condition to ensure that the time steps are small enough so that
    # waves do not interfere with each other 
    # rho, eta, zeta, tau_nu, cs
    courant_number = dx/np.max(local_propagation_speed(rho, eta, zeta, tau_nu, cs))

    if (np.finfo(float).eps > courant_number):
      print("slow update")

    dt  =  np.minimum(dtmax, 0.1*courant_number) 

    # choose the scheme to integrate(evolve over time) the system 
    if method == "Heuns":
      q = Heuns(q,C,dt,t)
    if method == "RK4":
      q = RK4(q,C,dt,t)
    if method == "Modified_RK":
      q = explicit_modified_RK(q,C,dt,t)

    #Apply Boundary conditions

    BC(q)


    t = t+dt
    
    if t >= dtmax*outputCount + t0:
      Q.append(q)
      #M.append(np.sum(rho*args[0]*args[0]))
      print('{:.2f}/{:.2f}'.format(t,tEnd))
      outputCount += 1


  return Q


def BC(q):
  rho = q[0:N]
  Momx = q[N:2*N]
  Momy = q[2*N:3*N]
  Pixx = q[3*N:4*N]
  Pixy = q[4*N:5*N]
  Piyx = q[5*N:6*N]
  Piyy = q[6*N:7*N]

  #rho[0]    = rho[1]
  #rho[-1]   = rho[-2]
  rho[:,0]  = rho[:,1]
  rho[:,-1] = rho[:,-2]


  #Momx[0]    = -Momx[1]
  #Momx[-1]   = -Momx[-2]
  Momx[:,0]  = Momx[:,1]
  Momx[:,-1] = Momx[:,-1]

  #Momy[0]    = 0
  #Momy[-1]   = 0
  Momy[:,0]  = -Momy[:,1]
  Momy[:,-1] = -Momy[:,-2]

  #Pixx[0]    = 0
  #Pixx[-1]   = 0
  Pixx[:,0]  = Pixx[:,1]
  Pixx[:,-1] = Pixx[:,-2]

  #Pixy[0]    = 0
  #Pixy[-1]   = 0
  Pixy[:,0]  = Pixy[:,1]
  Pixy[:,-1] = Piyx[:,-2]

  #Piyx[0]    = 0
  #Piyx[-1]   = 0
  Piyx[:,0]  = Piyx[:,1]
  Piyx[:,-1] = Piyx[:,-2]

  #Piyy[0]    = 0
  #Piyy[-1]   = 0
  Piyy[:,0]  = Piyy[:,1]
  Piyy[:,-1] = Piyy[:,-2]
  

t                      = 0    # s 
tEnd                   = 1    # time at the end
tOut                   = 0.01 # time of each output
N                      = 2**9 # resolution
boxsize                = 1.   # in some unit system l
gamma                  = 1.4  # adiabatic index
P0                     = 1    # pressure constant
zeta                   = 2    # bulk viscosity coefficient
eta                    = 2    # shear viscosity coefficient
tau_nu                 = 1    # relaxation time
theta                  = 1    # flux limiter parameter


# Define Mesh
dx = boxsize / N   # box size
dy = dx
vol = dx**2        # volume of each box
a = dx*0.5
b = (boxsize - 0.5*dx)
xlin = np.linspace(a, b, N)# simulation limits

parameters = [t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b]

Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
s = X.shape
R = np.sqrt(X**2 + Y**2)
R1 = np.sqrt((X+np.ones(s))**2 + Y**2)
Theta = np.arctan(Y/X)*(X>=0)*(Y>=0) + (np.pi/2 + np.arctan(Y/np.abs(X)))*(X<0)*(Y>=0) + (np.pi + np.arctan(Y/X))*(X<=0)*(Y<0) + (3*np.pi/2 + np.arctan(np.abs(Y)/X))*(X>0)*(Y<0)

''' initial condition of density'''

#rho = (2*(R1 <= 1) + 0.5*(R1 > 1))
#rho = ((1 - ((R)**2) )**4 )*(R < 1) + 1*np.ones(s) # Mauricio`s funtion advice 
#rho = (1/(R))*(R>0)*(R<1) + 0.1*np.ones(s)
#rho = 1*(Y >= 0) + 0.5*( Y < 0)

''' initial condition of velocity '''
#vx = np.zeros(s)
#vx = -1*np.sin(Theta)*(R < 1)
#vx = 2*np.ones(s)*(R1 < 1)
#vx = -2*(Y < 0) + 2*(Y >= 0)
#vx = np.abs((xlin - (boxsize-0.5*dx)*0.5)/16)

#vy = np.zeros(s)
#vy = 1*np.cos(Theta)*(R < 1)
#vy = 0.5*np.ones(xlin.shape)

#sigma = 0.05/np.sqrt(2.)
#vy = 0.01*np.sin(4*np.pi*X)*(np.exp(-(Y)**2/(2 * sigma**2)) + np.exp(-(Y)**2/(2*sigma**2)))

''' Kelvin Helmholtz'''
w0 = 0.1
sigma = 0.05/np.sqrt(2.)
rho = 1. + (np.abs(Y-0.5) < 0.25)
vx = -0.5 + (np.abs(Y-0.5) < 0.25)
vy = w0*np.sin(4*np.pi*X) * ( np.exp(-(Y-0.25)**2/(2 * sigma**2)) + np.exp(-(Y-0.75)**2/(2*sigma**2)) )

'''
w0 = 0.1
sigma = 0.05/np.sqrt(2.)
rho = 1. + (np.abs(Y-0.5) < 0.25)
vx = -0.5 + (np.abs(Y-0.5) < 0.25)
vy = w0*np.sin(4*np.pi*X) * ( np.exp(-(Y-0.25)**2/(2 * sigma**2)) + np.exp(-(Y-0.75)**2/(2*sigma**2)) )
'''

''' initial condition of Pi tensor '''
Pixx = np.zeros(s)
Pixy = np.zeros(s)
Piyx = np.zeros(s)
Piyy = np.zeros(s)


IC = np.vstack((rho,rho*vx,rho*vy,Pixx,Pixy,Piyx,Piyy)) # here the initial conditions are stacked in a vector 
                            # rho is IC[0:N] for example

sol = np.load("NonRelativisticISKelvinHelmholtz.npy")
t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt("NonRelativisticISKelvinHelmholtz_parameters",delimiter=',',unpack=True)
t=float(t) 
tEnd = float(t) + 2
tOut = float(tOut)
N = int(N)
boxsize = float(boxsize)
gamma = float(gamma)
zeta = float(zeta)
eta = float(eta)
tau_nu = float(tau_nu)
theta = float(theta)
a = float(a)
b = float(b)
IC = sol[0]

# input (dx, dy, xlin, gamma, zeta, tau_nu, BC, theta=1)
# output solution list of arrays that are 7N x N in the order (rho,rho*vx,rho*vy,Pixx,Pixy,Piyx,Piyy)
solution = integrator(KTschemeNonRelativisticIS, (t, tEnd), IC, 0.01, BC, method="Modified_RK", args=(dx, dy, N, gamma, zeta, tau_nu, eta, P0, theta))

np.savetxt('NonRelativisticISKH1_parameters',parameters)
np.save('NonRelativisticISKH1',solution)
