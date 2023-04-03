import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from KTalgorithm import *
from EoS import *
from RK_Heuns_integrator import *

def KTschemeNonRelativisticIS(t,IC, dx, xlin, gamma, zeta, tau_nu, BC, theta=1):

    """ Finite Volume simulation """
  
    #-----------------------------------------------------------------------------------------------------------------------------------#

    # Generate Initial Conditions  

    ''' Initial conditions for rho '''
    rho = IC[0:N]

    
    ''' Initial conditions for v'''
    vx = np.divide(IC[N:2*N] , rho, out=np.zeros_like(IC[N:2*N]), where=rho!=0)

    ''' Pressure due to equation of state '''
    P = (np.abs(rho))**gamma

    ''' B Constant '''
    B = zeta/tau_nu

    ''' Pi initial condition '''
    Pi = IC[2*N:]


   #-----------------------------------------------------------------------------------------------------------------------------------#


    # getSpeedOfSound(rho, gamma)
    cs = getSpeedOfSound(rho,gamma)
    

    #-----------------------------------------------------------------------------------------------------------------------------------#

    # get Conserved variables
    vol = dx*dx
    Mass, Momx, Pi_vx = getConserved( rho, vx, Pi, gamma, vol)


    # get Primitive variables
    rho, vx, P = getPrimitive( Mass, Momx, gamma, vol )

    # get Speed of sound
    cs = getSpeedOfSound(rho, gamma)

    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # calculate gradients
    # getGradient(f, dx, theta=1)

    rho_dx = getGradient(rho,  dx,theta)
    vx_dx  = getGradient(vx,   dx,theta)
    P_dx   = getGradient(P,    dx,theta)
    Pi_dx  = getGradient(Pi,   dx,theta)

    
    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # extrapolate in space to face centers
    # input extrapolateInSpaceToFace(q, q_dx, dx)
    # output qM_XL, qP_XL, qM_XR, qP_XR

    rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR = extrapolateInSpaceToFace(rho, rho_dx, dx)
    vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR  = extrapolateInSpaceToFace(vx,  vx_dx,  dx)
    PM_XL,   PP_XL,   PM_XR,   PP_XR   = extrapolateInSpaceToFace(P,   P_dx,   dx)
    PiM_XL,  PiP_XL,  PiM_XR,  PiP_XR  = extrapolateInSpaceToFace(Pi,  Pi_dx,  dx)


    #-----------------------------------------------------------------------------------------------------------------------------------# 
    
    # compute fluxes (local Kurganov-Tadmor)
    # getFlux(rho_P, rho_M, vx_P, vx_M, Pi_P, Pi_M, P_P, P_M, gamma, B)
    # output flux_Mass, flux_Momx, flux_Pi_v

    flux_Mass_XR, flux_Momx_XR, flux_Pi_vxR = getFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR, PiP_XR, PiM_XR, PP_XR, PM_XR, gamma, B)
    flux_Mass_XL, flux_Momx_XL, flux_Pi_vxL = getFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL, PiP_XL, PiM_XL, PP_XL, PM_XL, gamma, B)

    

    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # update solution

    J = - Pi
    timederivative_rho = applyFluxes( flux_Mass_XR,   flux_Mass_XL,   dx)
    timederivative_Momx = applyFluxes( flux_Momx_XR,   flux_Momx_XL,   dx) 
    timederivative_Pi  = applyFluxes( flux_Pi_vxR,    flux_Pi_vxL,    dx, J)

    #-----------------------------------------------------------------------------------------------------------------------------------#
    


    return np.hstack((timederivative_rho,timederivative_Momx,timederivative_Pi))
  
def integrator(scheme, time, y0, dt, BC, methodName = "Heuns", *args):

  t, tEnd = time

  y = y0 
  
  while t < tEnd: 

    C = scheme(t,y0,args)

    if methodName == "Heuns":
      y = Heuns(y,C,dt,t)
    if methodName == "RK4":
      y = RK4(y,C,dt,t)
    
    #Apply Boundary conditions

    BC(y)

    t = t+dt
    
  return y

def applyBC(y):

  rho = y[0:N] 
  vx = y[N:2N]/rho 
  Pi = y[2N:] 

  '''Absorbing boundary conditions'''
  rho[0]    = rho[1]
  rho[-1]   = rho[-2]

  vx[0]     = vx[1]      
  vx[-1]    =  vx[-2]

  Pi[0]     = 0    
  Pi[-1]    = 0



t                      = 0   # s 
tEnd                   = 2   # time at the end
tOut                   = 0.01 # time of each output

N                      = 400 # resolution
boxsize                = 1.  # in some unit system l
gamma                  = 2 # adiabatic index
zeta                   = 1 # bulk viscosity coefficient
tau_nu                 = 1
theta                  = 1


# Define Mesh
dx = boxsize / N   # box size
vol = dx**2        # volume of each box
xlin = np.linspace(0.5*dx, boxsize-0.5*dx, N) # simulation limits


#rho = ((1 - ((xlin - (boxsize-0.5*dx)*0.5)**2)/0.25 )**4 ) + 0.5*np.ones(xlin.shape) # Mauricio`s funtion advice    
rho = 1*(xlin < boxsize*0.5) + 0.125*(xlin >= boxsize*0.5)

#vx = np.zeros(xlin.shape)
vx = 0.5*np.ones(xlin.shape)
#vx = np.abs((xlin - (boxsize-0.5*dx)*0.5)/16)


Pi = np.zeros(xlin.shape)

IC = np.hstack((rho,rho*vx,Pi)) # here the initial conditions are stacked in a vector 
                            # rho is IC[0:N] for example




                      

solution = integrate.solve_ivp(KTschemeNonRelativisticIS, (t,tEnd), IC, args=(dx, xlin, gamma, zeta, tau_nu, applyBC, theta), vectorized=True)

print(solution.y.shape)
plt.plot(xlin,solution.y[0:N,0])
plt.show()