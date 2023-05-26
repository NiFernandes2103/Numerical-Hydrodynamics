import numpy as np
from EoS import *

'''
These are auxiliary functions in for the gradient 
'''

def minmod2(x,y):
  return (np.sign(x) + np.sign(y))*np.minimum(np.abs(x), np.abs(y))/2

def minmod3(x,y,z):
  return minmod2(x,minmod2(y,z))


def local_propagation_speed(rho, vx, Pi, gamma, B): 
  
   '''
    Get the local propagation speeds using the eigenvalues 
    of the flux matrix of the non relativistic IS equations

    rho          is a matrix of left-state  density
    vx           is a matrix of x-velocity
    Pi           is a matrix of bulk viscosity 
    cs           is the speed of sound
    '''

  
   C1 = np.abs(vx)
   cs = getSpeedOfSound(rho, gamma)
   try:
     C2 = np.abs(vx - np.sqrt((cs**2 +(Pi+B)/ rho)/gamma))
   except:
     C2 = np.zeros(vx.shape)
   try:
     C3 = np.abs(vx + np.sqrt((cs**2 + (Pi+B)/rho)/gamma))
   except:
     C3 = np.zeros(vx.shape)

   return np.maximum(C1,C2,C3)


def getGradient(f, dx, theta=1):
    """
    Calculate the gradients of a field
    f        is a matrix of the field
    dx       is the cell size
    f_dx     is a matrix of derivative of f in the x-direction
    theta    is the flux limiter 1 <= theta <= 2
    
    """
    # directions for np.roll() 
    R = -1   # right
    L = 1    # left
    
    df_dx = minmod3( theta * ( f - np.roll(f,L,axis=0) )/dx, (np.roll(f,R,axis=0) - np.roll(f,L,axis=0) ) / (2*dx),theta * ( np.roll(f,R,axis=0) - f ) / dx)

    return df_dx

def extrapolateInSpaceToFace(q, q_dx, dx):
  
  """
  Calculate the gradients of a field
  q        is a matrix of the field
  q_dx     is a matrix of the field x-derivatives
  dx       is the cell size
  q_XL     is a matrix of spatial-extrapolated values on `left' face along x-axis 
  q_XR     is a matrix of spatial-extrapolated values on `right' face along x-axis 
  """
  # directions for np.roll() 
  R = -1   # right
  L = 1    # left
  
  qP_XL = q - q_dx * dx/2 # minus
  qP_XR = np.roll(q,R,axis=0) - np.roll(q_dx,R,axis=0) * dx/2
  qM_XR = q + q_dx * dx/2
  qM_XL = np.roll(q,L,axis=0) + np.roll(q_dx,L,axis=0) * dx/2
  
  return qM_XL, qP_XL, qM_XR, qP_XR

def getFlux(rho_P, rho_M, vx_P, vx_M, Pi_P, Pi_M, P_P, P_M, gamma, B):

  """

  Calculate fluxes between 2 states with local Kurganov Tadmor rule 
  rho_P        is a matrix of left-state  density
  rho_M        is a matrix of right-state density
  vx_P         is a matrix of left-state  x-velocity
  vx_M         is a matrix of right-state x-velocity
  Pi_P         is a matrix of left-state bulk viscosity 
  Pi_M         is a matrix of right-state bulk viscosity
  P_P          is a matrix of left-state  pressure
  P_M          is a matrix of right-state pressure
  gamma        is the ideal gas gamma
  flux_Mass    is the matrix of mass fluxes
  flux_Momx    is the matrix of x-momentum fluxes
  flux_Pi_v    is the matrix of the bulk viscosity var

  """
  

  # compute (averaged) states over the left and right states
  momx_av  = 0.5*(rho_P * vx_P + rho_M * vx_M)
  Pi_av    = 0.5*(Pi_P + Pi_M)
  Pi_vx_av = 0.5*(Pi_P * vx_P + Pi_M * vx_M)
  P_av     = 0.5*(P_P + P_M)
  
  # compute fluxes (local Kurganov-Tadmor)

  flux_Mass   = momx_av
  flux_Momx   = 0.5*(rho_P*(vx_P)**2 + rho_M*(vx_M)**2) + (P_av + Pi_av)/gamma
  flux_Pi_v   = Pi_vx_av + B * (vx_P + vx_M)*0.5 
  
  # find wavespeeds

  C_P = local_propagation_speed(rho_P , vx_P, Pi_P, gamma, B) # max propagation speed from the left

  C_M = local_propagation_speed(rho_M , vx_M, Pi_M, gamma, B) # max propagation speed from the right

  C = np.maximum(C_M, C_P)

  # add stabilizing diffusive term
  flux_Mass   -= C * 0.5 * (rho_P - rho_M)
  flux_Momx   -= C * 0.5 * (rho_P * vx_P - rho_M * vx_M)
  flux_Pi_v   -= C * 0.5 * ( Pi_P - Pi_M )

  return flux_Mass, flux_Momx, flux_Pi_v

def applyFluxes(flux_H1_X, flux_H2_X , dx, J = 0):
    """
    Apply fluxes to conserved variables
    H         is a matrix of the conserved variable field
    flux_H1_X is a matrix of the x-dir fluxes from the right 
    flux_H2_X is a matrix of the x-dir fluxes from the left
    dx        is the cell size
    """
    C = 0

    # update solution
    C -= (flux_H1_X - flux_H2_X ) / dx
    C += J
    
    return C

def KTschemeNonRelativisticIS(t,IC, dx, xlin, gamma, zeta, tau_nu, BC = None, theta=1):

    """ Finite Volume simulation """

    if BC is not None:
      BC(IC)
  
    #-----------------------------------------------------------------------------------------------------------------------------------#

    # Generate Initial Conditions  

    N = xlin.shape[0]

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

    # get Conserved variables
    vol = dx*dx
    Mass, Momx, Pi_vx = getConserved( rho, vx, Pi, vol)


    # get Primitive variables
    rho, vx, P = getPrimitive( Mass, Momx, gamma, vol )


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
    # input getFlux(rho_P, rho_M, vx_P, vx_M, Pi_P, Pi_M, P_P, P_M, gamma, B)
    # output flux_Mass, flux_Momx, flux_Pi_v

    flux_Mass_XR, flux_Momx_XR, flux_Pi_vxR = getFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR, PiP_XR, PiM_XR, PP_XR, PM_XR, gamma, B)
    flux_Mass_XL, flux_Momx_XL, flux_Pi_vxL = getFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL, PiP_XL, PiM_XL, PP_XL, PM_XL, gamma, B)


    #-----------------------------------------------------------------------------------------------------------------------------------#
    
    # get time derivative 

    J = - Pi/tau_nu
    timederivative_rho = applyFluxes( flux_Mass_XR,   flux_Mass_XL,   dx)
    timederivative_Momx = applyFluxes( flux_Momx_XR,   flux_Momx_XL,   dx) 
    timederivative_Pi  = applyFluxes( flux_Pi_vxR,    flux_Pi_vxL,    dx, J)

    #-----------------------------------------------------------------------------------------------------------------------------------#
    


    return np.hstack((timederivative_rho,timederivative_Momx,timederivative_Pi))