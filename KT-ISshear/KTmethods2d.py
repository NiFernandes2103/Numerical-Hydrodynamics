import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import integrate


def getConserved( rho, vx, vy, gamma, vol ):

    """
    Calculate the conserved variable from the primitive
    rho      is matrix of cell densities
    vx       is matrix of cell x-velocity
    vy       is matrix of cell y-velocity
    gamma    is ideal gas gamma
    vol      is cell volume
    Mass     is matrix of mass in cells
    Momx     is matrix of x-momentum in cells
    Momy     is matrix of y-momentum in cells
    """

    Mass   = rho * vol
    Momx   = rho * vx
    Momy   = rho * vy  
    
    return Mass, Momx, Momy

def getPrimitive( Mass, Momx, Momy, gamma, vol):
  """
  Calculate the primitive variable from the conservative
  Mass     is matrix of mass in cells
  Momx     is matrix of x-momentum in cells
  Momy     is matrix of y-momentum in cells
  gamma    is ideal gas gamma
  vol      is cell volume
  rho      is matrix of cell densities
  vx       is matrix of cell x-velocity
  vy       is matrix of cell y-velocity
  P        is matrix of cell pressures
  """
  rho = Mass / vol
  vx  = np.divide(Momx , rho, out=np.zeros_like(Momx), where=rho!=0)
  vy  = np.divide(Momy , rho, out=np.zeros_like(Momy), where=rho!=0)
  P   = (np.abs(rho))**gamma
  
  return rho, vx, vy, P

def getSpeedOfSound(rho, gamma):
  '''
  find the speed of sound in the fluid
  rho

  '''
  cs = np.sqrt((gamma)*np.abs(rho)**(gamma-1))

  return cs
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


def local_propagation_speed(rho, eta, zeta, tau_nu, cs): 
  
   '''
    Get the local propagation speeds using the eigenvalues 
    of the flux matrix of the non relativistic IS equations

    rho          is a matrix of density
    eta
    zeta
    tau_nu
    cs           is the speed of sound
    '''

   eps= 10**(-6)

   C1 = np.sqrt(eta*tau_nu / (rho + eps))

   C2 = np.sqrt(cs**2 * (zeta + 4/3 * tau_nu)/ (rho*tau_nu + eps))

   return np.maximum(C1,C2)

def getXFlux(rho_P, rho_M, vx_P, vx_M, vy_P, vy_M, Pixx_P, Pixx_M, Pixy_P,
              Pixy_M, Piyx_P, Piyx_M, Piyy_P, Piyy_M, P_P, P_M, gamma, eta,
              zeta, tau_nu):

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

    # x fluxes
    
    rho_av   = 0.5*(rho_P + rho_M)
    momx_av  = 0.5*(rho_P * vx_P + rho_M * vx_M)
    Pixx_av  = 0.5*(Pixx_P + Pixx_M)    
    Piyx_av  = 0.5*(Piyx_P + Piyx_M)
    Pixx_vx_av = 0.5*(Pixx_P * vx_P + Pixx_M * vx_M)
    Pixy_vx_av = 0.5*(Pixy_P * vx_P + Piyx_M * vx_M)
    Piyx_vx_av = 0.5*(Piyx_P * vx_P + Pixy_M * vx_M)
    Piyy_vx_av = 0.5*(Piyy_P * vx_P + Piyy_M * vx_M)
    P_av     = 0.5*(P_P + P_M)
    
    # compute fluxes 

    B = eta/tau_nu
    A = zeta/tau_nu

    flux_Mass   = momx_av  
    flux_Momx   = 0.5*(rho_P*(vx_P)**2 + rho_M*(vx_M)**2) + P_av + (Pixx_av)/gamma
    flux_Momy   = 0.5*(rho_P*(vx_P*vy_P) + rho_M*(vx_M*vy_M)) + (Piyx_av)/gamma
    flux_Pixx_vx   = Pixx_vx_av + B * (vx_P + vx_M) + (A - 2/3 * B) * (vx_P + vx_M) * 0.5
    flux_Pixy_vx   = Pixy_vx_av + B * (vy_P + vy_M) * 0.5
    flux_Piyx_vx   = Piyx_vx_av + B * (vy_P + vy_M) * 0.5
    flux_Piyy_vx   = Piyy_vx_av + (A - 2/3 * B) * (vx_P + vx_M) * 0.5
    
    # find wavespeeds

    cs_P = getSpeedOfSound(rho_P,gamma)

    C_P = local_propagation_speed(rho_P , eta, zeta, tau_nu, cs_P) # max propagation speed from the left

    cs_M = getSpeedOfSound(rho_M,gamma)

    C_M = local_propagation_speed(rho_M , eta, zeta, tau_nu, cs_M) # max propagation speed from the right

    C = np.maximum(C_M, C_P)

    # add stabilizing diffusive term
    flux_Mass    -= C * 0.5 * (rho_P - rho_M)
    flux_Momx    -= C * 0.5 * (rho_P * vx_P - rho_M * vx_M)
    flux_Momy    -= C * 0.5 * (rho_P * vy_P - rho_M * vy_M)
    flux_Pixx_vx -= C * 0.5 * ( Pixx_P - Pixx_M )
    flux_Pixy_vx -= C * 0.5 * ( Pixy_P - Pixy_M )
    flux_Piyx_vx -= C * 0.5 * ( Piyx_P - Piyx_M )
    flux_Piyy_vx -= C * 0.5 * ( Piyy_P - Piyy_M )

    return flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vx, flux_Pixy_vx, flux_Piyx_vx, flux_Piyy_vx

def getYFlux(rho_P, rho_M, vx_P, vx_M, vy_P, vy_M, Pixx_P, Pixx_M, Pixy_P,
              Pixy_M, Piyx_P,Piyx_M, Piyy_P, Piyy_M, P_P, P_M, gamma, eta,
              zeta, tau_nu):
  

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

  # y fluxes
  
  rho_av   = 0.5*(rho_P + rho_M)
  momy_av  = 0.5*(rho_P * vy_P + rho_M * vy_M)
  Piyy_av  = 0.5*(Piyy_P + Piyy_M)
  Pixy_av  = 0.5*(Pixy_P + Pixy_M)
  Pixx_vy_av = 0.5*(Pixx_P * vy_P + Pixx_M * vy_M)
  Pixy_vy_av = 0.5*(Pixy_P * vy_P + Piyx_M * vy_M)
  Piyx_vy_av = 0.5*(Piyx_P * vy_P + Pixy_M * vy_M)
  Piyy_vy_av = 0.5*(Piyy_P * vy_P + Piyy_M * vy_M)
  P_av     = 0.5*(P_P + P_M)
  
  # compute fluxes 

  B = eta/tau_nu
  A = zeta/tau_nu

  flux_Mass   = momy_av 
  flux_Momx   = 0.5*(rho_P*(vx_P*vy_P) + rho_M*(vx_M*vy_M)) + (Pixy_av)/gamma
  flux_Momy   = 0.5*(rho_P*(vy_P)**2 + rho_M*(vy_M)**2) + (P_av) + (Piyy_av)/gamma
  flux_Pixx_vy   = Pixx_vy_av + (A - 2/3 * B) * (vy_P + vy_M) * 0.5
  flux_Pixy_vy   = Pixy_vy_av + B * (vx_P + vx_M) * 0.5
  flux_Piyx_vy   = Piyx_vy_av + B * (vx_P + vx_M) * 0.5
  flux_Piyy_vy   = Piyy_vy_av + B * (vy_P + vy_M) + (A - 2/3 * B) * (vy_P + vy_M) * 0.5

  
  # find wavespeeds

  cs_P = getSpeedOfSound(rho_P, gamma)

  C_P = local_propagation_speed(rho_P , eta, zeta, tau_nu, cs_P) # max propagation speed from the left

  cs_M = getSpeedOfSound(rho_M, gamma)

  C_M = local_propagation_speed(rho_M , eta, zeta, tau_nu, cs_M) # max propagation speed from the right

  C = np.maximum(C_M, C_P)

  # add stabilizing diffusive term
  flux_Mass    -= C * 0.5 * (rho_P - rho_M)
  flux_Momx    -= C * 0.5 * (rho_P * vx_P - rho_M * vx_M)
  flux_Momy    -= C * 0.5 * (rho_P * vy_P - rho_M * vy_M)
  flux_Pixx_vy -= C * 0.5 * ( Pixx_P - Pixx_M )
  flux_Pixy_vy -= C * 0.5 * ( Pixy_P - Pixy_M )
  flux_Piyx_vy -= C * 0.5 * ( Piyx_P - Piyx_M )
  flux_Piyy_vy -= C * 0.5 * ( Piyy_P - Piyy_M )

  return flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vy, flux_Pixy_vy, flux_Piyx_vy, flux_Piyy_vy


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

def Heuns(q,f,dt,t):

  k1 = dt*f(t,q)
  k2 = dt*f(t + dt,q + k1)

  return q + 1/2 * (k1 + k2)

def RK4(y0,f,h,t):
  
  k1 = h * (f(t, y0))
  k2 = h * (f((t+h/2), (y0+k1/2)))
  k3 = h * (f((t+h/2), (y0+k2/2)))
  k4 = h * (f((t+h), (y0+k3)))
  k = (k1+2*k2+2*k3+k4)/6
  yn = y0 + k

  return yn
