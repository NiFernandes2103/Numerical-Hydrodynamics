import numpy as np
import matplotlib.pyplot as plt

def Lorentz_factor(vx):
  return np.sqrt(1/(1-np.abs(vx)**2))


def polytrope(rho, gamma, P0):
   return P0*pow(rho,gamma)

def enthalpy_density(eps,P,rho):
   s = rho.shape
   return np.ones(s) + eps + P + rho

def f(rho,gamma,P0,P,Eos = polytrope):
   return Eos(rho,gamma,P0) - P

def fprime(vx, cs):
   s = vx.shape
   return vx*vx*cs*cs - np.ones(s) 

def getConserved( rho, vx, eps, gamma, P0, Eos = polytrope):
    """
    Calculate the conserved variable from the primitive
    rho      is the matrix of cell densities
    vx       is the matrix of cell x-velocity
    gamma    is the ideal gas gamma
    D        is the matrix of rest-mass density in cells
    Sx       is the matrix of x-momentum density in cells
    tau      is the energy density
    """

    P  = Eos(rho,gamma,P0)
    W  = Lorentz_factor(vx)
    W2 = W**2
    h  = enthalpy_density(eps,P,rho)
    D    = rho * W
    Sx   = rho * h * W2 * vx
    tau  = rho * h * W2 - P - D
    
    return (D, Sx, tau)

def getPrimitive( D, Sx, tau, gamma, P0, P, n, Eos=polytrope):
  
  """
  Calculate the primitive variable from the conservative
  Mass     is matrix of mass in cells
  Sx     is matrix of x-momentum in cells
  gamma    is ideal gas gamma
  rho      is the matrix of cell densities
  vx       is the matrix of cell x-velocity
  P        is the matrix of cell pressures
  """

  Pnew = P
  for i in range(n):
       vx = Sx/(tau + D + Pnew)
       W  = Lorentz_factor(vx)
       rho = D / np.where(vx == 0, 1, vx)       
       eps = (tau + D*(1-W) + Pnew*(1-W**2))/(D*W) 
       cs = getSpeedOfSound(rho,Pnew,gamma)
       Pnew = Pnew - f(rho,gamma,P0,P,Eos) / fprime(vx,cs)

  return (rho, vx, eps)

def getSpeedOfSound(rho, P, gamma):
  '''
  find the speed of sound in the fluid
  rho

  '''

  cs = np.sqrt((gamma)*P/rho)

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

    
    df_dx =  minmod3( theta * ( f - f[Km1] )/dx, (f[Kp1] - f[Km1] ) / (2*dx),theta * ( f[Kp1] - f ) / dx)
    
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

    n = q.shape[axis]

    K = np.arange(0, n)
    Kp1 = np.roll(K, -1)
    Km1 = np.roll(K, 1)


    qP_XL = np.zeros_like(q)
    qP_XR = np.zeros_like(q)
    qM_XR = np.zeros_like(q)
    qM_XL = np.zeros_like(q)
       
    qP_XL = q - q_dx * dx/2
    qP_XR = q[Kp1] - q_dx[Kp1] * dx/2
    qM_XR = q + q_dx * dx/2
    qM_XL = q[Km1] + q_dx[Km1] * dx/2

    
    return (qM_XL, qP_XL, qM_XR, qP_XR)


def local_propagation_speed(vx, cs): 
  
   '''
    Get the local propagation speeds using the eigenvalues 
    of the flux matrix of the non relativistic IS equations

    rho          is a matrix of density

    cs           is the speed of sound
    '''

   C = (1/(1-(vx*cs)**2))*(vx*(1-cs**2) + cs*(1-vx**2))

   return C

def getXFlux(D_P, D_M, Sx_P, Sx_M, tau_P, tau_M, gamma, P0, Eos=polytrope):

    """

    Calculate fluxes between 2 states with local Kurganov Tadmor rule 
    rho_P        is a matrix of left-state  density
    rho_M        is a matrix of right-state density
    vx_P         is a matrix of left-state  x-velocity
    vx_M         is a matrix of right-state x-velocity
    P_P          is a matrix of left-state  pressure
    P_M          is a matrix of right-state pressure
    gamma        is the ideal gas gamma
    flux_D       is the matrix of rest-mass density fluxes
    flux_Sx      is the matrix of x-momentum density fluxes
    flux_tau     is the matrix of the energy density fluxes

    """

    # compute conserved quatities

    rho_P,vx_P,eps_P = getPrimitive(D_P,Sx_P,tau_P, gamma, P0, 0.5*np.ones(D_P.shape), 100, Eos)
    rho_M,vx_M,eps_M = getPrimitive(D_M,Sx_M,tau_M, gamma, P0, 0.5*np.ones(D_P.shape), 100, Eos)

    P_P = Eos(rho_P,gamma,P0)
    P_M = Eos(rho_M,gamma,P0)

    # compute (averaged) states over the left and right states

    # x fluxes
    
    D_v_av   = 0.5*(D_P * vx_P + D_M * vx_M)
    Sx_av    = 0.5*(Sx_P + Sx_M)
    P_av     = 0.5*(P_P + P_M)    
    
    # compute fluxes 


    flux_D    = D_v_av  
    flux_Sx   = 0.5*(Sx_P*vx_P + Sx_M*vx_M) + P_av
    flux_tau  = Sx_av - D_v_av
    
    # find wavespeeds (use c =1)
    cs_P = getSpeedOfSound(rho_P,P_P,gamma)
    C_P = local_propagation_speed(vx_P,cs_P)

    cs_M = getSpeedOfSound(rho_M,P_M,gamma)
    C_M = local_propagation_speed(vx_M,cs_M)

    C = np.maximum(C_P,C_M)

    # add stabilizing diffusive term
    flux_D       -= C * 0.5 * (D_P - D_M)
    flux_Sx      -= C * 0.5 * (Sx_P - Sx_M)
    flux_tau     -= C * 0.5 * (tau_P - tau_M)

    return (flux_D, flux_Sx, flux_tau)


def applyFluxes(flux_H1_X, flux_H2_X, dx, J = 0):
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
    C += J
    
    return C

def Heuns(q,f,dt,t):

  k1 = dt*f(t,q)
  k2 = dt*f(t + dt,q + k1)

  return q + 0.5 * (k1 + k2)

def Euler_step(q,f,t,dt):
   return q + dt*f(t,q)

def explicit_modified_RK(q,f,dt,t):
   q1 = Euler_step(q,f,t,dt)
   q2 = 1/2 * q1 + 1/2 * (Euler_step(q1,f,t,dt))

   return q2

def RK4(y0,f,h,t):
  
  k1 = h * (f(t, y0))
  k2 = h * (f((t+h/2), (y0+k1/2)))
  k3 = h * (f((t+h/2), (y0+k2/2)))
  k4 = h * (f((t+h), (y0+k3)))
  k = (k1+2*k2+2*k3+k4)/6
  yn = y0 + k

  return yn
