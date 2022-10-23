
from xml.etree.ElementTree import PI
import numpy as np

def minmod2(x,y):
  return (sign(x) + sign(y))*np.minimum(np.abs(x), np.abs(y))
  

def minmod3(x,y,z):
  return minmod2(x,minmod2(y,z))

def sign(x):
  return np.sign(x)

def local_propagation_speed(rho, vx, Pi, cs): 
   '''
    rho          is a matrix of left-state  density
    vx_P         is a matrix of x-velocity
    Pi           is a matrix of bulk viscosity 
    cs           is the speed of sound
    '''

   C1 = np.abs(vx)

   C2 = np.abs(vx - np.sqrt( cs^2 + np.divide(Pi , rho , out=np.zeros_like(Pi), where=vx!=0 )  ) )

   C3 = np.abs(vx + np.sqrt( cs^2 + np.divide(Pi , rho , out=np.zeros_like(Pi), where=vx!=0 )  ) )

   return np.maximum(C1,C2,C3)

def getGradient(f, dx, theta=1):
	"""
    Calculate the gradients of a field
	f        is a matrix of the field
	dx       is the cell size
	f_dx     is a matrix of derivative of f in the x-direction
	
	"""
	# directions for np.roll() 
	R = -1   # right
	L = 1    # left
	
	df_dx = minmod3( theta * ( f - np.roll(f,L,axis=0) )/dx, theta * (np.roll(f,R,axis=0) - np.roll(f,L,axis=0) ) / (2*dx), theta * ( f - np.roll(f,L,axis=0) ) / dx)

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
  
  qP_XL = q - q_dx * dx/2
  qP_XR = np.roll(q,R,axis=0) - np.roll(q_dx,L,axis=0) * dx/2
  qM_XR = q + q_dx * dx/2
  qM_XL = np.roll(q,L,axis=0) + np.roll(q_dx,L,axis=0) * dx/2
  
  return qM_XL, qP_XL, qM_XR, qP_XR

def getFlux(rho_P, rho_M, vx_P, vx_M, Pi_P, Pi_M, P_P, P_M, gamma, cs):
  """
  Calculate fluxed between 2 states with local Kurganov Tadmor rule 
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
  
  # left and right energies
  en_P = P_P/(gamma-1)+0.5*rho_P * (vx_P**2)
  en_M = P_M/(gamma-1)+0.5*rho_M * (vx_M**2)

  # compute star (averaged) states
  rho_star  = 0.5*(rho_P + rho_M)
  momx_star = 0.5*(rho_P * vx_P + rho_M * vx_M)
  en_star   = 0.5*(en_P + en_M)
  Pi_star   = 0.5*(Pi_P + Pi_M)
  
  # Pressure equation of state
  P_star = (gamma-1)*(en_star-np.divide(0.5*(momx_star**2),rho_star, out=np.zeros_like(0.5*(momx_star**2)), where=rho_star!=0)) 
  
  # compute fluxes (local Kurganov-Tadmor)
  flux_Mass   = momx_star
  flux_Momx   = np.divide(momx_star**2,rho_star, out=np.zeros_like(momx_star**2), where=rho_star!=0) + P_star + Pi_star
  flux_Pi_v   = Pi_star * np.divide(momx_star,rho_star, out=np.zeros_like(momx_star), where=rho_star!=0)
  
  # find wavespeeds

  C_P = local_propagation_speed(rho_P , vx_P, Pi_P, cs) # max propagation speed from the left

  C_M = local_propagation_speed(rho_P , vx_P, Pi_P, cs) # max propagation speed from the right

  C = np.maximum(C_M, C_P)

  # add stabilizing diffusive term
  flux_Mass   -= C * 0.5 * (rho_P - rho_M)
  flux_Momx   -= C * 0.5 * (rho_P * vx_P - rho_M * vx_M)
  flux_Pi_v   -= C * 0.5 * ( Pi_P - Pi_M )

  return flux_Mass, flux_Momx, flux_Pi_v



def applyFluxes(H, flux_H1_X, flux_H2_X , dx, J = 0):
    """
    Apply fluxes to conserved variables
    H         is a matrix of the conserved variable field
    flux_H1_X is a matrix of the x-dir fluxes
    flux_H2_X is a matrix of the x-dir fluxes
    dx        is the cell size
    """
    
    # update solution
    H += - (flux_H1_X - flux_H2_X ) / dx
    H += J
    
    return H

