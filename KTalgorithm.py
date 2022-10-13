import numpy as np

def minmod2(x,y):
  return (sign(x) + sign(y))*np.min(abs(x)*abs(y))

def minmod3(x,y,z):
  return minmod2(x,minmod2(y,z))

def sign(x):
  return np.sign(x)

def q_r(q,i,dr,theta):
  return minmod3(theta*(q(i) - q(i-1))/dr,theta*(q(i+1) - q(i-1))/dr,theta*(q(i) - q(i-1))/(2*dr),theta*(q(i+1) - q(i))/dr)

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
	
	df_dx = minmod3( theta * ( f - np.roll(f,L,axis=0) )/dx, theta * (np.roll(f,R,axis=0) - np.roll(f,L,axis=0) ) / (2*dx),theta * ( f - np.roll(f,L,axis=0) ) / dx)

	return df_dx

def getTimederivative(f, dtau):

	# directions for np.roll() 
	R = -1   # right
	L = 1    # left
	
	df_dtau = ( np.roll(f,R,axis=0) - np.roll(f,L,axis=0) ) / (dtau)

	return df_dtau

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
  flux_Energy  is the matrix of energy fluxes
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
  
  # compute fluxes (local Kurganov-Tadamor)
  flux_Mass   = momx_star
  flux_Momx   = np.divide(momx_star**2,rho_star, out=np.zeros_like(momx_star**2), where=rho_star!=0) + P_star + P_star
  flux_Energy = (en_star+P_star) * np.divide(momx_star,rho_star, out=np.zeros_like(momx_star), where=rho_star!=0)
  flux_Pi_v   = Pi_star * np.divide(momx_star,rho_star, out=np.zeros_like(momx_star), where=rho_star!=0)
  
  # find wavespeeds

  C1_P = vx_P
  C1_M = vx_M
  C1 = np.maximum( C1_P, C1_M )

  C2_P = rho_P*vx_P - np.sqrt(np.divide(Pi_P * rho_P + cs**2 * rho_P**2 - rho_P * vx_P**2 +  rho_P**2 * vx_P**2 ,rho_P,
   out=np.zeros_like(Pi_P * rho_P + cs**2 * rho_P**2 - rho_P * vx_P**2 +  rho_P**2 * vx_P**2), where=rho_P!=0))
  C2_M = rho_M*vx_M - np.sqrt(np.divide(Pi_M * rho_M + cs**2 * rho_M**2 - rho_M * vx_M**2 +  rho_M**2 * vx_M**2 ,rho_M,
   out=np.zeros_like(Pi_M * rho_M + cs**2 * rho_M**2 - rho_M * vx_M**2 +  rho_M**2 * vx_M**2), where=rho_M!=0))
  C2 = np.maximum( C2_P, C2_M )

  C3_P = rho_P*vx_P + np.sqrt(np.divide(Pi_P * rho_P + cs**2 * rho_P**2 - rho_P * vx_P**2 +  rho_P**2 * vx_P**2 ,rho_P,
   out=np.zeros_like(Pi_P * rho_P + cs**2 * rho_P**2 - rho_P * vx_P**2 +  rho_P**2 * vx_P**2), where=rho_P!=0))
  C3_M = rho_M*vx_M - np.sqrt(np.divide(Pi_M * rho_M + cs**2 * rho_M**2 - rho_M * vx_M**2 +  rho_M**2 * vx_M**2 ,rho_M,
   out=np.zeros_like(Pi_M * rho_M + cs**2 * rho_M**2 - rho_M * vx_M**2 +  rho_M**2 * vx_M**2), where=rho_M!=0))
  C3 = np.maximum( C3_P, C3_M )

  C = np.maximum( C1, C2, C3 )
  
  # add stabilizing diffusive term
  flux_Mass   -= C * 0.5 * (rho_P - rho_M)
  flux_Momx   -= C * 0.5 * (rho_P * vx_P - rho_M * vx_M)
  flux_Energy -= C * 0.5 * ( en_P - en_M )
  flux_Pi_v   -= C * 0.5 * ( Pi_P - Pi_M )

  return flux_Mass, flux_Momx, flux_Energy, flux_Pi_v

def applyFluxes(H, flux_H_X, dx, J = 0):
	"""
    Apply fluxes to conserved variables
	H        is a matrix of the conserved variable field
	flux_H_X is a matrix of the x-dir fluxes
	dx       is the cell size
	"""
	# directions for np.roll() 
	R = -1   # right
	L = 1    # left
	
	# update solution
	H += - flux_H_X / dx
	H += J
	
	return H