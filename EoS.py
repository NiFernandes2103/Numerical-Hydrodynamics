import numpy as np


def getConserved( rho, vx, P, gamma, vol ):
	"""
    Calculate the conserved variable from the primitive
	rho      is matrix of cell densities
	vx       is matrix of cell x-velocity
	P        is matrix of cell pressures
	gamma    is ideal gas gamma
	vol      is cell volume
	Mass     is matrix of mass in cells
	Momx     is matrix of r-momentum in cells
	Energy   is matrix of energy in cells
	"""
	Mass   = rho * vol
	Momr   = rho * vx * vol
	Energy = (P/(gamma-1) + 0.5* rho *(vx**2)*vol)
	
	return Mass, Momr, Energy

def getPrimitive( Mass, Momx, Energy, gamma, vol ):
  """
  Calculate the primitive variable from the conservative
  Mass     is matrix of mass in cells
  Momx     is matrix of x-momentum in cells
  Energy   is matrix of energy in cells
  gamma    is ideal gas gamma
  vol      is cell volume
  rho      is matrix of cell densities
  vx       is matrix of cell x-velocity
  P        is matrix of cell pressures
  """
  rho = Mass / vol
  vx  = np.divide(Momx , rho, out=np.zeros_like(Momx), where=rho!=0) / vol
  P   = (Energy/vol - 0.5*rho * (vx)) * (gamma-1)
  
  return rho, vx, P