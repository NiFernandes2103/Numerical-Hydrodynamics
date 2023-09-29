import numpy as np


def getConserved( rho, vx, Pi, vol ):
    """
    Calculate the conserved variable from the primitive
    rho      is the matrix of cell densities
    vx       is the matrix of cell x-velocity
    gamma    is the ideal gas gamma
    vol      is the cell volume
    Mass     is the matrix of mass in cells
    Momx     is the matrix of x-momentum in cells
    Pi_vx    is the matrix of Pi * vx fluxes
    """
    Mass   = rho * vol
    Momx   = rho * vx 
    Pi_vx  = Pi  * vx

    return Mass, Momx, Pi_vx
  
def getPrimitive( Mass, Momx, gamma, vol):
  """
  Calculate the primitive variable from the conservative
  Mass     is the matrix of mass in cells
  Momx     is the matrix of x-momentum in cells
  gamma    is the ideal gas gamma
  vol      is the cell volume
  rho      is the matrix of cell densities
  vx       is the matrix of cell x-velocity
  P        is the matrix of cell pressures
  """
  rho = Mass / vol
  vx  = np.divide(Momx , rho, out=np.zeros_like(Momx), where=rho!=0)
  P   = (np.abs(rho))**gamma
  
  return rho, vx, P

def getSpeedOfSound(rho, gamma):
  """
  Calculate the speed of sound based on
  equation of state
  """
  cs = np.sqrt((gamma)*rho**(gamma-1))

  return cs

  