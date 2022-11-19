import numpy as np


def getConserved( rho, vx, vol ):
    """
    Calculate the conserved variable from the primitive
    rho      is matrix of cell densities
    vx       is matrix of cell x-velocity
    gamma    is ideal gas gamma
    vol      is cell volume
    Mass     is matrix of mass in cells
    Momx     is matrix of x-momentum in cells
    """
    Mass   = rho * vol
    Momr   = rho * vx * vol

    return Mass, Momr
  
def getPrimitive( Mass, Momx, gamma, vol):
  """
  Calculate the primitive variable from the conservative
  Mass     is matrix of mass in cells
  Momx     is matrix of x-momentum in cells
  gamma    is ideal gas gamma
  vol      is cell volume
  rho      is matrix of cell densities
  vx       is matrix of cell r-velocity
  P        is matrix of cell pressures
  """
  rho = Mass / vol
  vx  = np.divide(Momx , rho, out=np.zeros_like(Momx), where=rho!=0)
  P   = (np.abs(rho))**gamma
  
  return rho, vx, P

def getSpeedOfSound(rho, gamma):
  
  cs = np.sqrt((gamma)*np.abs(rho)**(gamma-1))

  return cs


def regParams(Pi,rho_max):

  r = np.amax(np.abs(Pi)/rho_max)

  return np.divide(np.tanh(r) * Pi , r, out=np.zeros_like(np.tanh(r) * Pi), where=r!=0)
  