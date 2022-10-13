import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d
from KTalgorithm import *
from EoS import *


def main():
    """ Finite Volume simulation """
    
    # Simulation parameters
    N                      = 128 # resolution
    boxsize                = 1.
    gamma                  = 5/3 # ideal gas gamma
    zeta                   = 1   # bulk viscosity coefficient
    tau_nu                 = 1
    courant_fac            = 0.4
    cs                     = 1/3
    t                      = 0
    tEnd                   = 2
    tOut                   = 0.02
    plotRealTime = True # switch on for plotting as the simulation goes along
    
    # Mesh
    dx = boxsize / N
    vol = dx**2
    xlin = np.linspace(0.5*dx, boxsize-0.5*dx, N)
    Y, X = np.meshgrid( xlin, xlin )
    
    # Generate Initial Conditions 
    w0 = 0.1
    sigma = 0.05/np.sqrt(2.)
    rho = np.random.uniform(0.1,0.5, (X).shape)
    vx = w0*np.sin(4*np.pi*X) * ( np.exp(-(X-0.25)**2/(2 * sigma**2)) + np.exp(-(X-0.75)**2/(2*sigma**2)))
    P = 2.5 * np.ones(X.shape)
    Pi = np.zeros(X.shape)
    Pi_vx = Pi*vx

    # Get conserved variables
    Mass, Momx, Energy = getConserved( rho, vx, P, gamma, vol)

    # prep figure
    fig = plt.figure(figsize=(4,4), dpi=80)
    outputCount = 1

    # Simulation Main Loop
    while t < tEnd:
        
        # get Primitive variables
        rho, vx, P = getPrimitive( Mass, Momx, Energy, gamma, vol )
        
        # get time step (CFL) = dx / max signal speed
        dt = courant_fac * np.min( np.divide(dx , (np.sqrt( np.divide( gamma*P, rho, out=np.zeros_like(gamma*P), where=rho!=0)) + np.sqrt(vx**2)),
		  out=np.zeros_like((np.sqrt( np.divide( gamma*P, rho, out=np.zeros_like(gamma*P), where=rho!=0)) + np.sqrt(vx**2))),
		   where=(np.sqrt( np.divide( gamma*P, rho, out=np.zeros_like(gamma*P), where=rho!=0)) + np.sqrt(vx**2))!=0) )
        plotThisTurn = False
        if t + dt > outputCount*tOut:
            dt = outputCount*tOut - t
            plotThisTurn = True
        
        # calculate gradients
        rho_dx = getGradient(rho, dx)
        vx_dx = getGradient(vx,   dx)
        P_dx  = getGradient(P,    dx)
        Pi_dx = getGradient(Pi,   dx)
        
        # extrapolate half-step in time
        rho_prime = rho - 0.5*dt * ( vx * rho_dx + rho * vx_dx)
        vx_prime  = vx  - 0.5*dt * ( vx * vx_dx + np.divide(1 , rho, out=np.zeros_like(rho), where=rho!=0) * P_dx )
        P_prime   = P   - 0.5*dt * ( gamma*P * (vx_dx)  + vx * P_dx  )
        Pi_prime  = Pi  - 0.5*dt * (Pi_dx)
        
        # extrapolate in space to face centers
        # qM_XL, qP_XL, qM_XR, qP_XR

        rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR = extrapolateInSpaceToFace(rho_prime, rho_dx, dx)
        vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR  = extrapolateInSpaceToFace(vx_prime,  vx_dx,  dx)
        PM_XL,   PP_XL,   PM_XR,   PP_XR   = extrapolateInSpaceToFace(P_prime,   P_dx,   dx)
        PiM_XL,  PiP_XL,  PiM_XR,  PiP_XR  = extrapolateInSpaceToFace(Pi_prime,  Pi_dx,  dx)
        
        # compute fluxes (local Kurganov-Tadmor)
        # getFlux(rho_L, rho_R, vx_L, vx_R, Pi_L, Pi_R, P_L, P_R, gamma) inputs
        # flux_Mass, flux_Momx, flux_Energy, flux_Pi_v outputs 

        flux_Mass_XR, flux_Momx_XR, flux_Energy_XR, flux_Pi_vxR = getFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR, PiM_XR, PiP_XR, PP_XR, PM_XR, gamma, cs)
        flux_Mass_XL, flux_Momx_XL, flux_Energy_XL, flux_Pi_vxL = getFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL, PiM_XL, PiP_XL, PP_XL, PM_XL, gamma, cs)

        
        # update solution
        # applyFluxes(H, flux_H_X, dx)

        J = Pi*vx_dx - (zeta / tau_nu * vx_dx + Pi / tau_nu)

        Mass   = applyFluxes(Mass,   flux_Mass_XR,   dx) - applyFluxes(Mass, flux_Mass_XL, dx)
        Momx   = applyFluxes(Momx,   flux_Momx_XR,   dx) - applyFluxes(Momx, flux_Momx_XL, dx)
        Energy = applyFluxes(Energy, flux_Energy_XR, dx) - applyFluxes(Energy, flux_Energy_XL, dx)
        Pi_vx  = applyFluxes(Pi_vx,  flux_Pi_vxR,  dx, J) - applyFluxes(Pi_vx,  flux_Pi_vxL,  dx, J)
        
        # update time
        t += dt
        
        # plot in real time - color 1/2 particles blue, other half red
        if (plotRealTime and plotThisTurn) or (t >= tEnd):
            plt.cla()
            plt.imshow(rho.T)
            plt.clim(0.8, 2.2)
            ax = plt.gca()
            ax.invert_yaxis()
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)   
            ax.set_aspect('equal')  
            plt.pause(0.001)
            outputCount += 1
  
    
    # Save figure
    plt.savefig('finitevolume.png',dpi=240)

    plt.show()

    return 0



if __name__== "__main__":
  main()