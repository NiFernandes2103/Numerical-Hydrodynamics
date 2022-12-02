import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d
from KTalgorithm import *
from EoS import *
from RK_Heuns_integrator import *


def main():
    """ Finite Volume simulation """
    
    # Simulation parameters
    N                      = 128 # resolution
    boxsize                = 1.  # in some unit system l
    gamma                  = 5/3 # adiabatic index
    zeta                   = 0 # bulk viscosity coefficient
    tau_nu                 = 200
    t                      = 0   # s 
    tEnd                   = 20   # time at the end
    tOut                   = 2 # time of each output
    plotRealTime = True  # switch on for plotting as the simulation goes along

    #test parameters
    flimit                 = 10**4  # Warns if fluxes are larger than this value
    runFluxLimitTests      = False  # Switch to conduct tests

    a = 0  #count
    
    #-----------------------------------------------------------------------------------------------------------------------------------#

    # Define Mesh
    dx = boxsize / N   # box size
    vol = dx**2        # volume of each box
    xlin = np.linspace(0.5*dx, boxsize-0.5*dx, N, dtype=np.float64 ) # simulation limits


    #-----------------------------------------------------------------------------------------------------------------------------------#

    # Generate Initial Conditions  


    #rho0 = np.ones(xlin.shape)
    rho0 = ((1 - ((xlin - (boxsize-0.5*dx)*0.5)**2)/0.25 )**4 ) # Mauricio`s funtion advice   
    #rho0 = 1*(xlin < boxsize*0.5) + 0.125*(xlin >= boxsize*0.5)
    

    # put Cubic spline in place of Gaussian
  

  
    #vx0 = np.abs((xlin - (boxsize-0.5*dx)*0.5)/16)
    vx0 = np.zeros(xlin.shape)
    #vx0 = np.sin(xlin)*np.ones(xlin.shape)

    P0 = (np.abs(rho0))**gamma
    Pi0 = np.zeros(xlin.shape)

    #-----------------------------------------------------------------------------------------------------------------------------------#

    rho = rho0
    vx = vx0
    P = P0
    Pi = Pi0
    cs = getSpeedOfSound(rho,gamma)
    

    #-----------------------------------------------------------------------------------------------------------------------------------#
  

    # Get mean rho 
    
    time  = np.zeros(int(tEnd/0.01)+1)
    #rho_mean = np.mean(rho[:0]) * np.ones(int(tEnd/0.01)+1)

    #-----------------------------------------------------------------------------------------------------------------------------------#


    # Get conserved variables
    Mass, Momx = getConserved( rho, vx, vol)

    #-----------------------------------------------------------------------------------------------------------------------------------#


    # prep figure
    fig = plt.figure(figsize=(4,4), dpi=80)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.title.set_text('Density')
    ax2.title.set_text('Velocity')

    ax1.plot(xlin, rho0)
    ax2.plot(xlin, vx0)

    outputCount = 1

    #-----------------------------------------------------------------------------------------------------------------------------------#


    # Simulation Main Loop
    while t < tEnd:

        #-----------------------------------------------------------------------------------------------------------------------------------#

        # get Conserved variables
        Mass, Momx = getConserved( rho, vx, vol)

        if np.any(Mass < 0):
          print("Mass is negative")
          exit()

        # get Primitive variables
        rho, vx, P = getPrimitive( Mass, Momx, gamma, vol )

        #-----------------------------------------------------------------------------------------------------------------------------------#

        # get mean values (to plot)
        #a +=1
        #time[a] = t 
        #rho_mean[a] = t
  
        cs = getSpeedOfSound(rho, gamma)

        #-----------------------------------------------------------------------------------------------------------------------------------#

        # get time step

        dt = (dx/(np.amax(cs)))
      
        plotThisTurn = False

        if t + dt > outputCount*tOut or t == 0 :
            plotThisTurn = True
        
        #-----------------------------------------------------------------------------------------------------------------------------------#
        
        # calculate gradients
        rho_dx = getGradient(rho, dx)
        vx_dx = getGradient(vx,   dx)
        P_dx  = getGradient(P,    dx)
        Pi_dx = getGradient(Pi,   dx)

       
        #-----------------------------------------------------------------------------------------------------------------------------------#
        
        # extrapolate in space to face centers

        rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR = extrapolateInSpaceToFace(rho, rho_dx, dx)
        vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR  = extrapolateInSpaceToFace(vx,  vx_dx,  dx)
        PM_XL,   PP_XL,   PM_XR,   PP_XR   = extrapolateInSpaceToFace(P,   P_dx,   dx)
        PiM_XL,  PiP_XL,  PiM_XR,  PiP_XR  = extrapolateInSpaceToFace(Pi,  Pi_dx,  dx)

        #-----------------------------------------------------------------------------------------------------------------------------------# 
        
        # compute fluxes (local Kurganov-Tadmor)

        flux_Mass_XR, flux_Momx_XR, flux_Pi_vxR = getFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR, PiP_XR, PiM_XR, PP_XR, PM_XR, gamma)
        flux_Mass_XL, flux_Momx_XL, flux_Pi_vxL = getFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL, PiP_XL, PiM_XL, PP_XL, PM_XL, gamma)

        # tests for flux limitation
        if runFluxLimitTests:
          print(np.any(flux_Mass_XR > flimit), np.any(flux_Momx_XR > flimit), np.any(flux_Pi_vxR > flimit))
          print(np.any(flux_Mass_XL > flimit), np.any(flux_Momx_XL > flimit), np.any(flux_Pi_vxL > flimit))

        #-----------------------------------------------------------------------------------------------------------------------------------#
        
        # update solution

        J = (Pi * vx_dx - (zeta / tau_nu * vx_dx + Pi / tau_nu))
        

        rho    = modified_RungeKutta(rho,   applyFluxes( flux_Mass_XR,   flux_Mass_XL,   dx),    dt)
        vx     = np.divide(modified_RungeKutta(Momx,  applyFluxes( flux_Momx_XR,   flux_Momx_XL,   dx),    dt), rho,
                         out=np.zeros_like(modified_RungeKutta(Momx,  applyFluxes( flux_Momx_XR,   flux_Momx_XL,   dx), dt)), where=rho!=0)
        Pi     = modified_RungeKutta(Pi, applyFluxes( flux_Pi_vxR,    flux_Pi_vxL,    dx, J), dt)

        Pi = regParams(Pi,1)
        

        #-----------------------------------------------------------------------------------------------------------------------------------#

        # Boundary conditions
        

        # Periodic boundary conditions
      
        rho[0]    = rho[3]
        rho[1]    = rho[3]        
        rho[N-1]  = rho[N-4]
        rho[N-2]  = rho[N-4]

        vx[0]    = -vx[3]
        vx[1]    = -vx[3]        
        vx[N-1]  = -vx[N-4]
        vx[N-2]  = -vx[N-4]
        
        Pi[0]    = 0
        Pi[1]    = 0
        Pi[N-1]  = 0
        Pi[N-2]  = 0

        #-----------------------------------------------------------------------------------------------------------------------------------#

        # update time
        t += dt

        #-----------------------------------------------------------------------------------------------------------------------------------#
        
        # plot in real time 
        if (plotRealTime and plotThisTurn) or (t >= tEnd):
            ax1.plot(xlin,rho)
            ax1.set_ylim((0,1+0.1))
            ax2.plot(xlin,vx)
            plt.xlim((0,boxsize))
            plt.pause(0.001)
            outputCount += 1
            
       #-----------------------------------------------------------------------------------------------------------------------------------#    
  
           
    #plt.plot(rho_mean, time)
        
    # Save figure
  
    
    plt.savefig('finitevolume_PolynomialDensity_.png',dpi=240)

    plt.show()

    return 0

  

if __name__== "__main__":
  main()