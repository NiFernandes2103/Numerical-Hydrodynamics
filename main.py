import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d
from KTalgorithm import *
from EoS import *
from RK_Heuns_integrator import *

def plot_density(x,y,z):

    h = plt.contourf(x, y, z)
    plt.axis('scaled')
    plt.colorbar()
    plt.show()



def main():
    """ Finite Volume simulation """
    
    # Simulation parameters
    N                      = 128 # resolution
    boxsize                = 1.  # in some unit system l
    gamma                  = 5/3 # adiabatic index
    zeta                   = 0.2 # bulk viscosity coefficient
    tau_nu                 = 0.1
    cs                     = 1   # l/s
    t                      = 0   # s 
    tEnd                   = 2
    tOut                   = 0.02
    plotRealTime = True # switch on for plotting as the simulation goes along

    a = 0  #count
    
    # Mesh
    dx = boxsize / N
    vol = dx**2
    xlin = np.linspace(0.5*dx, boxsize-0.5*dx, N, dtype=np.float64 )
    Y, X = np.meshgrid( xlin, xlin)
    
    # Generate Initial Conditions  

    w0 = 0.1
    sigma = 5/np.sqrt(2.)
    rho = (1 - (X**2 + Y**2)/0.25 )**2 > 0
    # put Cubic spline in place of Gaussian
    plot_density(X,Y,rho)
    vx = np.zeros(X.shape)
    plot_density(X,Y,vx)
    P = 0.5 * np.ones(X.shape)
    Pi = np.zeros(X.shape)
    Pi_vx = Pi*vx

    # Get mean rho 
    times  = np.zeros(int(tEnd/tOut))
    rho_mean = np.mean(rho[:0]) * np.ones(int(tEnd/tOut))

    # Get conserved variables
    Mass, Momx, Energy = getConserved( rho, vx, P, gamma, vol)

    # prep figure
    fig = plt.figure(figsize=(4,4), dpi=80)
    outputCount = 1

    # Simulation Main Loop
    while t < tEnd:

        # get Primitive variables
        rho, vx, P = getPrimitive( Mass, Momx, Energy, gamma, vol )

        # get mean values (to plot)
        print(a,rho)
        a +=1
        times[a] = t 
        rho_mean[a] = t

        # define timestep
        dt = dx/cs
      
        plotThisTurn = False
        if t + dt > outputCount*tOut:
            dt = outputCount*tOut - t
            plotThisTurn = True
        
        # calculate gradients
        rho_dx = getGradient(rho, dx)
        print(rho_dx)
        vx_dx = getGradient(vx,   dx)
        P_dx  = getGradient(P,    dx)
        Pi_dx = getGradient(Pi,   dx)
        
        # extrapolate half-step in time
        rho_prime = rho - 0.5*dt * ( vx * rho_dx + rho * vx_dx)
        vx_prime  = vx  - 0.5*dt * np.divide(( vx * vx_dx * rho + P_dx ), rho , out=np.zeros_like(vx * vx_dx * rho + P_dx), where=rho!=0)
        P_prime   = P   - 0.5*dt * ( gamma*P * (vx_dx)  + vx * P_dx  )
        Pi_prime  = Pi  - 0.5*dt * (Pi_dx)
        
        # extrapolate in space to face centers

        rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR = extrapolateInSpaceToFace(rho_prime, rho_dx, dx)
        vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR  = extrapolateInSpaceToFace(vx_prime,  vx_dx,  dx)
        PM_XL,   PP_XL,   PM_XR,   PP_XR   = extrapolateInSpaceToFace(P_prime,   P_dx,   dx)
        PiM_XL,  PiP_XL,  PiM_XR,  PiP_XR  = extrapolateInSpaceToFace(Pi_prime,  Pi_dx,  dx)
        
        # compute fluxes (local Kurganov-Tadmor)

        flux_Mass_XR, flux_Momx_XR, flux_Pi_vxR = getFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR, PiM_XR, PiP_XR, PP_XR, PM_XR, gamma, cs)
        flux_Mass_XL, flux_Momx_XL, flux_Pi_vxL = getFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL, PiM_XL, PiP_XL, PP_XL, PM_XL, gamma, cs)
        
        # update solution

        J = Pi*vx_dx - (zeta / tau_nu * vx_dx + Pi / tau_nu)

        Mass   =  modified_RungeKutta( Mass, vol* applyFluxes(Mass,   flux_Mass_XR,   flux_Mass_XL,   dx), dt)
        SMass = np.sum(Mass)
        print('Mass',SMass)
        Momx   =  modified_RungeKutta( Momx, vol * applyFluxes(Momx,   flux_Momx_XR,   flux_Momx_XL,   dx), dt)  
        SMomx = np.sum(Momx)
        print('Momentum',SMomx)
        Pi_vx  =  modified_RungeKutta( Pi_dx , vol * applyFluxes(Pi_vx,  flux_Pi_vxR,    flux_Pi_vxL,    dx, J))
        

        # effectively zero terms
        threshold = 10^-10
        Mass = Mass > threshold
        Momx  = Momx  > threshold
        Pi_vx  = Pi_vx  > threshold
        Energy = Energy > threshold

        # Boundary conditions
        
        '''
        rho[0][:] = rho[3][:]
        rho[1][:] = rho[3][:]        
        rho[N-1][:] = rho[N-4][:]
        rho[N-2][:] = rho[N-4][:]
        rho[:][0] = rho[:][3]
        rho[:][1] = rho[:][3]        
        rho[:][N-1] = rho[:][N-4]
        rho[:][N-2] = rho[:][N-4]
        '''

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
            
        else:
           plt.plot(rho_mean, times)
  
    
    # Save figure
    plt.savefig('finitevolume.png',dpi=240)

    plt.show()

    return 0

  

if __name__== "__main__":
  main()