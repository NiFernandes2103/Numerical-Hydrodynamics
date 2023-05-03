import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import animation, rc
from IPython.display import HTML
from scipy import integrate
from KTalgorithm import *
from EoS import *
from RK_Heuns_integrator import *

def applyBC(y):

  rho = y[0:N] 
  vx = y[N:2*N]/rho 
  Pi = y[2*N:] 

  '''Absorbing boundary conditions'''
  rho[0]    = rho[1]
  rho[-1]   = rho[-2]

  vx[0]     = vx[1]      
  vx[-1]    = vx[-2]

  Pi[0]     = Pi[1]   
  Pi[-1]    = Pi[-2]

t                      = 0   # s 
tEnd                   = 2   # time at the end
tOut                   = 0.01 # time of each output

N                      = 400 # resolution
boxsize                = 1.  # in some unit system l
gamma                  = 4 # adiabatic index
zeta                   = 1 # bulk viscosity coefficient
tau_nu                 = 1
theta                  = 1


# Define Mesh
dx = boxsize / N   # box size
vol = dx**2        # volume of each box
xlin = np.linspace(0.5*dx, boxsize-0.5*dx, N) # simulation limits


rho = ((1 - ((xlin - (boxsize-0.5*dx)*0.5)**2)/0.25 )**4 ) + 0.5*np.ones(xlin.shape) # Mauricio`s funtion advice
#rho  = 0.5*np.sin(2*np.pi*xlin) + 1*np.ones(xlin.shape)    
#rho = 1*(xlin < boxsize*0.5) + 0.125*(xlin >= boxsize*0.5)


vx = np.zeros(xlin.shape)
#vx = 0.5*np.ones(xlin.shape)
#vx = np.abs((xlin - (boxsize-0.5*dx)*0.5)/16)


Pi = np.zeros(xlin.shape)

IC = np.hstack((rho,rho*vx,Pi)) # here the initial conditions are stacked in a vector 
                            # rho is IC[0:N] for example


solution = integrator(KTschemeNonRelativisticIS, (t,tEnd), IC, 0.01, applyBC, args=(dx, xlin, gamma, zeta, tau_nu, None, theta))


figure = plt.figure()
ax1 = plt.subplot(1,1,1)

# set ax boundaries
ax1.set_xlim((0,boxsize))
ax1.set_ylim((-3, 3))
line, = ax1.plot([], [], lw=2)
ax1.set_xlabel('x/x_0')
ax1.set_ylabel('rho/rho_0')
ax1.set_title('non-relativistic Israel-Stewart equation')

def init():
    line.set_data([], [])
    return (line,)


def animate_density(i):
  x = xlin
  y = solution[i][2*N:]
  line.set_data(x, y)
  return (line,)


ani = animation.FuncAnimation(figure, animate_density, init_func=init,
                               frames=1000, interval=20, blit=True)

ani.save("nonRelativisticIS.gif")

#HTML(ani.to_html5_video())