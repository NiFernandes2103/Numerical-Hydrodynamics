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

  Pi[0]     = 0
  Pi[-1]    = 0

plotfinalstate = 1

for i in range(10):
  t                      = 0   # s 
  tEnd                   = 2   # time at the end
  tOut                   = 0.01 # time of each output

  if i == 1:
    N                      = 1000 # resolution
  else:
    N                      = 400 # resolution
    
  boxsize                = 10.  # in some unit system l
  gamma                  = 1 # adiabatic index
  zeta                   = 1 # bulk viscosity coefficient
  tau_nu                 = 1
  theta                  = 1
  

  # Define Mesh
  dx = boxsize / N   # box size
  vol = dx**2        # volume of each box
  xlin = np.linspace(0.5*dx, boxsize-0.5*dx, N) # simulation limits
  if i == 0 or i == 1 or i == 4 or i == 7:
    rho = ((1 - ((xlin - (boxsize-0.5*dx)*0.5)**2)/0.1 )**4 )*(abs((xlin-(boxsize-0.5*dx)*0.5)) < 0.2) + 0.5*np.ones(xlin.shape) # Mauricio`s funtion advice
  elif i == 2 or i == 5 or i == 8:
    rho  = 0.5*np.sin(2*np.pi*xlin) + 1*np.ones(xlin.shape)    
  else:
    rho = 1*(xlin < boxsize*0.5) + 0.125*(xlin >= boxsize*0.5)

  if i == 0 or i == 1 or i == 2 or i == 3:
    vx = np.zeros(xlin.shape)
  if i == 4 or i == 5 or i == 6:
    vx = 0.5*np.ones(xlin.shape)
  else:
    vx = np.abs((xlin - (boxsize-0.5*dx)*0.5)/16)


  Pi = np.zeros(xlin.shape)

  IC = np.hstack((rho,rho*vx,Pi)) # here the initial conditions are stacked in a vector 
                              # rho is IC[0:N] for example


  solution = integrator(KTschemeNonRelativisticIS, (t,tEnd), IC, 0.01, applyBC, args=(dx, xlin, gamma, zeta, tau_nu, None, theta))


  figure, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

  # set ax boundaries
  ax1.set_xlim((4,6))
  ax1.set_ylim((0, 3))
  line1, = ax1.plot([], [], lw=2)
  ax1.set_ylabel('rho/rho_0')
  ax1.set_title('Density')

  ax2.set_xlim((4,6))
  ax2.set_ylim((-1, 2))
  line2, = ax2.plot([], [], lw=2)
  ax2.set_title('Velocity')

  ax3.set_xlim((4,6))
  ax3.set_ylim((-1.5, 1.5))

  line3, = ax3.plot([], [], lw=2)
  ax3.set_xlabel('x/x_0')
  ax3.set_ylabel('Pi/P_0')

  line4, = ax4.plot([], [], lw=2)
  ax4.set_xlabel('x/x_0')
  ax4.set_ylabel('P/P_0')

  if plotfinalstate == 1:
    x = xlin
    a = solution[0][:N]
    b = solution[0][N:2*N]
    c = solution[0][2*N:]
    d = (np.abs(a))**gamma
    line1.set_data(x, a)
    line2.set_data(x, b)
    line3.set_data(x, c)
    line4.set_data(x, d)
    filename = "FinalStatenonRelativisticIS{}.png".format(i)
    plt.savefig(str(filename))

  def init():
      line1.set_data([], [])
      line2.set_data([], [])
      line3.set_data([], [])
      line4.set_data([], [])
      return (line1,line2,line3,line4)


  def animate(i):
    x = xlin
    a = solution[i][:N]
    b = solution[i][N:2*N]
    c = solution[i][2*N:]
    d = (np.abs(a))**gamma
    line1.set_data(x, a)
    line2.set_data(x, b)
    line3.set_data(x, c)
    line4.set_data(x, d)
    return (line1,line2,line3,line4)


  #ani = animation.FuncAnimation(figure, animate, init_func=init,
  #                              frames=200, interval=20, blit=True)

  #filename = "nonRelativisticIS{}.png".format(i)
  #plt.save(str(filename))

