import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

sol = np.load("PerfectFluidRelativistic.npy")

print(np.sqrt(1-0.99*2))
t,tEnd,tOut,N,boxsize,gamma,theta,a,b = np.loadtxt("PerfectFluidRelativistic_parameters",delimiter=',',unpack=True)

N = int(N)
xlin = np.linspace(float(a),float(b),N)
print(sol.shape)
i=0

Dic = sol[0][0:N]
Sxic = sol[0][N:2*N]
tauic = sol[0][2*N:3*N]
D = sol[i][0:N]
Sx = sol[i][N:2*N]
tau = sol[i][2*N:3*N]
'''
plt.plot(xlin,Dic,label='t = 0.0', linestyle='dashed')
plt.plot(xlin, D,label='t = {:.2f}'.format(float(i*tOut)))
plt.title('rest-mass density')
plt.xlim((float(a),float(b)))
plt.ylim(0,10)
plt.xlabel('x/x_0')
plt.ylabel('m/m_0')
plt.legend()
plt.savefig('PerfectFluidRestMassDensity.png')
'''
'''
plt.plot(xlin,Sxic,label='t = 0.0', linestyle = 'dashed')
plt.plot(xlin,Sx,label='t = {:.2f}'.format(float(i*tOut)))
plt.title('x-Momentum density')
plt.xlabel('x/x_0')
plt.ylabel('Sx/Sx_0')
plt.legend()
plt.savefig('PerfectFluidxMomentumDensity.png')
'''
'''
plt.plot(xlin,tauic, label='t = 0.0', linestyle='dashed')
plt.plot(xlin,tau, label='t = {:.2f}'.format(float(i*tOut)))
plt.title('Rest-Energy Density')
plt.xlabel('x/x_0')
plt.ylabel('v/c')
plt.legend()
plt.savefig('PerfectFluidRestEnergyDensity.png')
'''



figure, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

# set ax boundaries
ax1.set_ylim(0,10)
ax1.set_xlim((float(a),float(b)))
line1, = ax1.plot([], [], lw=2)
ax1.set_ylabel('D/D_0')
ax1.set_title('Rest-mass density')

ax2.set_ylim(-10,500)
ax2.set_xlim((float(a),float(b)))
line2, = ax2.plot([], [], lw=2)
ax2.set_title('x-Momentum density')

ax3.set_ylim(-10,60)
ax3.set_xlim((float(a),float(b)))
line3, = ax3.plot([], [], lw=2)
ax3.set_xlabel('x/x_0')
ax3.set_ylabel('tau/tau_0')


x = xlin
a = sol[0][:N]
b = sol[0][N:2*N]
c = sol[0][2*N:3*N]
line1.set_data(x, a)
line2.set_data(x, b)
line3.set_data(x, c)
filename = "InitialStateRelativisticPerfectFluid{}.png".format(i)
plt.savefig(str(filename))

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    return (line1,line2,line3)


def animate(i):
    x = xlin
    a = sol[i][:N]
    b = sol[i][N:2*N]
    c = sol[i][2*N:3*N]
    line1.set_data(x, a)
    line2.set_data(x, b)
    line3.set_data(x, c)
    return (line1,line2,line3)


ani = animation.FuncAnimation(figure, animate, init_func=init,
                    frames=100, interval=20, blit=True)

filename = "RelativisticPerfectFluid{}.gif".format(i)
ani.save(str(filename)) 