#Visualization code 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_ic_csv(file, parameters_file):

    parameters = pd.read_csv(parameters_file)
    N = int(parameters['N'])
    boxsize = int(parameters['boxsize'])
    dx = boxsize/N
    xlin = np.linspace(0.5*dx, boxsize-0.5*dx,N)

    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    s = X.shape

    rho = np.zeros(s)
    Momx = np.zeros(s)
    Momy = np.zeros(s)  
    Pixx = np.zeros(s)
    Pixy = np.zeros(s)
    Piyx = np.zeros(s)
    Piyy = np.zeros(s)

    IC = pd.read_csv(file, header=None)

    ic = IC.to_numpy()[0]

    # get initial conditions
    for i in range(N):
        for j in range(N):
            rho[i][j]  = float(ic[7*(N*i + j)])
            Momx[i][j] = float(ic[7*(N*i + j) + 1])
            Momy[i][j] = float(ic[7*(N*i + j) + 2])
            Pixx[i][j] = float(ic[7*(N*i + j) + 3])
            Pixy[i][j] = float(ic[7*(N*i + j) + 4])
            Piyx[i][j] = float(ic[7*(N*i + j) + 5])
            Piyy[i][j] = float(ic[7*(N*i + j) + 6])
    
    plt.imshow(rho.T)
    plt.show()



def plot_solution_csv(file, parameters_file):

    parameters = pd.read_csv(parameters_file)
    N = int(parameters['N'])
    boxsize = int(parameters['boxsize'])
    dx = boxsize/N
    xlin = np.linspace(0.5*dx, boxsize-0.5*dx,N)

    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    s = X.shape

    rho = np.zeros(s)
    Momx = np.zeros(s)
    Momy = np.zeros(s)  
    Pixx = np.zeros(s)
    Pixy = np.zeros(s)
    Piyx = np.zeros(s)
    Piyy = np.zeros(s)

    print("starting...")
    solution = pd.read_csv(file, header=None)
    print("solution Dataframe created")

    s = 0
    outputcount = 0
    print("reading solution...")
    while s < solution.size():
        for i in range(N):
            for j in range(N):
                rho[i][j]  = float(solution[7*(N*i + j) + s])
                Momx[i][j] = float(solution[7*(N*i + j) + 1 + s])
                Momy[i][j] = float(solution[7*(N*i + j) + 2 + s])
                Pixx[i][j] = float(solution[7*(N*i + j) + 3 + s])
                Pixy[i][j] = float(solution[7*(N*i + j) + 4 + s])
                Piyx[i][j] = float(solution[7*(N*i + j) + 5 + s])
                Piyy[i][j] = float(solution[7*(N*i + j) + 6 + s])


        plt.imshow(rho.T)
        plt.show()
        outputcount += 1
        print(outputcount)
        s = 7*(N*i + j + 1)




def plot_each_csv(file, parameters_file):
    
    parameters = pd.read_csv(parameters_file)
    N = int(parameters['N'])
    boxsize = int(parameters['boxsize'])
    dx = boxsize/N
    xlin = np.linspace(0.5*dx, boxsize-0.5*dx,N)

    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    S = X.shape

    v = np.zeros(S)
    

    print("starting...")
    solution = pd.read_csv(file, header=None)
    print("solution Dataframe created")

    sol = solution.to_numpy()
    print(sol.shape)

    s = 0
    outputcount = 0
    print("reading solution...")
    while s < sol.shape[0]:
        for i in range(N):
            for j in range(N):
                v[i][j]  = float(sol[s][(N*i + j)])


        plt.imshow(v.T)
        plt.show()
        outputcount += 1
        print(outputcount)
        s +=1 



#plot_ic_csv('KT-ISshear\C++\initial_state.csv','KT-ISshear\C++\parameters.csv')
#plot_each_csv('KT-ISshear\C++\density_solution.csv','KT-ISshear\C++\parameters.csv')

sol = np.array(np.load("ShearKelvinHelmholtz.npy"))

print(sol.shape)

t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt("ShearKelvinHelmholtz_parameters", delimiter=',', unpack=True)

plotfinalstate = 1
N = int(N)
xlin = np.linspace(float(a),float(b),N)
#plt.imshow(sol[0][3*N:4*N].T)
#plt.show()

figure, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

# set ax boundaries
#ax1.set_xlim((4,6))
#ax1.set_ylim((0, 3))
line1, = ax1.plot([], [], lw=2)
ax1.set_ylabel('rho/rho_0')
ax1.set_title('Density')

#ax2.set_xlim((4,6))
#ax2.set_ylim((-1, 2))
line2, = ax2.plot([], [], lw=2)
ax2.set_title('Velocity')

#ax3.set_xlim((4,6))
#ax3.set_ylim((-1.5, 1.5))

line3, = ax3.plot([], [], lw=2)
ax3.set_xlabel('x/x_0')
ax3.set_ylabel('Pi/P_0')

line4, = ax4.plot([], [], lw=2)
ax4.set_xlabel('x/x_0')
ax4.set_ylabel('P/P_0')

if plotfinalstate == 1:
    x = xlin
    a = sol[-1][:N]
    b = sol[-1][N:2*N]
    c = sol[-1][2*N:3*N]
    d = (a)**gamma
    line1.set_data(x, a)
    line2.set_data(x, b)
    line3.set_data(x, c)
    line4.set_data(x, d)
    filename = "FinalStatenonRelativisticISwithShear.png"
    plt.savefig(str(filename))

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    return (line1,line2,line3,line4)

def animate(i):
    x = xlin
    a = sol[i][:N]
    b = sol[i][N:2*N]
    c = sol[i][2*N:3*N]
    d = (a)**gamma
    line1.set_data(x, a)
    line2.set_data(x, b)
    line3.set_data(x, c)
    line4.set_data(x, d)
    return (line1,line2,line3,line4)


ani = animation.FuncAnimation(figure, animate, init_func=init,
                                frames=200, interval=20, blit=True)

filename = "nonRelativisticISwithShear.gif"
plt.save(str(filename))

