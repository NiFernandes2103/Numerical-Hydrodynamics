#Visualization code 

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
        plt.clim(0.8,2.2)
        plt.show()
        outputcount += 1
        print(outputcount)
        s +=1 



#plot_ic_csv('KT-ISshear\C++\initial_state.csv','KT-ISshear\C++\parameters.csv')
#plot_each_csv('KT-ISshear\C++\Pixx_solution.csv','KT-ISshear\C++\parameters.csv')
#plot_each_csv('KT-ISshear\C++\Piyy_solution.csv','KT-ISshear\C++\parameters.csv')


def show_2dsolution_static(file,parameters_file,i,n):

    sol = np.load(file)

    t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt(parameters_file,delimiter=',',unpack=True)

    N = int(N)
    a = sol[i][n*N:(n+1)*N].T
    
    plt.imshow(a)
    plt.show()
    

def show_solution_static_slice(file,parameters_file,i,n):

    sol = np.load(file)

    t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt(parameters_file,delimiter=',',unpack=True)

    N = int(N)
    a = sol[i][n*N:(n+1)*N].T
    
    plt.imshow(a[int(N/2)])
    plt.show()

def animate_solution_gif(file,parameters_file,gif_file):

    sol = np.load(file)

    t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt(parameters_file,delimiter=',',unpack=True)

    N = int(N)

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes()
    #line, = ax.plot([], [], lw=2)
    a = sol[0][n*N:(n+1)*N].T

    im=plt.imshow(a,interpolation='none')

    # initialization function: plot the background of each frame
    def init():
        im.set_data(a)
        return [im]

    # animation function.  This is called sequentially
    def animate(i):
        im.set_array(sol[i][n*N:(n+1)*N].T)
        return [im]


    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=1000, interval=20, blit=True)


    anim.save(gif_file, fps=30)



sol = np.load("NonRelativisticISgamma1.npy")

i=10
rho = sol[i][:200].T
vx = sol[i][200:2*200].T
vy = sol[i][2*200:3*200].T

plt.imshow(rho)
plt.show()
plt.plot(rho[int(200/2)])
plt.show()
plt.imshow(vx)
plt.show()
plt.plot(vx[int(200/2)])
plt.show()
plt.imshow(vy)
plt.show()
plt.plot(vy[int(200/2)])
plt.show()

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes()
#line, = ax.plot([], [], lw=2)
im=plt.imshow(sol[0][:200].T,interpolation='none')

# initialization function: plot the background of each frame
def init():
    im.set_data(sol[0][:200].T)
    return [im]

# animation function.  This is called sequentially
def animate(i):
    im.set_array(sol[i][:200].T)
    return [im]


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=2000, interval=20, blit=True)


anim.save('NonRelativisticISgamma1.gif', fps=30)

plt.show()