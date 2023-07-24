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
    a = float(parameters['a'])
    b = float(parameters['b'])
    dx = boxsize/N
    xlin = np.linspace(a, b,N)

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
        plt.plot(xlin, v[int(N/2)].T)
        plt.show()
        outputcount += 1
        print(outputcount)
        s += 10



#plot_ic_csv('KT_ISshear\C++\initial_state.csv','KT_ISshear\C++\parameters.csv')
plot_each_csv('KT_ISshear\C++\density_solution.csv','KT_ISshear\C++\parameters.csv')
#plot_each_csv('KT_ISshear\C++\Pixy_solution.csv','KT_ISshear\C++\parameters.csv')

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

def animate_solution_gif(file,parameters_file,gif_file,n):

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


def animate_numpy_solution_gif(sol,N):


    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes()
    line, = ax.plot([], [], lw=2)
    s = 0
    outputcount = 0
    v = np.zeros((N,N))
    b = []
    print("reading solution...")
    while s < sol.shape[0]:
        for i in range(N):
            for j in range(N):
                v[i][j]  = sol[s][(N*i + j)]

        b.append(v)
        s += 1
    print("finished reading")

    a = b[0].T
    b = np.array(b)
    print(b.shape)

    im=plt.imshow(a,interpolation='none')

    # initialization function: plot the background of each frame
    def init():
        im.set_data(a)
        return [im]

    # animation function.  This is called sequentially
    def animate(i):
        im.set_array((b[i]).T)
        return [im]


    return animation.FuncAnimation(fig, animate, init_func=init,
                                frames=100, interval=20, blit=True)

    



def animate_numpy_solution_slice_gif(sol,N,xlin):


    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes()
    line, = ax.plot([], [], lw=2)
    s = 0
    outputcount = 0
    v = np.zeros((N,N))
    b = []
    print("reading solution...")
    while s < sol.shape[0]:
        for i in range(N):
            for j in range(N):
                v[i][j]  = sol[s][(N*i + j)]

        b.append(v)
        s += 1
    print("finished reading")

    b = np.array(b)
    a = b[0][int(N/2)].T


    # initialization function: plot the background of each frame
    def init():
        line.set_data(xlin,a)
        return (line,)

    # animation function.  This is called sequentially
    def animate(i):
        line.set_array(xlin,b[i][int(N/2)].T)
        return (line,)


    return animation.FuncAnimation(fig, animate, init_func=init,
                                frames=40, interval=100, blit=True)

    
                            


'''
parameters = pd.read_csv('KT_ISshear\C++\parameters.csv')
N = int(parameters['N'])
boxsize = int(parameters['boxsize'])
a = float(parameters['a'])
b = float(parameters['b'])

dx = boxsize/N
xlin = np.linspace(a, b, N)

Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
S = X.shape

v = np.zeros(S)

sol = (pd.read_csv('KT_ISshear\C++\density_solution.csv', header=None)).to_numpy()


# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes()
#line, = ax.plot([], [], lw=2)
s = 0
outputcount = 0

b = []
print("reading solution...")
while s < sol.shape[0]:
    for i in range(N):
        for j in range(N):
            v[i][j]  = sol[s][(N*i + j)]

    b.append(v)
    s += 1
print("finished reading")

b = np.array(b)
print(b.shape)

im=plt.imshow(b[0].T,interpolation='none')

# initialization function: plot the background of each frame
def init():
    im.set_data(b[0].T)
    return [im]

# animation function.  This is called sequentially
def animate(i):
    im.set_array((b[i]).T)
    return [im]


ani = animation.FuncAnimation(fig, animate, init_func=init,
                            frames=100, interval=20, blit=True)


ani.save("NonRelativisticISC++.gif",fps=30)
'''

'''
sol = np.load("NonRelativisticISHeuns.npy")

t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt("NonRelativisticISHeuns_parameters",delimiter=',',unpack=True)

xlin = np.linspace(float(a),float(b),int(N))
print(sol.shape)

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
    im.set_data(sol[0][:00].T)
    return [im]

# animation function.  This is called sequentially
def animate(i):
    im.set_array(sol[i][:200].T)
    return [im]


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=100, interval=20, blit=True)


anim.save('NonRelativisticISHeuns.gif', fps=30)


# First set up the figure, the axis, and the plot element we want to animate
fig_slice = plt.figure()
ax_slice = plt.axes()
ax_slice.set_xlim((-2,2))
ax_slice.set_ylim((1, 2))
line, = ax_slice.plot([], [], lw=2)




# initialization function: plot the background of each frame
def init_slice():
    line.set_data(xlin,sol[0][150].T)
    return (line,)

# animation function.  This is called sequentially
def animate_slice(i):
    line.set_data(xlin,sol[i][150].T)
    return (line,)


anim_slice = animation.FuncAnimation(fig_slice, animate_slice, init_func=init_slice,
                            frames=100, interval=20, blit=True)

anim_slice.save('NonRelativisticISHeuns_density_slice.gif', fps=30)
'''