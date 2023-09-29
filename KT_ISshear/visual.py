#Visualization code 

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

#plt.rcParams['text.usetex'] = True    # used to display Latex

''' Methods for reading and ploting for C++ code'''

def plot_ic_csv(file, parameters_file,variable_number):

    #plot initial conditions for a variable with number = variable_number

    # read parameter file
    parameters = pd.read_csv(parameters_file+'.csv')

    #intialize parameters
    N = int(parameters['N'])
    boxsize = int(parameters['boxsize'])
    a = float(parameters['a'])
    b = float(parameters['b'])
    dx = boxsize/N
    xlin = np.linspace(a,b,N)
    gamma = float(parameters['gamma'])
    zeta = float(parameters['zeta'])
    eta = float(parameters['eta'])
    tau_nu = float(parameters['tau_nu'])
    theta = float(parameters['theta'])
    xlin = np.linspace(a, b,N)
    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    s = X.shape

    #initialize variables
    rho = np.zeros(s)
    Momx = np.zeros(s)
    Momy = np.zeros(s)  
    Pixx = np.zeros(s)
    Pixy = np.zeros(s)
    Piyx = np.zeros(s)
    Piyy = np.zeros(s)

    #read initial conditions file
    IC = pd.read_csv(file+'.csv', header=None)

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
    
    if variable_number == 0:
        plt.imshow(rho.T)
    elif variable_number == 1:
        plt.imshow(Momx.T/rho.T)
    elif variable_number == 2:
        plt.imshow(Momy.T/rho.T)
    elif variable_number == 3:
        plt.imshow(Pixx.T)
    elif variable_number == 4:
        plt.imshow(Pixy.T)
    elif variable_number == 5:
        plt.imshow(Piyx.T)
    elif variable_number == 6:
        plt.imshow(Piyy.T)

    plt.show()



def plot_solution_csv(file, parameters_file, variable_number):

    #plot solution each frame for a variable with number = variable_number

    # read parameter file
    parameters = pd.read_csv(parameters_file+'.csv')
    
    #intialize parameters
    N = int(parameters['N'])
    boxsize = int(parameters['boxsize'])
    a = float(parameters['a'])
    b = float(parameters['b'])
    dx = boxsize/N
    xlin = np.linspace(a,b,N)
    gamma = float(parameters['gamma'])
    zeta = float(parameters['zeta'])
    eta = float(parameters['eta'])
    tau_nu = float(parameters['tau_nu'])
    theta = float(parameters['theta'])
    xlin = np.linspace(a, b,N)
    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    s = X.shape

    #initialize variables
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

        if variable_number == 0:
            plt.imshow(rho.T)
        elif variable_number == 1:
            plt.imshow(Momx.T/rho.T)
        elif variable_number == 2:
            plt.imshow(Momy.T/rho.T)
        elif variable_number == 3:
            plt.imshow(Pixx.T)
        elif variable_number == 4:
            plt.imshow(Pixy.T)
        elif variable_number == 5:
            plt.imshow(Piyx.T)
        elif variable_number == 6:
            plt.imshow(Piyy.T)
        plt.show()
        outputcount += 1
        print(outputcount)
        s = 7*(N*i + j + 1)



'''
def plot_each_csv(file, parameters_file):
    
    parameters = pd.read_csv(parameters_file)
    N = int(parameters['N'])
    boxsize = int(parameters['boxsize'])
    a = float(parameters['a'])
    b = float(parameters['b'])
    dx = boxsize/N
    xlin = np.linspace(a,b,N)

    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    S = X.shape
    previous = np.zeros(S)
    v = np.zeros(S)
    after = np.zeros(S)
    

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
                previous[i][j]  = float(sol[s-10][(N*i + j)])
                v[i][j]  = float(sol[s][(N*i + j)])
                after[i][j]  = float(sol[s+10][(N*i + j)])


        plt.imshow(v.T)
        plt.clim(0.8,2.2)
        plt.show()
        plt.plot(xlin, v[int(N/2)].T)
        plt.show()
        outputcount += 1
        print(outputcount)
        s += 10
'''
def plot_image(v, ax, fontsize=12, hide_labels=False):
        # auxiliary function used for plotting
        pc = ax.pcolormesh(v, vmin=1, vmax=2)
        if not hide_labels:
            ax.set_xlabel('x/x_0', fontsize=fontsize)
            ax.set_ylabel('rho/rho_0', fontsize=fontsize)
            ax.set_title('density', fontsize=fontsize)
        return pc


def plot_csv_static(file1, file2, file3, file4, file5, file6, file7, parameters_file, frame, nameOfFigure, hide_labels = False):

    # plot using separate files
    
    parameters = pd.read_csv(parameters_file + ".csv")

    #intialize parameters
    N = int(parameters['N'])
    boxsize = int(parameters['boxsize'])
    a = float(parameters['a'])
    b = float(parameters['b'])
    dx = boxsize/N
    xlin = np.linspace(a,b,N)
    gamma = float(parameters['gamma'])
    zeta = float(parameters['zeta'])
    eta = float(parameters['eta'])
    tau_nu = float(parameters['tau_nu'])
    theta = float(parameters['theta'])

    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    R = np.sqrt(X*X+Y*Y)
    S = X.shape

    # initialize variables used to store solutions
    init1 = np.zeros(S)
    v1    = np.zeros(S)
    init2 = np.zeros(S)
    v2    = np.zeros(S)
    init3 = np.zeros(S)
    v3    = np.zeros(S)
    init4 = np.zeros(S)
    v4    = np.zeros(S)
    init5 = np.zeros(S)
    v5    = np.zeros(S)
    init6 = np.zeros(S)
    v6    = np.zeros(S)
    init7 = np.zeros(S)
    v7    = np.zeros(S)

    print("starting...")
    solution1 = pd.read_csv(file1+".csv", header=None)
    solution2 = pd.read_csv(file2+".csv", header=None)
    solution3 = pd.read_csv(file3+".csv", header=None)
    solution4 = pd.read_csv(file4+".csv", header=None)
    solution5 = pd.read_csv(file5+".csv", header=None)
    solution6 = pd.read_csv(file6+".csv", header=None)
    solution7 = pd.read_csv(file7+".csv", header=None)
    print("solution Dataframe created")

    # convert to numpy array
    sol1 = solution1.to_numpy()
    sol2 = solution2.to_numpy()
    sol3 = solution3.to_numpy()
    sol4 = solution4.to_numpy()
    sol5 = solution5.to_numpy()
    sol6 = solution6.to_numpy()
    sol7 = solution7.to_numpy()

    print("solution format:{}".format(sol1.shape))

    print("reading solution...")
    for i in range(N):
        for j in range(N):
            init1[i][j] = float(sol1[0][N*i + j])
            v1[i][j]  = float(sol1[frame][N*i + j])
            init2[i][j] = float(sol2[0][N*i + j])
            v2[i][j]  = float(sol2[frame][N*i + j])
            init3[i][j] = float(sol3[0][N*i + j])
            v3[i][j]  = float(sol3[frame][N*i + j])
            init4[i][j] = float(sol4[0][N*i + j])
            v4[i][j]  = float(sol4[frame][N*i + j])
            init5[i][j] = float(sol5[0][N*i + j])
            v5[i][j]  = float(sol5[frame][N*i + j])
            init6[i][j] = float(sol6[0][N*i + j])
            v6[i][j]  = float(sol6[frame][N*i + j])
            init7[i][j] = float(sol7[0][N*i + j])
            v7[i][j]  = float(sol7[frame][N*i + j])
            
            
    print('finished reading')

    fig, axs = plt.subplots(2, 4, figsize=(10, 6))
    gridspec = axs[0, 0].get_subplotspec().get_gridspec()
    print("created the figure")

    # clear the left column for the subfigure:

    # plot data in remaining axes:
    count = 0
    print("plotting subfigures")
    for a in axs[:, :].flat:
        if count == 0:
            a.plot(xlin, init1.T[int(N/2)]) 
            a.plot(xlin, v1.T[int(N/2)])
            a.set_title("Density")
            count += 1
        elif count == 1:

            urinit = np.sqrt((init3.T[int(N/2)] / init1.T[int(N/2)])**2 + (init2.T[int(N/2)] / init1.T[int(N/2)])**2)
            ur = np.sqrt((v3.T[int(N/2)] / v1.T[int(N/2)])**2 + (v2.T[int(N/2)] / v1.T[int(N/2)])**2)

            a.plot(xlin, urinit) 
            a.plot(xlin, ur)
            a.set_title("Radial Velocity" )
            count += 1
        elif count == 2:
            uphyinit = (X[:][int(N/2)]*init3.T[int(N/2)]/ init1.T[int(N/2)] - Y[:][int(N/2)]*init3.T[int(N/2)] / init1.T[int(N/2)])/(R[:][int(N/2)])
            uphy = (X[:][int(N/2)]*v3.T[int(N/2)] / v1.T[int(N/2)] - Y[:][int(N/2)]*v3.T[int(N/2)]/ v1.T[int(N/2)])/(R[:][int(N/2)])
            a.plot(xlin, uphyinit)
            a.plot(xlin, uphy)
            a.set_title("Angular Velocity" )
            count += 1
        elif count == 3:
            a.plot(xlin, init4.T[int(N/2)])
            a.plot(xlin, v4.T[int(N/2)])
            a.set_title("$Pixx$" )
            count += 1
        elif count == 4:
            a.plot(xlin, init5.T[int(N/2)])
            a.plot(xlin, v5.T[int(N/2)])
            a.set_title("$Pixy$" )
            count += 1
        elif count == 5:
            a.plot(xlin, init6.T[int(N/2)])
            a.plot(xlin, v6.T[int(N/2)])
            a.set_title("$Piyx$" )
            count += 1
        elif count == 6:
            a.plot(xlin, init7.T[int(N/2)])
            a.plot(xlin, v7.T[int(N/2)])
            a.set_title("$Piyy$" )
            count += 1

    print("Done")
    fig.delaxes(axs[-1,-1])
    fig.suptitle('NonrelativisticShearIS(N={}, gamma = {:.2f}, zeta = {:.2f}, eta = {:.2f}, tau_nu = {:.2f}, theta = {:.2f})'.format(N,gamma,zeta,eta,tau_nu,theta), fontsize='xx-large')
    plt.tight_layout()
    plt.savefig(nameOfFigure + ".png")
    plt.show()

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
# this is an example of how to call this function above. It will show the variables in 2d at frame 50 selected

plot_csv_static('KT_ISshear\C++\density_solution.csv', 'KT_ISshear\C++\momentx_solution.csv', 
'KT_ISshear\C++\momenty_solution.csv', 'KT_ISshear\C++\Pixx_solution.csv', 'KT_ISshear\C++\Pixy_solution.csv',
 'KT_ISshear\C++\Piyx_solution.csv', 'KT_ISshear\C++\Piyy_solution.csv',
  'KT_ISshear\C++\parameters.csv', 50 ,"NonRelativisticIS")
'''


''' Methods for reading and ploting for python code '''


def show_2dsolution_static(file,parameters_file,i,n):

    sol = np.load(file)

    t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt(parameters_file,delimiter=',',unpack=True)

    N = int(N)
    a = sol[i][n*N:(n+1)*N].T
    
    plt.imshow(a)
    plt.show()
    

def show_solution_static_slice(file,parameters_file,i,n,title,label):

    sol = np.load(file)

    t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt(parameters_file,delimiter=',',unpack=True)

    N = int(N)
    a = sol[i][n*N:(n+1)*N].T
    
    plt.plot(a[int(N/2)], label=label)
    plt.show()

    plt.savefig("static_slice_{}".format(file))


def animate_numpy_solution_gif_after_loading(sol,N):


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

    



def animate_numpy_solution_slice_gif_after_loading_solution(sol,N,xlin):


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

    
                            




def gif2d_python_numpyfile(filename,parameter_filename,frames, name_of_figure):

    #names should be strings

    #read the solution and parameters from the especified files
    sol = np.load(filename + ".npy")
    t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt(parameter_filename + ".txt",delimiter=',',unpack=True)
    N = int(N)
    print(sol.shape)

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes()
    #line, = ax.plot([], [], lw=2)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    im=plt.imshow(sol[0][:N].T,interpolation='none')

    # initialization function: plot the background of each frame
    def init():
        im.set_data(sol[0][:N].T)
        return [im]

    # animation function.  This is called sequentially
    def animate(i):
        im.set_array(sol[i][:N].T)
        return [im]

    # create animation
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=frames, interval=20, blit=True)

    anim.save(name_of_figure + ".gif", fps=60)

def gif2d_python_numpyfile(filename,parameter_filename,frames, name_of_figure):

    #names should be strings 

    #read the solution and parameters from the especified files
    sol = np.load(filename)
    t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt(parameter_filename,delimiter=',',unpack=True)
    N = int(N)
    print(sol.shape)
    xlin = np.linspace(float(a),float(b),N)

    # First set up the figure, the axis, and the plot element we want to animate
    fig_slice = plt.figure()
    ax_slice = plt.axes()
    ax_slice.set_xlim((float(a),float(b)))
    ax_slice.set_ylim((float(a),float(b)))
    line, = ax_slice.plot([], [], lw=2)

    # initialization function: plot the background of each frame
    def init_slice():
        line.set_data(xlin,sol[0][int(N/2)].T)
        
        return (line,)

    # animation function.  This is called sequentially
    def animate_slice(i):
        line.set_data(xlin,sol[i][int(N/2)].T)
        return (line,)

    anim_slice = animation.FuncAnimation(fig_slice, animate_slice, init_func=init_slice,
                                frames=frames, interval=20, blit=True)

    # save animation as name_of_figure_slice.gif
    anim_slice.save(name_of_figure+'_slice.gif', fps=60)


def plot_all_variables_static(filename,parameter_filename, name_of_figure):

    # create figure and axis with 7 subplots for all 7 variables
    fig, axs = plt.subplots(2, 4, figsize=(10, 6))
    gridspec = axs[0, 0].get_subplotspec().get_gridspec()
    print("created the figure")

    #read the solution and parameters from the specified files
    sol = np.load(filename + ".npy")
    t,tEnd,tOut,N,boxsize,gamma,zeta,eta,tau_nu,theta,a,b = np.loadtxt(parameter_filename + ".txt",delimiter=',',unpack=True)
    N = int(N)

    # plot data in remaining axes:
    i = 100
    count = 0
    print("plotting subfigures")
    for a in axs[:, :].flat:
        if count == 0:
            a.imshow(sol[i][0:N].T)
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
            a.set_title("Density")
            count += 1
        elif count == 1:
            ux  = (sol[i][N:2*N].T/ sol[-1][0:N].T) 
            a.imshow(ux)
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
            a.set_title("x Velocity" )
            count += 1
        elif count == 2:
            uy  = (sol[i][2*N:3*N].T / sol[-1][0:N].T) 
            a.imshow(uy)
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
            a.set_title("y Velocity" )
            count += 1
        elif count == 3:
            a.imshow(sol[i][3*N:4*N].T)
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
            a.set_title("$Pixx$" )
            count += 1
        elif count == 4:
            a.imshow(sol[i][4*N:5*N].T)
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
            a.set_title("$Pixy$" )
            count += 1
        elif count == 5:
            a.imshow(sol[i][5*N:6*N].T)
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
            a.set_title("$Piyx$" )
            count += 1
        elif count == 6:
            a.imshow(sol[i][6*N:7*N].T)
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
            a.set_title("$Piyy$" )
            count += 1

    print("Done")
    fig.delaxes(axs[-1,-1])
    fig.suptitle('NonrelativisticShearIS(N={}, gamma = {:.2f}, zeta = {:.2f}, eta = {:.2f}, tau_nu = {:.2f}, theta = {:.2f})'.format(N,gamma,zeta,eta,tau_nu,theta), fontsize='xx-large')
    plt.tight_layout()
    plt.savefig(name_of_figure +'_t:' + str(i*tOut) + ".png")
