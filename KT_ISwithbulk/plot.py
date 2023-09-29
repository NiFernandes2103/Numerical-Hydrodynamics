import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def plot_csv_static(file1, file2, file3, file4, parameters_file, s, nameOfFigure, hide_labels = False):
    
    parameters = pd.read_csv(parameters_file)
    N = int(parameters['N'])
    boxsize = int(parameters['boxsize'])
    a = float(parameters['a'])
    b = float(parameters['b'])
    dx = boxsize/N
    xlin = np.linspace(a,b,N)

    gamma = float(parameters['gamma'])
    zeta = float(parameters['zeta'])
    tau_nu = float(parameters['tau_nu'])
    theta = float(parameters['theta'])

    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    R = np.sqrt(X*X+Y*Y)
    S = X.shape

    init1 = np.zeros(S)
    v1    = np.zeros(S)
    init2 = np.zeros(S)
    v2    = np.zeros(S)
    init3 = np.zeros(S)
    v3    = np.zeros(S)
    init4 = np.zeros(S)
    v4    = np.zeros(S)


    print("starting...")
    solution1 = pd.read_csv(file1, header=None)
    solution2 = pd.read_csv(file2, header=None)
    solution3 = pd.read_csv(file3, header=None)
    solution4 = pd.read_csv(file4, header=None)
  
    print("solution Dataframe created")

    sol1 = solution1.to_numpy()
    sol2 = solution2.to_numpy()
    sol3 = solution3.to_numpy()
    sol4 = solution4.to_numpy()
  

    print("solution format:{}".format(sol1.shape))

    print("reading solution...")
    for i in range(N):
        for j in range(N):
            init1[i][j] = float(sol1[0][N*i + j])
            v1[i][j]  = float(sol1[s][N*i + j])
            init2[i][j] = float(sol2[0][N*i + j])
            v2[i][j]  = float(sol2[s][N*i + j])
            init3[i][j] = float(sol3[0][N*i + j])
            v3[i][j]  = float(sol3[s][N*i + j])
            init4[i][j] = float(sol4[0][N*i + j])
            v4[i][j]  = float(sol4[s][N*i + j])
           
            
            
    print('finished reading')

    fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    gridspec = axs[0, 0].get_subplotspec().get_gridspec()
    print("created the figure")

    # clear the left column for the subfigure:

    # plot data in remaining axes:
    count = 0
    print("plotting subfigures")
    for a in axs[:, :].flat:
        if count == 0:
            a.plot(xlin, init1[int(N/2)]) 
            a.plot(xlin, v1[int(N/2)])
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
            a.plot(xlin, init4[int(N/2)])
            a.plot(xlin, v4[int(N/2)])
            a.set_title("Pi")
            count += 1
    
    print("Done")
    fig.suptitle('NonrelativisticIS(N = {}, gamma = {:.2f}, zeta ={:.2f}, tau_nu = {:.2f}, theta = {:.2f})'.format(N,gamma,zeta,tau_nu,theta), fontsize='xx-large')
    plt.tight_layout()
    plt.savefig(nameOfFigure + ".png")
    plt.show()


plot_csv_static('KT_ISwithbulk\C++\density_solution.csv', 'KT_ISwithbulk\C++\momentx_solution.csv', 
'KT_ISwithbulk\C++\momenty_solution.csv', 'KT_ISwithbulk\C++\Pi_solution.csv',
  'KT_ISwithbulk\C++\parameters.csv', 5 ,"NonRelativisticIS")