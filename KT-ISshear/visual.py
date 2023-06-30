#Visualization code 

import pandas as pd
import matplotlib.pyplot as pt
import numpy as np

def plot_ic_csv(file):

    N = 128
    boxsize = 1
    dx = boxsize/N
    xlin = np.linspace(0.5*dx, boxsize-0.5*dx,N)

    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    s = X.shape
    R = np.sqrt(X**2 + Y**2)

    rho = np.zeros(s)
    Momx = np.zeros(s)
    Momy = np.zeros(s)  
    Pixx = np.zeros(s)
    Pixy = np.zeros(s)
    Piyx = np.zeros(s)
    Piyy = np.zeros(s)

    IC = pd.read_csv(file, header=None)

    # get initial conditions
    for i in range(N):
        for j in range(N):
            rho[i][j]  = float(IC[7*(N*i + j)])
            Momx[i][j] = float(IC[7*(N*i + j) + 1])
            Momy[i][j] = float(IC[7*(N*i + j) + 2])
            Pixx[i][j] = float(IC[7*(N*i + j) + 3])
            Pixy[i][j] = float(IC[7*(N*i + j) + 4])
            Piyx[i][j] = float(IC[7*(N*i + j) + 5])
            Piyy[i][j] = float(IC[7*(N*i + j) + 6])
    
    plt.imshow(rho.T)
    plt.show()



def plot_solution_csv(file):

    N = 128
    boxsize = 1
    dx = boxsize/N
    xlin = np.linspace(0.5*dx, boxsize-0.5*dx,N)

    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    s = X.shape
    R = np.sqrt(X**2 + Y**2)

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




def plot_each_csv(file):
    
    N = 200
    boxsize = 1
    dx = boxsize/N
    xlin = np.linspace(0.5*dx, boxsize-0.5*dx,N)

    Y, X = np.meshgrid( xlin, xlin ) # define the mesh grid
    S = X.shape
    R = np.sqrt(X**2 + Y**2)

    v = np.zeros(S)
    

    print("starting...")
    solution = pd.read_csv(file, header=None)
    print("solution Dataframe created")

    sol = solution.to_numpy

    s = 0
    outputcount = 0
    print("reading solution...")
    while s < sol.shape[1]:
        for i in range(N):
            for j in range(N):
                v[i][j]  = float(sol[(N*i + j)][s])


        plt.imshow(v.T)
        plt.show()
        outputcount += 1
        print(outputcount)
        s +=1 


plot_each_csv('./density_solution.csv')