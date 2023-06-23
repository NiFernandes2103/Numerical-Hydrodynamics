
#include <cmath>
#include <vector>
#include <tuple>
#include <list>
#include <iostream>
#include <map> 
#include "fileFunc.h"
#include "fileFunc.cpp"
#include "KTmethods2d.h"
#include "KTmethods2d.cpp"
#include "nonRelativisticISwithShear.h"
#include "nonRelativisticISwithShear.cpp"


using namespace std;



int main() {

    double t = 0.0;  // s 
    double tEnd = 0.1;  // time at the end
    double tOut = 0.01;  // time of each output

    int N = 128;  // resolution
    double boxsize = 1.0;  // in some unit system l
    double gamma = 5/3;  // adiabatic index
    double zeta = 1.0;  // bulk viscosity coefficient
    double eta = 1.0;  // shear viscosity coefficient
    double tau_nu = 1.0;  // relaxation time
    double theta = 1.0;  // flux limiter parameter

    double dx = boxsize / N;  // box size
    double dy = dx;
    double vol = dx * dx;  // volume of each box
    vector<double> xlin(N);
    for (int i = 0; i < N; i++) {
        xlin[i] = 0.5 * dx + (boxsize - 0.5 * dx) * i / (N - 1);  // simulation limits
    }
    vector<vector<double>> Y(N, vector<double>(N, 0.0));
    vector<vector<double>> X(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Y[i][j] = xlin[j];
            X[i][j] = xlin[i];
        }
    }
    int s = X.size();
    vector<vector<double>> R(s, vector<double>(s, 0.0));
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            R[i][j] = sqrt(X[i][j] * X[i][j] + Y[i][j] * Y[i][j]);
        }
    }

    /* initial condition of density */

    //rho = (1.5*(R <= 0.25) + 1*(R > 0.25))
    //rho = ((1 - ((R - (boxsize-0.5*dx)*0.5)**2)/0.25 )**4 )*(R < 0.5) + 0.1*np.ones(R.shape) // Mauricio`s funtion advice    
    //rho = 1*(X < 0) + 0.125*(X >= 0)

    /* initial condition of velocity */
    //vx = np.zeros(s)
    //vx = 0.5*np.ones(xlin.shape)
    //vx = 3*(Y < 0) - 0*(Y >= 0)
    //vx = np.abs((xlin - (boxsize-0.5*dx)*0.5)/16)

    //vy = np.zeros(s)
    //vy = 0.5*np.ones(xlin.shape)

    double w0 = 0.1;
    double sigma = 0.05 / sqrt(2.0);
    vector<vector<double>> rho(s, vector<double>(s, 0.0));
    vector<vector<double>> vx(s, vector<double>(s, 0.0));
    vector<vector<double>> vy(s, vector<double>(s, 0.0));
    vector<vector<double>> Momx(s, vector<double>(s, 0.0));
    vector<vector<double>> Momy(s, vector<double>(s, 0.0));
    vector<vector<double>> Pixx(s, vector<double>(s, 0.0));
    vector<vector<double>> Pixy(s, vector<double>(s, 0.0));
    vector<vector<double>> Piyx(s, vector<double>(s, 0.0));
    vector<vector<double>> Piyy(s, vector<double>(s, 0.0));

    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            rho[i][j] = 1.0 + (abs(Y[i][j] - 0.5) < 0.25);
            vx[i][j] = -0.5 + (abs(Y[i][j] - 0.5) < 0.25);
            vy[i][j] = w0 * sin(4 * M_PI * X[i][j]) * (exp(-(Y[i][j] - 0.25) * (Y[i][j] - 0.25) / (2 * sigma * sigma)) + exp(-(Y[i][j] - 0.75) * (Y[i][j] - 0.75) / (2 * sigma * sigma)));
            Momx[i][j] = vx[i][j]*rho[i][j];
            Momy[i][j] = vy[i][j]*rho[i][j];
        }
    }

    state IC = {rho, Momx, Momy, Pixx, Pixy, Piyx, Piyy};

    create(IC, "initial_state.csv");

    map<double, state> initial_state = {{t, IC}};

    map<double, state> solution = integrator(KTschemeNonRelativisticIS, make_tuple(t, tEnd), initial_state, tOut, make_tuple(dx, dy, N, gamma, zeta, tau_nu, eta, theta));

    write(solution, "solution.csv");
    
}

