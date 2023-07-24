
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
    double tEnd = 1;  // time at the end
    double tOut = 0.01;  // time of each output

    int N = 200;  // resolution
    double boxsize = 4.0;  // in some unit system l
    double gamma = 1;  // adiabatic index
    double zeta = 1.0;  // bulk viscosity coefficient
    double eta = 1.0;  // shear viscosity coefficient
    double tau_nu = 1.0;  // relaxation time
    double theta = 1.0;  // flux limiter parameter
    double dx = boxsize / N;  // box size
    double dy = dx;
    //double vol = dx * dx;  // volume of each box
    double a = 0.5*(0.5 * dx - boxsize); 
    double b = 0.5*(boxsize - 0.5 * dx);

    parameters_csv(t,tEnd,tOut,N,boxsize,a,b,gamma,zeta,eta,tau_nu,theta,"parameters.csv");

    vector<double> xlin(N);
    for (int i = 0; i < N; i++) {
        //xlin[i] = 0.5 * dx + (boxsize - 0.5 * dx) * i / (N - 1);  // simulation limits

        xlin[i] = a +  (b-a)* i / (N - 1);
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

    //double w0 = 0.1;
    //double sigma = 0.05 / sqrt(2.0);
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
            //rho[i][j] = 1.0 + (abs(Y[i][j] - 0.5) < 0.25);
            //vx[i][j] = -0.5 + (abs(Y[i][j] - 0.5) < 0.25);
            //vy[i][j] = w0 * sin(4 * M_PI * X[i][j]) * (exp(-(Y[i][j] - 0.25) * (Y[i][j] - 0.25) / (2 * sigma * sigma)) + exp(-(Y[i][j] - 0.75) * (Y[i][j] - 0.75) / (2 * sigma * sigma)));
            //Momx[i][j] = vx[i][j]*rho[i][j];
            //Momy[i][j] = vy[i][j]*rho[i][j];

            rho[i][j] = (pow((1 - (R[i][j])*(R[i][j])),4))*(R[i][j] < 1) + 1; // Mauricio's function advice

        }
    }

    state IC = {rho, Momx, Momy, Pixx, Pixy, Piyx, Piyy};

    create(IC, "initial_state.csv");

    map<double, state> initial_state = {{t, IC}};

map<double, state> solution = integrator(KTschemeNonRelativisticIS, make_tuple(t, tEnd), initial_state, tOut, make_tuple(dx, dy, N, gamma, zeta, tau_nu, eta, theta), "Heuns");

    write_each(solution, "density_solution.csv", 0);
    write_each(solution, "momentx_solution.csv", 1);
    write_each(solution, "momenty_solution.csv", 2);
    write_each(solution, "Pixx_solution.csv", 3);
    write_each(solution, "Pixy_solution.csv", 4);
    write_each(solution, "Piyx_solution.csv", 5);
    write_each(solution, "Piyy_solution.csv", 6);
    
}

