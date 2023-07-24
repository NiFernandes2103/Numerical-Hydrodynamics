#include <cmath>
#include <tuple>
#include <list>
#include <iostream>
#include "fileFunc.h"
#include "fileFunc.cpp"
#include "KTmethods2d.h"
#include "KTmethods2d.cpp"
#include "nonRelativisticISwithShear.h"
#include "nonRelativisticISwithShear.cpp"

using namespace std;

int main() {

    double t = 0.0;  // s 
    double tEnd = 0.5;  // time at the end
    double tOut = 0.01;  // time of each output

    int N = 200;  // resolution
    double boxsize = 1.0;  // in some unit system l
    double gamma = 5/3;  // adiabatic index
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

    double xlin[N];
    for (int i = 0; i < N; i++) {
        //xlin[i] = 0.5 * dx + (boxsize - 0.5 * dx) * i / (N - 1);  // simulation limits

        xlin[i] = a +  (b-a)* i / (N - 1);
    }


    smatrix Y(N);
    smatrix X(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Y.set(xlin[j],i,j);
            X.set(xlin[i],i,j);

            
        }
    }

    smatrix R(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R.set(sqrt(X.get(i,j) * X.get(i,j) + Y.get(i,j) * Y.get(i,j)), i,j);
        }
    }

    //double w0 = 0.1;
    //double sigma = 0.05 / sqrt(2.0);
    smatrix rho(N);
    smatrix vx(N);
    smatrix vy(N);
    smatrix Momx(N);
    smatrix Momy(N);
    smatrix Pixx(N);
    smatrix Pixy(N);
    smatrix Piyx(N);
    smatrix Piyy(N);


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
           
            rho.set((pow((1 - (R.get(i,j))*(R.get(i,j))),4))*(R.get(i,j) < 1) + 1, i,j); // Mauricio's function advice

        }
    }

    state IC = {rho, Momx, Momy, Pixx, Pixy, Piyx, Piyy};

    create(IC, "initial_state.csv");

    list<state> initial_state = {IC};

list<state> solution = integrator(KTschemeNonRelativisticIS, make_tuple(t, tEnd), initial_state, tOut, make_tuple(dx, dy, N, gamma, zeta, tau_nu, eta, theta), "Heuns");

    write_each(solution, "density_solution.csv", 0);
    write_each(solution, "momentx_solution.csv", 1);
    write_each(solution, "momenty_solution.csv", 2);
    write_each(solution, "Pixx_solution.csv", 3);
    write_each(solution, "Pixy_solution.csv", 4);
    write_each(solution, "Piyx_solution.csv", 5);
    write_each(solution, "Piyy_solution.csv", 6);
    
}

