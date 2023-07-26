
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

    double* xlin = new double[N];
    for (int i = 0; i < N; i++) {
        //xlin[i] = 0.5 * dx + (boxsize - 0.5 * dx) * i / (N - 1);  // simulation limits

        xlin[i] = a +  (b-a)* i / (N - 1);
    }


    double** Y = new double*[N];
    double** X = new double*[N];
    double** R = new double*[N];
    for (int i = 0; i < N; i++) {
        Y[i] = new double[N];
        X[i] = new double[N];
        R[i] = new double[N];

        for (int j = 0; j < N; j++) {
            Y[i][j] = xlin[j];
            X[i][j] = xlin[i];
            R[i][j] = sqrt(X[i][j] * X[i][j] + Y[i][j] * Y[i][j]);

        }
    }
    delete[] xlin;
    

    //double w0 = 0.1;
    //double sigma = 0.05 / sqrt(2.0);
    double** rho = new double*[N];
    double** vx = new double*[N];
    double** vy = new double*[N];
    double** Momx = new double*[N];
    double** Momy = new double*[N];
    double** Pixx = new double*[N];
    double** Pixy = new double*[N];
    double** Piyx = new double*[N];
    double** Piyy = new double*[N];

    for (int i = 0; i < N; i++) {
        rho[i]  = new double[N];
        Momx[i] = new double[N];
        Momy[i] = new double[N];
        Pixx[i] = new double[N];
        Pixy[i] = new double[N];
        Piyx[i] = new double[N];
        Piyy[i] = new double[N];
        for (int j = 0; j < N; j++) {
            //rho[i][j] = 1.0 + (abs(Y[i][j] - 0.5) < 0.25);
            //vx[i][j] = -0.5 + (abs(Y[i][j] - 0.5) < 0.25);
            //vy[i][j] = w0 * sin(4 * M_PI * X[i][j]) * (exp(-(Y[i][j] - 0.25) * (Y[i][j] - 0.25) / (2 * sigma * sigma)) + exp(-(Y[i][j] - 0.75) * (Y[i][j] - 0.75) / (2 * sigma * sigma)));
            //Momx[i][j] = vx[i][j]*rho[i][j];
            //Momy[i][j] = vy[i][j]*rho[i][j];

            rho[i][j] = (pow((1 - (R[i][j])*(R[i][j])),4))*(R[i][j] < 1) + 1; // Mauricio's function advice
            Momx[i][j] = 0;
            Momy[i][j] = 0;
            Pixx[i][j] = 0;
            Pixy[i][j] = 0;
            Piyx[i][j] = 0;
            Piyy[i][j] = 0;
        }
    }

    state IC = {rho, Momx, Momy, Pixx, Pixy, Piyx, Piyy};

    create(IC, "initial_state.csv", N, N);

    list<state> initial_state = {IC};

list<state> solution = integrator(KTschemeNonRelativisticIS, make_tuple(t, tEnd), initial_state, tOut, make_tuple(dx, dy, N, gamma, zeta, tau_nu, eta, theta), "Heuns");

    write_each(solution, "density_solution.csv", 0, N, N);
    write_each(solution, "momentx_solution.csv", 1, N, N);
    write_each(solution, "momenty_solution.csv", 2, N, N);
    write_each(solution, "Pixx_solution.csv", 3, N, N);
    write_each(solution, "Pixy_solution.csv", 4, N, N);
    write_each(solution, "Piyx_solution.csv", 5, N, N);
    write_each(solution, "Piyy_solution.csv", 6, N, N);

 for (int i = 0; i < N; i++) {
        delete[] Y[i];
        delete[] X[i];
        delete[] R[i];
    }

    delete[] Y;
    delete[] X;
    delete[] R;

    for (int i = 0; i < N; i++) {
        delete[] rho[i];
        delete[] vx[i];
        delete[] vy[i];
        delete[] Momx[i];
        delete[] Momy[i];
        delete[] Pixx[i];
        delete[] Pixy[i];
        delete[] Piyx[i];
        delete[] Piyy[i];
    }

    delete[] rho;
    delete[] vx;
    delete[] vy;
    delete[] Momx;
    delete[] Momy;
    delete[] Pixx;
    delete[] Pixy;
    delete[] Piyx;
    delete[] Piyy;
      
}

