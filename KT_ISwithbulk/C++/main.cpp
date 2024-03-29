
#include <cmath>
#include <vector>
#include <tuple>
#include <list>
#include <iostream>
#include <string>
#include "file_handler.h"
#include "file_handler.cpp"
#include "KTmethods2dBulk.h"
#include "KTmethods2dBulk.cpp"
#include "nonRelativisticISBulk.h"
#include "nonRelativisticISBulk.cpp"


int main() {

    double t = 0.0;  // s 
    double tEnd = 1;  // time at the end
    double tOut = 0.01;  // time of each output

    int N = 200;  // resolution
    double boxsize = 4.0;  // in some unit system l
    double gamma = 2;  // adiabatic index
    double zeta = 1.0;  // bulk viscosity coefficient
    double tau_nu = 1.0;  // relaxation time
    double theta = 1.0;  // flux limiter parameter
    double dx = boxsize / N;  // box size
    double dy = dx;
    //double vol = dx * dx;  // volume of each box
    double a = 0.5*(0.5 * dx - boxsize); 
    double b = 0.5*(boxsize - 0.5 * dx);

    parameters_csv(t,tEnd,tOut,N,boxsize,a,b,gamma,zeta,tau_nu,theta,"parameters.csv");

    std::vector<double> xlin(N);
    for (int i = 0; i < N; i++) {
        //xlin[i] = 0.5 * dx + (boxsize - 0.5 * dx) * i / (N - 1);  // simulation limits

        xlin[i] = a +  (b-a)* i / (N - 1);
    }


    std::vector<std::vector<double>> Y(N, std::vector<double>(N, 0.0));
    std::vector<std::vector<double>> X(N, std::vector<double>(N, 0.0));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Y[i][j] = xlin[j];
            X[i][j] = xlin[i];
        }
    }
    int s = X.size();
    std::vector<std::vector<double>> R(s, std::vector<double>(s, 0.0));
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            R[i][j] = sqrt(X[i][j] * X[i][j] + Y[i][j] * Y[i][j]);
        }
    }

    //double w0 = 0.1;
    //double sigma = 0.05 / sqrt(2.0);
    std::vector<std::vector<double>> rho(s, std::vector<double>(s, 0.0));
    std::vector<std::vector<double>> vx(s, std::vector<double>(s, 0.0));
    std::vector<std::vector<double>> vy(s, std::vector<double>(s, 0.0));
    std::vector<std::vector<double>> Momx(s, std::vector<double>(s, 0.0));
    std::vector<std::vector<double>> Momy(s, std::vector<double>(s, 0.0));
    std::vector<std::vector<double>> Pi(s, std::vector<double>(s, 0.0));
    
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            

            rho[i][j] = 2*pow((1 - (R[i][j])*(R[i][j]) / 0.25) ,4)*(R[i][j] < 0.5) + 2*pow((1 - (R[i][j])*(R[i][j]) / 0.25) ,4)*(R[i][j] < 0.5) + 1; // Mauricio's function advice
            
        }
    }

    stateb IC = {rho, Momx, Momy, Pi};

    create(IC, "initial_state.csv");

    std::list<stateb> initial_state = {IC};

std::list<stateb> solution = integrator(KTschemeNonRelativisticIS, std::make_tuple(t, tEnd), initial_state, tOut, std::make_tuple(dx, dy, N, gamma, zeta, tau_nu, theta), "Heuns");

    //write(solution, "NonrelativisticISwithsmoothIC.csv");
    
    write_each(solution, "density_solution.csv", 0);
    write_each(solution, "momentx_solution.csv", 1);
    write_each(solution, "momenty_solution.csv", 2);
    write_each(solution, "Pi_solution.csv", 3);

    
}

