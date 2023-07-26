#include <cmath>
#include <vector>
#include <tuple>
#include <list>
#include <iostream>
#include "KTmethods2d.h"
#include "nonRelativisticISwithShear.h"
using namespace std;

state KTschemeNonRelativisticIS(double t,  state& IC, double dx, double dy, int N, double gamma, double zeta, double tau_nu, double eta, double theta = 1) {
     
    /* Finite Volume simulation */

    // Generate Initial Conditions

    /* Initial conditions for rho */ 
    double** rho;
    double** P = new double*[N];

    /* Initial conditions for v */
    double** Momx;
    double** Momy;
    double** vx = new double*[N];
    double** vy = new double*[N];

    /* Pi initial condition */
    double** Pixx;
    double** Pixy;
    double** Piyx;
    double** Piyy;

    rho  = IC.get(0);
    Momx = IC.get(1);
    Momy = IC.get(2);
    Pixx = IC.get(3);
    Pixy = IC.get(4);
    Piyx = IC.get(5);
    Piyy = IC.get(6);
    

    for (int i = 0; i < N; i++) {
        vx[i] = new double[N];
        vy[i] = new double[N];
        P[i]  = new double[N];
        for (int j = 0; j < N; j++) {

            // set boundary conditions
            rho[0][j] = rho[1][j];
            rho[i][0] = rho[i][1];
            rho[N-1][j] = rho[N-2][j];
            rho[i][N-1] = rho[i][N-2];

            Momx[0][j] = -Momx[1][j];
            Momx[i][0] = 0;
            Momx[N-1][j] = -Momx[N-2][j];
            Momx[i][N-1] = 0;

            Momy[0][j] = 0;
            Momy[i][0] = -Momy[i][1];
            Momy[N-1][j] = 0;
            Momy[i][N-1] = -Momy[i][N-2];

            Pixx[0][j] = Pixx[1][j];
            Pixx[i][0] = Pixx[i][1];
            Pixx[N-1][j] = Pixx[N-2][j];
            Pixx[i][N-1] = Pixx[i][N-2];

            Pixy[0][j] = Pixy[1][j];
            Pixy[i][0] = Pixy[i][1];
            Pixy[N-1][j] = Pixy[N-2][j];
            Pixy[i][N-1] = Pixy[i][N-2];
            
            Piyx[0][j] = Piyx[1][j];
            Piyx[i][0] = Piyx[i][1];
            Piyx[N-1][j] = Piyx[N-2][j];
            Piyx[i][N-1] = Piyx[i][N-2];

            Piyy[0][j] = Piyy[1][j];
            Piyy[i][0] = Piyy[i][1];
            Piyy[N-1][j] = Piyy[N-2][j];
            Piyy[i][N-1] = Piyy[i][N-2];

            vx[i][j] = Momx[i][j] / rho[i][j];
            vy[i][j] = Momy[i][j] / rho[i][j];

             /* Pressure from equation of state */
            P[i][j] = pow(rho[i][j], gamma);
        }
    }        



    // calculate speed of sound
    // getSpeedOfSound(rho, gamma)
    //double** cs = getSpeedOfSound(rho, gamma, N, N);

    // calculate gradients
    double** rho_dx = getGradient(rho, dx, 0, N, N, theta);
    double** vx_dx = getGradient(vx, dx, 0, N, N, theta);
    double** vy_dx = getGradient(vy, dx, 0, N, N, theta);
    double** P_dx = getGradient(P, dx, 0, N, N, theta);
    double** Pixx_dx = getGradient(Pixx, dx, 0, N, N, theta);
    double** Pixy_dx = getGradient(Pixy, dx, 0, N, N, theta);
    double** Piyx_dx = getGradient(Piyx, dx, 0, N, N, theta);
    double** Piyy_dx = getGradient(Piyy, dx, 0, N, N, theta);

    double** rho_dy = getGradient(rho, dy, 1, N, N, theta);
    double** vx_dy = getGradient(vx, dy, 1, N, N, theta);
    double** vy_dy = getGradient(vy, dy, 1, N, N, theta);
    double** P_dy = getGradient(P, dy, 1, N, N, theta);
    double** Pixx_dy = getGradient(Pixx, dy, 1, N, N, theta);
    double** Pixy_dy = getGradient(Pixy, dy, 1, N, N, theta);
    double** Piyx_dy = getGradient(Piyx, dy, 1, N, N, theta);
    double** Piyy_dy = getGradient(Piyy, dy, 1, N, N, theta);
    

    // extrapolate fluxes

    double** rhoM_XL, **rhoP_XL, **rhoM_XR, **rhoP_XR;
    std::tie(rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR) = extrapolateInSpaceToFace(rho, rho_dx, dx, 0, N, N);
    double** vxM_XL, **vxP_XL,  **vxM_XR,  **vxP_XR;   
    std::tie(vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR)  = extrapolateInSpaceToFace(vx,  vx_dx,   dx, 0, N, N);
    double** vyM_XL, **vyP_XL,  **vyM_XR,  **vyP_XR;   
    std::tie(vyM_XL,  vyP_XL,  vyM_XR,  vyP_XR)  = extrapolateInSpaceToFace(vy,  vy_dx,   dx, 0, N, N);
    double** PM_XL,   **PP_XL,   **PM_XR,   **PP_XR;   
    std::tie(PM_XL,   PP_XL,   PM_XR,   PP_XR)   = extrapolateInSpaceToFace(P,   P_dx,   dx, 0, N, N);
    double** PixxM_XL,  **PixxP_XL,  **PixxM_XR,  **PixxP_XR;   
    std::tie(PixxM_XL,  PixxP_XL,  PixxM_XR,  PixxP_XR)  = extrapolateInSpaceToFace(Pixx,  Pixx_dx,  dx, 0, N, N);
    double** PixyM_XL,  **PixyP_XL,  **PixyM_XR,  **PixyP_XR;   
    std::tie(PixyM_XL,  PixyP_XL,  PixyM_XR,  PixyP_XR)  = extrapolateInSpaceToFace(Pixy,  Pixy_dx,  dx, 0, N, N);
    double** PiyxM_XL,  **PiyxP_XL,  **PiyxM_XR,  **PiyxP_XR;   
    std::tie(PiyxM_XL,  PiyxP_XL,  PiyxM_XR,  PiyxP_XR)  = extrapolateInSpaceToFace(Piyx,  Piyx_dx,  dx, 0, N, N);
    double** PiyyM_XL,  **PiyyP_XL,  **PiyyM_XR,  **PiyyP_XR;   
    std::tie(PiyyM_XL,  PiyyP_XL,  PiyyM_XR,  PiyyP_XR)  = extrapolateInSpaceToFace(Piyy,  Piyy_dx,  dx, 0, N, N);

    double** rhoM_YL, **rhoP_YL, **rhoM_YR, **rhoP_YR;
    std::tie(rhoM_YL, rhoP_YL, rhoM_YR, rhoP_YR) = extrapolateInSpaceToFace(rho, rho_dy, dy, 1, N, N);
    double** vxM_YL,  **vxP_YL,  **vxM_YR,  **vxP_YR;   
    std::tie(vxM_YL,  vxP_YL,  vxM_YR,  vxP_YR)  = extrapolateInSpaceToFace(vx,  vx_dy,  dy, 1, N, N);
    double** vyM_YL,  **vyP_YL,  **vyM_YR,  **vyP_YR;   
    std::tie(vyM_YL,  vyP_YL,  vyM_YR,  vyP_YR)  = extrapolateInSpaceToFace(vy,  vy_dy,  dy, 1, N, N);
    double** PM_YL,   **PP_YL,   **PM_YR,   **PP_YR;   
    std::tie(PM_YL,   PP_YL,   PM_YR,   PP_YR)   = extrapolateInSpaceToFace(P,   P_dy,   dy, 1, N, N);
    double** PixxM_YL,  **PixxP_YL,  **PixxM_YR,  **PixxP_YR;   
    std::tie(PixxM_YL,  PixxP_YL,  PixxM_YR,  PixxP_YR)  = extrapolateInSpaceToFace(Pixx,  Pixx_dy,  dy, 1, N, N);
    double** PixyM_YL,  **PixyP_YL,  **PixyM_YR,  **PixyP_YR;   
    std::tie(PixyM_YL,  PixyP_YL,  PixyM_YR,  PixyP_YR)  = extrapolateInSpaceToFace(Pixy,  Pixy_dy,  dy, 1, N, N);
    double** PiyxM_YL,  **PiyxP_YL,  **PiyxM_YR,  **PiyxP_YR;   
    std::tie(PiyxM_YL,  PiyxP_YL,  PiyxM_YR,  PiyxP_YR)  = extrapolateInSpaceToFace(Piyx,  Piyx_dy,  dy, 1, N, N);
    double** PiyyM_YL,  **PiyyP_YL,  **PiyyM_YR,  **PiyyP_YR;   
    std::tie(PiyyM_YL,  PiyyP_YL,  PiyyM_YR,  PiyyP_YR)  = extrapolateInSpaceToFace(Piyy,  Piyy_dy,  dy, 1, N, N);

    // Compute fluxes

    double** flux_Mass_XR, **flux_Momx_XR, **flux_Momy_XR, **flux_Pixx_vxR, **flux_Pixy_vxR, **flux_Piyx_vxR, **flux_Piyy_vxR;
    std::tie(flux_Mass_XR, flux_Momx_XR, flux_Momy_XR, flux_Pixx_vxR, flux_Pixy_vxR, flux_Piyx_vxR, flux_Piyy_vxR) = getXFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR,
                                            vyP_XR, vyM_XR, PixxP_XR, PixxM_XR,
                                            PixyP_XR, PixyM_XR, PiyxP_XR, PiyxM_XR,
                                            PiyyP_XR, PiyyM_XR, PP_XR, PM_XR, gamma,
                                            eta, zeta, tau_nu, N, N);

    double** flux_Mass_XL, **flux_Momx_XL, **flux_Momy_XL, **flux_Pixx_vxL, **flux_Pixy_vxL, **flux_Piyx_vxL, **flux_Piyy_vxL;
    std::tie(flux_Mass_XL, flux_Momx_XL, flux_Momy_XL, flux_Pixx_vxL, flux_Pixy_vxL,flux_Piyx_vxL, flux_Piyy_vxL) = getXFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL,
                                            vyP_XL, vyM_XL, PixxP_XL, PixxM_XL,
                                            PixyP_XL, PixyM_XL, PiyxP_XL, PiyxM_XL,
                                            PiyyP_XL, PiyyM_XL, PP_XL, PM_XL, gamma,
                                            eta, zeta, tau_nu, N, N);

    double** flux_Mass_YR, **flux_Momx_YR, **flux_Momy_YR, **flux_Pixx_vyR, **flux_Pixy_vyR, **flux_Piyx_vyR, **flux_Piyy_vyR;
    std::tie(flux_Mass_YR, flux_Momx_YR, flux_Momy_YR, flux_Pixx_vyR, flux_Pixy_vyR, flux_Piyx_vyR, flux_Piyy_vyR) = getYFlux(rhoP_YR, rhoM_YR, vxP_YR, vxM_YR,
                                            vyP_YR, vyM_YR, PixxP_YR, PixxM_YR,
                                            PixyP_YR, PixyM_YR, PiyxP_YR, PiyxM_YR,
                                            PiyyP_YR, PiyyM_YR, PP_YR, PM_YR, gamma,
                                            eta, zeta, tau_nu, N, N);

    double** flux_Mass_YL, **flux_Momx_YL, **flux_Momy_YL, **flux_Pixx_vyL, **flux_Pixy_vyL, **flux_Piyx_vyL, **flux_Piyy_vyL;
    std::tie(flux_Mass_YL, flux_Momx_YL, flux_Momy_YL, flux_Pixx_vyL, flux_Pixy_vyL, flux_Piyx_vyL, flux_Piyy_vyL) = getYFlux(rhoP_YL, rhoM_YL, vxP_YL, vxM_YL,
                                            vyP_YL, vyM_YL, PixxP_YL, PixxM_YL,
                                            PixyP_YL, PixyM_YL, PiyxP_YL, PiyxM_YL,
                                            PiyyP_YL, PiyyM_YL, PP_YL, PM_YL, gamma,
                                            eta, zeta, tau_nu, N, N);
    
    
    
    // Update conservative variables

    double** J = new double*[N];
    double** Jxx = new double*[N];
    double** Jxy = new double*[N];
    double** Jyx = new double*[N];
    double** Jyy = new double*[N];

    for (int i = 0; i < N; i++) {
        J[i] = new double[N];
        Jxx[i] = new double[N];
        Jxy[i] = new double[N];
        Jyx[i] = new double[N];
        Jyy[i] = new double[N];        
        for (int j = 0; j < N; j++) {
            Jxx[i][j] = -Pixx[i][j]/tau_nu;
            Jxy[i][j] = -Pixy[i][j]/tau_nu;
            Jyx[i][j] = -Piyx[i][j]/tau_nu;
            Jyy[i][j] = -Piyy[i][j]/tau_nu;
        }
    }

    double** timederivative_rho, **timederivative_Momx, **timederivative_Momy;
    double** timederivative_Pixx, **timederivative_Pixy, **timederivative_Piyx, **timederivative_Piyy;

    timederivative_rho   = applyFluxes( flux_Mass_XR,   flux_Mass_XL, flux_Mass_YR,   flux_Mass_YL,  dx, dy, J, N, N);
    timederivative_Momx  = applyFluxes( flux_Momx_XR,   flux_Momx_XL, flux_Momx_YR,   flux_Momx_YL, dx, dy, J, N, N); 
    timederivative_Momy  = applyFluxes( flux_Momy_XR,   flux_Momy_XL, flux_Momy_YR,   flux_Momy_YL, dx, dy, J, N, N);
    timederivative_Pixx  = applyFluxes( flux_Pixx_vxR,    flux_Pixx_vxL, flux_Pixx_vyR,    flux_Pixx_vyL, dx, dy, Jxx, N, N);
    timederivative_Pixy  = applyFluxes( flux_Pixy_vxR,    flux_Pixy_vxL, flux_Pixy_vyR,    flux_Pixy_vyL, dx, dy, Jxy, N, N);
    timederivative_Piyx  = applyFluxes( flux_Piyx_vxR,    flux_Piyx_vxL, flux_Piyx_vyR,    flux_Piyx_vyL, dx, dy, Jyx, N, N);
    timederivative_Piyy  = applyFluxes( flux_Piyy_vxR,    flux_Piyy_vxL, flux_Piyy_vyR,    flux_Piyy_vyL, dx, dy, Jyy, N, N);
    
    for (int i = 0; i < N; i++) {
        delete[] vx[i];
        delete[] vy[i];
        delete[] P[i];
        delete[] Jxx[i];
        delete[] Jxy[i];
        delete[] Jyx[i];
        delete[] Jyy[i];
    }

    delete[] rho;
    delete[] P;
    delete[] Momx;
    delete[] Momy;
    delete[] vx;
    delete[] vy;
    delete[] Pixx;
    delete[] Pixy;
    delete[] Piyx;
    delete[] Piyy;
    delete[] J;
    delete[] Jxx;
    delete[] Jxy;
    delete[] Jyx;
    delete[] Jyy;

    return {timederivative_rho,timederivative_Momx,timederivative_Momy,timederivative_Pixx,timederivative_Pixy,timederivative_Piyx,timederivative_Piyy};
    
}


list<state> integrator(state (*scheme)(double, state&, double, double, int, double, double, double, double, double), tuple<double,double> time, list<state> Q, double dtmax,  tuple<double, double, int, double, double, double, double, double> args, string method)
{
    /*
    This is an integrator that evolves a

    scheme     is the method to get dy/dt e.g. KTscheme
    time       is the current time
    q0         is the current state
    dtmax      is the upperbound of dt set by the user
    BC         is a function that enforces the boundary conditions
    method     is the method used in the integrator
    args       are additional arguments for scheme
    */

    cout.precision(2);

    double t = get<0>(time);
    double tEnd = get<1>(time);
    int outCount = 1;

    state q = Q.back();
    state qprime;

    int N;
    double dx, dy, gamma, zeta, tau_nu, eta, theta;
    std::tie(dx, dy, N, gamma, zeta, tau_nu, eta, theta) = args;
    


    double** rho = new double*[N];
    double** vx = new double*[N];
    double** vy = new double*[N];
    double** Momx = new double*[N];
    double** Momy = new double*[N];
    double** cs = new double*[N];

    auto C = [&](double t, state y) {return scheme(t, q, dx, dy, N, gamma, zeta, tau_nu, eta, theta);};

    while (t < tEnd) {

        rho = q.get(0);
        Momx  = q.get(1);
        Momy  = q.get(2);
        cs  = getSpeedOfSound(rho, gamma, N, N);

        for (int i = 0; i < N; i++) {
            vx[i] = new double[N];
            vy[i] = new double[N];
            for (int j = 0; j < N; j++) {
                vx[i][j] = Momx[i][j] / rho[i][j];
                vy[i][j] = Momy[i][j] / rho[i][j];

        }
    }        

        // condition to ensure that the time steps are small enough so that
        // waves do not interfere with each other
        double** propagation_speed = local_propagation_speed(rho, vx, vy, eta, zeta, tau_nu, cs, N, N);
        double courant_number;

        
        courant_number = dx / max_value(propagation_speed, N, N);
            
        double dt = std::min(dtmax, 0.4 * (courant_number));

        if (method == "Heuns") {
            qprime = heuns(q, C, dt, t, N, N);
        } if (method == "RK4") {
            qprime = rK4(q, C, dt, t, N, N);
        } 
                
        //cout << dt << endl;
        t = t + dt;
        q = qprime;

        if (t > outCount*dtmax) {
            Q.push_back(qprime);
            cout << t << '/' << tEnd << endl;
            ++outCount;
        }
    }

    return Q;
}


