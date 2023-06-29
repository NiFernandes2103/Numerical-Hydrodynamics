
#include <cmath>
#include <vector>
#include <tuple>
#include <list>
#include <iostream>
#include <map> 
#include "KTmethods2d.h"
#include "nonRelativisticISwithShear.h"
using namespace std;

state KTschemeNonRelativisticIS(double t,  state& IC, double dx, double dy, int N, double gamma, double zeta, double tau_nu, double eta, double theta = 1) {
     
    /* Finite Volume simulation */

    // Generate Initial Conditions

    /* Initial conditions for rho */ 
    vector<vector<double>> rho(N,vector<double>(N,0.0));
    vector<vector<double>> P(N,vector<double>(N,0.0));

    /* Initial conditions for v */
    vector<vector<double>> Momx(N, vector<double>(N,0.0));
    vector<vector<double>> Momy(N, vector<double>(N,0.0));
    vector<vector<double>> vx(N, vector<double>(N,0.0));
    vector<vector<double>> vy(N, vector<double>(N,0.0));

    /* Pi initial condition */
    vector<vector<double>> Pixx(N,vector<double>(N,0.0));
    vector<vector<double>> Pixy(N,vector<double>(N,0.0));
    vector<vector<double>> Piyx(N,vector<double>(N,0.0));
    vector<vector<double>> Piyy(N,vector<double>(N,0.0));

    rho  = IC.get(0);
    Momx = IC.get(1);
    Momy = IC.get(2);
    Pixx = IC.get(3);
    Pixy = IC.get(4);
    Piyx = IC.get(5);
    Piyy = IC.get(6);
    
   
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            vx[i][j] = Momx[i][j] / rho[i][j];
            vy[i][j] = Momy[i][j] / rho[i][j];

             /* Pressure from equation of state */
            P[i][j] = pow(fabs(rho[i][j]), gamma);

        }
    }        

    /* B Constant */
    double B = zeta / tau_nu;

     // getSpeedOfSound(rho, gamma)
    vector<vector<double>> cs = getSpeedOfSound(rho, gamma);

    // calculate gradients
    vector<vector<double>> rho_dx = getGradient(rho, dx, 0, theta);
    vector<vector<double>> vx_dx = getGradient(vx, dx, 0, theta);
    vector<vector<double>> vy_dx = getGradient(vy, dx, 0, theta);
    vector<vector<double>> P_dx = getGradient(P, dx, 0, theta);
    vector<vector<double>> Pixx_dx = getGradient(Pixx, dx, 0, theta);
    vector<vector<double>> Pixy_dx = getGradient(Pixy, dx, 0, theta);
    vector<vector<double>> Piyx_dx = getGradient(Piyx, dx, 0, theta);
    vector<vector<double>> Piyy_dx = getGradient(Piyy, dx, 0, theta);

    vector<vector<double>> rho_dy = getGradient(rho, dy, 1, theta);
    vector<vector<double>> vx_dy = getGradient(vx, dy, 1, theta);
    vector<vector<double>> vy_dy = getGradient(vy, dy, 1, theta);
    vector<vector<double>> P_dy = getGradient(P, dy, 1, theta);
    vector<vector<double>> Pixx_dy = getGradient(Pixx, dy, 1, theta);
    vector<vector<double>> Pixy_dy = getGradient(Pixy, dy, 1, theta);
    vector<vector<double>> Piyx_dy = getGradient(Piyx, dy, 1, theta);
    vector<vector<double>> Piyy_dy = getGradient(Piyy, dy, 1, theta);

    // extrapolate fluxes
    
    vector<vector<double>> rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR;
    tie(rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR) = extrapolateInSpaceToFace(rho, rho_dx, dx, 0);
    vector<vector<double>> vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR;   
    tie(vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR)  = extrapolateInSpaceToFace(vx,  vx_dx,   dx, 0);
    vector<vector<double>> vyM_XL,  vyP_XL,  vyM_XR,  vyP_XR;   
    tie(vyM_XL,  vyP_XL,  vyM_XR,  vyP_XR)  = extrapolateInSpaceToFace(vy,  vy_dx,   dx, 0);
    vector<vector<double>> PM_XL,   PP_XL,   PM_XR,   PP_XR;   
    tie(PM_XL,   PP_XL,   PM_XR,   PP_XR)   = extrapolateInSpaceToFace(P,   P_dx,   dx, 0);
    vector<vector<double>> PixxM_XL,  PixxP_XL,  PixxM_XR,  PixxP_XR;   
    tie(PixxM_XL,  PixxP_XL,  PixxM_XR,  PixxP_XR)  = extrapolateInSpaceToFace(Pixx,  Pixx_dx,  dx, 0);
    vector<vector<double>> PixyM_XL,  PixyP_XL,  PixyM_XR,  PixyP_XR;   
    tie(PixyM_XL,  PixyP_XL,  PixyM_XR,  PixyP_XR)  = extrapolateInSpaceToFace(Pixy,  Pixy_dx,  dx, 0);
    vector<vector<double>> PiyxM_XL,  PiyxP_XL,  PiyxM_XR,  PiyxP_XR;   
    tie(PiyxM_XL,  PiyxP_XL,  PiyxM_XR,  PiyxP_XR)  = extrapolateInSpaceToFace(Piyx,  Piyx_dx,  dx, 0);
    vector<vector<double>> PiyyM_XL,  PiyyP_XL,  PiyyM_XR,  PiyyP_XR;   
    tie(PiyyM_XL,  PiyyP_XL,  PiyyM_XR,  PiyyP_XR)  = extrapolateInSpaceToFace(Piyy,  Piyy_dx,  dx, 0);

    vector<vector<double>> rhoM_YL, rhoP_YL, rhoM_YR, rhoP_YR;
    tie(rhoM_YL, rhoP_YL, rhoM_YR, rhoP_YR) = extrapolateInSpaceToFace(rho, rho_dy, dy, 1);
    vector<vector<double>> vxM_YL,  vxP_YL,  vxM_YR,  vxP_YR;   
    tie(vxM_YL,  vxP_YL,  vxM_YR,  vxP_YR)  = extrapolateInSpaceToFace(vx,  vx_dy,  dy, 1);
    vector<vector<double>> vyM_YL,  vyP_YL,  vyM_YR,  vyP_YR;   
    tie(vyM_YL,  vyP_YL,  vyM_YR,  vyP_YR)  = extrapolateInSpaceToFace(vy,  vy_dy,  dy, 1);
    vector<vector<double>> PM_YL,   PP_YL,   PM_YR,   PP_YR;   
    tie(PM_YL,   PP_YL,   PM_YR,   PP_YR)   = extrapolateInSpaceToFace(P,   P_dy,   dy, 1);
    vector<vector<double>> PixxM_YL,  PixxP_YL,  PixxM_YR,  PixxP_YR;   
    tie(PixxM_YL,  PixxP_YL,  PixxM_YR,  PixxP_YR)  = extrapolateInSpaceToFace(Pixx,  Pixx_dy,  dy, 1);
    vector<vector<double>> PixyM_YL,  PixyP_YL,  PixyM_YR,  PixyP_YR;   
    tie(PixyM_YL,  PixyP_YL,  PixyM_YR,  PixyP_YR)  = extrapolateInSpaceToFace(Pixy,  Pixy_dy,  dy, 1);
    vector<vector<double>> PiyxM_YL,  PiyxP_YL,  PiyxM_YR,  PiyxP_YR;   
    tie(PiyxM_YL,  PiyxP_YL,  PiyxM_YR,  PiyxP_YR)  = extrapolateInSpaceToFace(Piyx,  Piyx_dy,  dy, 1);
    vector<vector<double>> PiyyM_YL,  PiyyP_YL,  PiyyM_YR,  PiyyP_YR;   
    tie(PiyyM_YL,  PiyyP_YL,  PiyyM_YR,  PiyyP_YR)  = extrapolateInSpaceToFace(Piyy,  Piyy_dy,  dy, 1);

    // Compute fluxes

    vector<vector<double>> flux_Mass_XR, flux_Momx_XR, flux_Momy_XR, flux_Pixx_vxR, flux_Pixy_vxR, flux_Piyx_vxR, flux_Piyy_vxR;
    tie(flux_Mass_XR, flux_Momx_XR, flux_Momy_XR, flux_Pixx_vxR, flux_Pixy_vxR, flux_Piyx_vxR, flux_Piyy_vxR) = getXFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR,
                                            vyP_XR, vyM_XR, PixxP_XR, PixxM_XR,
                                            PixyP_XR, PixyM_XR, PiyxP_XR, PiyxM_XR,
                                            PiyyP_XR, PiyyM_XR, PP_XR, PM_XR, gamma,
                                            eta, zeta, tau_nu);

    vector<vector<double>> flux_Mass_XL, flux_Momx_XL, flux_Momy_XL, flux_Pixx_vxL, flux_Pixy_vxL,flux_Piyx_vxL, flux_Piyy_vxL;
    tie(flux_Mass_XL, flux_Momx_XL, flux_Momy_XL, flux_Pixx_vxL, flux_Pixy_vxL,flux_Piyx_vxL, flux_Piyy_vxL) = getXFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL,
                                            vyP_XL, vyM_XL, PixxP_XL, PixxM_XL,
                                            PixyP_XL, PixyM_XL, PiyxP_XL, PiyxM_XL,
                                            PiyyP_XL, PiyyM_XL, PP_XL, PM_XL, gamma,
                                            eta, zeta, tau_nu);

    vector<vector<double>> flux_Mass_YR, flux_Momx_YR, flux_Momy_YR, flux_Pixx_vyR, flux_Pixy_vyR, flux_Piyx_vyR, flux_Piyy_vyR;
    tie(flux_Mass_YR, flux_Momx_YR, flux_Momy_YR, flux_Pixx_vyR, flux_Pixy_vyR, flux_Piyx_vyR, flux_Piyy_vyR) = getYFlux(rhoP_YR, rhoM_YR, vxP_YR, vxM_YR,
                                            vyP_YR, vyM_YR, PixxP_YR, PixxM_YR,
                                            PixyP_YR, PixyM_YR, PiyxP_YR, PiyxM_YR,
                                            PiyyP_YR, PiyyM_YR, PP_YR, PM_YR, gamma,
                                            eta, zeta, tau_nu);

    vector<vector<double>> flux_Mass_YL, flux_Momx_YL, flux_Momy_YL, flux_Pixx_vyL, flux_Pixy_vyL, flux_Piyx_vyL, flux_Piyy_vyL;
    tie(flux_Mass_YL, flux_Momx_YL, flux_Momy_YL, flux_Pixx_vyL, flux_Pixy_vyL, flux_Piyx_vyL, flux_Piyy_vyL) = getYFlux(rhoP_YL, rhoM_YL, vxP_YL, vxM_YL,
                                            vyP_YL, vyM_YL, PixxP_YL, PixxM_YL,
                                            PixyP_YL, PixyM_YL, PiyxP_YL, PiyxM_YL,
                                            PiyyP_YL, PiyyM_YL, PP_YL, PM_YL, gamma,
                                            eta, zeta, tau_nu);


    // Update conservative variables

    vector<vector<double>> J(N,vector<double>(N,0.0));
    vector<vector<double>>  Jxx(N,vector<double>(N,0.0)),Jxy(N,vector<double>(N,0.0)),Jyx(N,vector<double>(N,0.0)),Jyy(N,vector<double>(N,0.0));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
             Jxx[i][j] = -Pixx[i][j]/tau_nu;
            Jxy[i][j] = -Pixy[i][j]/tau_nu;
            Jyx[i][j] = -Piyx[i][j]/tau_nu;
            Jyy[i][j] = -Piyy[i][j]/tau_nu;
        }
    }

    vector<vector<double>> timederivative_rho, timederivative_Momx, timederivative_Momy;
    vector<vector<double>> timederivative_Pixx, timederivative_Pixy, timederivative_Piyx, timederivative_Piyy;

    timederivative_rho = applyFluxes( flux_Mass_XR,   flux_Mass_XL, flux_Mass_YR,   flux_Mass_YL,  dx, dy, J);
    timederivative_Momx = applyFluxes( flux_Momx_XR,   flux_Momx_XL, flux_Momx_YR,   flux_Momx_YL, dx, dy, J); 
    timederivative_Momy = applyFluxes( flux_Momy_XR,   flux_Momy_XL, flux_Momy_YR,   flux_Momy_YL, dx, dy, J);
    timederivative_Pixx  = applyFluxes( flux_Pixx_vxR,    flux_Pixx_vxL, flux_Pixx_vyR,    flux_Pixx_vyL, dx, dy, Jxx);
    timederivative_Pixy  = applyFluxes( flux_Pixy_vxR,    flux_Pixy_vxL, flux_Pixy_vyR,    flux_Pixy_vyL, dx, dy, Jxy);
    timederivative_Piyx  = applyFluxes( flux_Piyx_vxR,    flux_Piyx_vxL, flux_Piyx_vyR,    flux_Piyx_vyL, dx, dy, Jyx);
    timederivative_Piyy  = applyFluxes( flux_Piyy_vxR,    flux_Piyy_vxL, flux_Piyy_vyR,    flux_Piyy_vyL, dx, dy, Jyy);

    return {timederivative_rho,timederivative_Momx,timederivative_Momy,timederivative_Pixx,timederivative_Pixy,timederivative_Piyx,timederivative_Piyy};
    
}


map<double,state> integrator(state (*scheme)(double, state&, double, double, int, double, double, double, double, double), tuple<double,double> time, map<double, state> Q, double dtmax,  tuple<double, double, int, double, double, double, double, double> args, string method)
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

    double t = get<0>(time);
    double tEnd = get<1>(time);
    int outCount = 1;

    state q = Q[t];

    int N;
    double dx, dy, gamma, zeta, tau_nu, eta, theta;
    tie(dx, dy, N, gamma, zeta, tau_nu, eta, theta) = args;

    vector<vector<double>> rho(N,vector<double>(N,0.0));
    vector<vector<double>> vx(N, vector<double>(N,0.0));
    vector<vector<double>> vy(N, vector<double>(N,0.0));
    vector<vector<double>> Momx(N, vector<double>(N,0.0));
    vector<vector<double>> Momy(N, vector<double>(N,0.0));
    vector<vector<double>> cs(N, vector<double>(N,0.0));


    auto C = [&](double t, state y) {return scheme(t, q, dx, dy, N, gamma, zeta, tau_nu, eta, theta);};

    while (t < tEnd) {

        rho = q.get(0);
        Momx  = q.get(1);
        Momy  = q.get(2);
        cs  = getSpeedOfSound(rho, gamma);

        for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            vx[i][j] = Momx[i][j] / rho[i][j];
            vy[i][j] = Momy[i][j] / rho[i][j];



        }
    }        

        // condition to ensure that the time steps are small enough so that
        // waves do not interfere with each other
        vector<vector<double>> propagation_speed = local_propagation_speed(rho, vx, vy, eta, zeta, tau_nu, cs);
        vector<vector<double>> courant_number(N, vector<double>(N, 0.0));

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (propagation_speed[i][j] != 0) {
                    courant_number[i][j] = dx / propagation_speed[i][j];
                }
            }
        }


        double dt = min(dtmax, 0.4 * (max(courant_number)));

        cout << "dt: " << dt << endl;

        if (method == "Heuns") {
            q = Heuns(q, C, dt, t);
        } 

        // BC(q);/

        t = t + dt;

        if (t > outCount*dtmax) {
            Q[t] = q;
            cout << t << endl;
            ++outCount;
        }
    }

    return Q;
}


