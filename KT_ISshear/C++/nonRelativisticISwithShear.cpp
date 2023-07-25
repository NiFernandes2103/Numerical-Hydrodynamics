#include <cmath>
#include <tuple>
#include <list>
#include <functional>
#include <iostream>
#include "matrix.h"
#include "KTmethods2d.h"
#include "nonRelativisticISwithShear.h"
using namespace std;

state KTschemeNonRelativisticIS(double t,  state IC, double dx, double dy, int N, double gamma, double zeta, double tau_nu, double eta, double theta = 1) {
     
    /* Finite Volume simulation */

    // Generate Initial Conditions

    /* Initial conditions for rho */ 
    smatrix rho(N);
    smatrix P(N);

    /* Initial conditions for v */
    smatrix Momx(N);
    smatrix Momy(N);
    smatrix vx(N);
    smatrix vy(N);

    /* Pi initial condition */
    smatrix Pixx(N);
    smatrix Pixy(N);
    smatrix Piyx(N);
    smatrix Piyy(N);

    rho  = IC.get(0);
    Momx = IC.get(1);
    Momy = IC.get(2);
    Pixx = IC.get(3);
    Pixy = IC.get(4);
    Piyx = IC.get(5);
    Piyy = IC.get(6);
    
   
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            vx.set(Momx.get(i,j) / rho.get(i,j), i,j);
            vy.set(Momy.get(i,j) / rho.get(i,j), i,j);

             /* Pressure from equation of state */
            P.set(pow(rho.get(i,j), gamma),i,j);

        }
    }        



    // calculate gradients
    smatrix rho_dx = getGradient(rho, dx, 0, theta);
    smatrix vx_dx = getGradient(vx, dx, 0, theta);
    smatrix vy_dx = getGradient(vy, dx, 0, theta);
    smatrix P_dx = getGradient(P, dx, 0, theta);
    smatrix Pixx_dx = getGradient(Pixx, dx, 0, theta);
    smatrix Pixy_dx = getGradient(Pixy, dx, 0, theta);
    smatrix Piyx_dx = getGradient(Piyx, dx, 0, theta);
    smatrix Piyy_dx = getGradient(Piyy, dx, 0, theta);

    smatrix rho_dy = getGradient(rho, dy, 1, theta);
    smatrix vx_dy = getGradient(vx, dy, 1, theta);
    smatrix vy_dy = getGradient(vy, dy, 1, theta);
    smatrix P_dy = getGradient(P, dy, 1, theta);
    smatrix Pixx_dy = getGradient(Pixx, dy, 1, theta);
    smatrix Pixy_dy = getGradient(Pixy, dy, 1, theta);
    smatrix Piyx_dy = getGradient(Piyx, dy, 1, theta);
    smatrix Piyy_dy = getGradient(Piyy, dy, 1, theta);

    
    
    
    // extrapolate fluxes

    smatrix rhoM_XL(N), rhoP_XL(N), rhoM_XR(N), rhoP_XR(N);
    tie(rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR) = extrapolateInSpaceToFace(rho, rho_dx, dx, 0);
    smatrix vxM_XL(N),  vxP_XL(N),  vxM_XR(N),  vxP_XR(N);   
    tie(vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR)  = extrapolateInSpaceToFace(vx,  vx_dx,   dx, 0);
    smatrix vyM_XL(N),  vyP_XL(N),  vyM_XR(N),  vyP_XR(N);   
    tie(vyM_XL,  vyP_XL,  vyM_XR,  vyP_XR)  = extrapolateInSpaceToFace(vy,  vy_dx,   dx, 0);
    smatrix PM_XL(N),   PP_XL(N),   PM_XR(N),   PP_XR(N);   
    tie(PM_XL,   PP_XL,   PM_XR,   PP_XR)   = extrapolateInSpaceToFace(P,   P_dx,   dx, 0);
    smatrix PixxM_XL(N),  PixxP_XL(N),  PixxM_XR(N),  PixxP_XR(N);   
    tie(PixxM_XL,  PixxP_XL,  PixxM_XR,  PixxP_XR)  = extrapolateInSpaceToFace(Pixx,  Pixx_dx,  dx, 0);
    smatrix PixyM_XL(N),  PixyP_XL(N),  PixyM_XR(N),  PixyP_XR(N);   
    tie(PixyM_XL,  PixyP_XL,  PixyM_XR,  PixyP_XR)  = extrapolateInSpaceToFace(Pixy,  Pixy_dx,  dx, 0);
    smatrix PiyxM_XL(N),  PiyxP_XL(N),  PiyxM_XR(N),  PiyxP_XR(N);   
    tie(PiyxM_XL,  PiyxP_XL,  PiyxM_XR,  PiyxP_XR)  = extrapolateInSpaceToFace(Piyx,  Piyx_dx,  dx, 0);
    smatrix PiyyM_XL(N),  PiyyP_XL(N),  PiyyM_XR(N),  PiyyP_XR(N);   
    tie(PiyyM_XL,  PiyyP_XL,  PiyyM_XR,  PiyyP_XR)  = extrapolateInSpaceToFace(Piyy,  Piyy_dx,  dx, 0);

    smatrix rhoM_YL(N), rhoP_YL(N), rhoM_YR(N), rhoP_YR(N);
    tie(rhoM_YL, rhoP_YL, rhoM_YR, rhoP_YR) = extrapolateInSpaceToFace(rho, rho_dy, dy, 1);
    smatrix vxM_YL(N),  vxP_YL(N),  vxM_YR(N),  vxP_YR(N);   
    tie(vxM_YL,  vxP_YL,  vxM_YR,  vxP_YR)  = extrapolateInSpaceToFace(vx,  vx_dy,  dy, 1);
    smatrix vyM_YL(N),  vyP_YL(N),  vyM_YR(N),  vyP_YR(N);   
    tie(vyM_YL,  vyP_YL,  vyM_YR,  vyP_YR)  = extrapolateInSpaceToFace(vy,  vy_dy,  dy, 1);
    smatrix PM_YL(N),   PP_YL(N),   PM_YR(N),   PP_YR(N);   
    tie(PM_YL,   PP_YL,   PM_YR,   PP_YR)   = extrapolateInSpaceToFace(P,   P_dy,   dy, 1);
    smatrix PixxM_YL(N),  PixxP_YL(N),  PixxM_YR(N),  PixxP_YR(N);   
    tie(PixxM_YL,  PixxP_YL,  PixxM_YR,  PixxP_YR)  = extrapolateInSpaceToFace(Pixx,  Pixx_dy,  dy, 1);
    smatrix PixyM_YL(N),  PixyP_YL(N),  PixyM_YR(N),  PixyP_YR(N);   
    tie(PixyM_YL,  PixyP_YL,  PixyM_YR,  PixyP_YR)  = extrapolateInSpaceToFace(Pixy,  Pixy_dy,  dy, 1);
    smatrix PiyxM_YL(N),  PiyxP_YL(N),  PiyxM_YR(N),  PiyxP_YR(N);   
    tie(PiyxM_YL,  PiyxP_YL,  PiyxM_YR,  PiyxP_YR)  = extrapolateInSpaceToFace(Piyx,  Piyx_dy,  dy, 1);
    smatrix PiyyM_YL(N),  PiyyP_YL(N),  PiyyM_YR(N),  PiyyP_YR(N);   
    tie(PiyyM_YL,  PiyyP_YL,  PiyyM_YR,  PiyyP_YR)  = extrapolateInSpaceToFace(Piyy,  Piyy_dy,  dy, 1);

    // Compute fluxes

    smatrix flux_Mass_XR(N), flux_Momx_XR(N), flux_Momy_XR(N), flux_Pixx_vxR(N), flux_Pixy_vxR(N), flux_Piyx_vxR(N), flux_Piyy_vxR(N);
    tie(flux_Mass_XR, flux_Momx_XR, flux_Momy_XR, flux_Pixx_vxR, flux_Pixy_vxR, flux_Piyx_vxR, flux_Piyy_vxR) = getXFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR,
                                            vyP_XR, vyM_XR, PixxP_XR, PixxM_XR,
                                            PixyP_XR, PixyM_XR, PiyxP_XR, PiyxM_XR,
                                            PiyyP_XR, PiyyM_XR, PP_XR, PM_XR, gamma,
                                            eta, zeta, tau_nu);

    smatrix flux_Mass_XL(N), flux_Momx_XL(N), flux_Momy_XL(N), flux_Pixx_vxL(N), flux_Pixy_vxL(N), flux_Piyx_vxL(N), flux_Piyy_vxL(N);
    tie(flux_Mass_XL, flux_Momx_XL, flux_Momy_XL, flux_Pixx_vxL, flux_Pixy_vxL,flux_Piyx_vxL, flux_Piyy_vxL) = getXFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL,
                                            vyP_XL, vyM_XL, PixxP_XL, PixxM_XL,
                                            PixyP_XL, PixyM_XL, PiyxP_XL, PiyxM_XL,
                                            PiyyP_XL, PiyyM_XL, PP_XL, PM_XL, gamma,
                                            eta, zeta, tau_nu);

    smatrix flux_Mass_YR(N), flux_Momx_YR(N), flux_Momy_YR(N), flux_Pixx_vyR(N), flux_Pixy_vyR(N), flux_Piyx_vyR(N), flux_Piyy_vyR(N);
    tie(flux_Mass_YR, flux_Momx_YR, flux_Momy_YR, flux_Pixx_vyR, flux_Pixy_vyR, flux_Piyx_vyR, flux_Piyy_vyR) = getYFlux(rhoP_YR, rhoM_YR, vxP_YR, vxM_YR,
                                            vyP_YR, vyM_YR, PixxP_YR, PixxM_YR,
                                            PixyP_YR, PixyM_YR, PiyxP_YR, PiyxM_YR,
                                            PiyyP_YR, PiyyM_YR, PP_YR, PM_YR, gamma,
                                            eta, zeta, tau_nu);

    smatrix flux_Mass_YL(N), flux_Momx_YL(N), flux_Momy_YL(N), flux_Pixx_vyL(N), flux_Pixy_vyL(N), flux_Piyx_vyL(N), flux_Piyy_vyL(N);
    tie(flux_Mass_YL, flux_Momx_YL, flux_Momy_YL, flux_Pixx_vyL, flux_Pixy_vyL, flux_Piyx_vyL, flux_Piyy_vyL) = getYFlux(rhoP_YL, rhoM_YL, vxP_YL, vxM_YL,
                                            vyP_YL, vyM_YL, PixxP_YL, PixxM_YL,
                                            PixyP_YL, PixyM_YL, PiyxP_YL, PiyxM_YL,
                                            PiyyP_YL, PiyyM_YL, PP_YL, PM_YL, gamma,
                                            eta, zeta, tau_nu);
    
    
    
    // Update conservative variables

    smatrix J(N);
    smatrix  Jxx(N),Jxy(N),Jyx(N),Jyy(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Jxx.set(-Pixx.get(i,j)/tau_nu, i,j);
            Jxy.set(-Pixy.get(i,j)/tau_nu, i,j);
            Jyx.set(-Piyx.get(i,j)/tau_nu, i,j);
            Jyy.set(-Piyy.get(i,j)/tau_nu, i,j);
        }
    }

    smatrix timederivative_rho(N), timederivative_Momx(N), timederivative_Momy(N);
    smatrix timederivative_Pixx(N), timederivative_Pixy(N), timederivative_Piyx(N), timederivative_Piyy(N);

    timederivative_rho = applyFluxes( flux_Mass_XR,   flux_Mass_XL, flux_Mass_YR,   flux_Mass_YL,  dx, dy, J);
    timederivative_Momx = applyFluxes( flux_Momx_XR,   flux_Momx_XL, flux_Momx_YR,   flux_Momx_YL, dx, dy, J); 
    timederivative_Momy = applyFluxes( flux_Momy_XR,   flux_Momy_XL, flux_Momy_YR,   flux_Momy_YL, dx, dy, J);
    timederivative_Pixx  = applyFluxes( flux_Pixx_vxR,    flux_Pixx_vxL, flux_Pixx_vyR,    flux_Pixx_vyL, dx, dy, Jxx);
    timederivative_Pixy  = applyFluxes( flux_Pixy_vxR,    flux_Pixy_vxL, flux_Pixy_vyR,    flux_Pixy_vyL, dx, dy, Jxy);
    timederivative_Piyx  = applyFluxes( flux_Piyx_vxR,    flux_Piyx_vxL, flux_Piyx_vyR,    flux_Piyx_vyL, dx, dy, Jyx);
    timederivative_Piyy  = applyFluxes( flux_Piyy_vxR,    flux_Piyy_vxL, flux_Piyy_vyR,    flux_Piyy_vyL, dx, dy, Jyy);

   
    
    return {timederivative_rho,timederivative_Momx,timederivative_Momy,timederivative_Pixx,timederivative_Pixy,timederivative_Piyx,timederivative_Piyy};
    
}


list<state> integrator(state (*scheme)(double, state, double, double, int, double, double, double, double, double), tuple<double,double> time, list<state>& Q, double dtmax,  tuple<double, double, int, double, double, double, double, double> args, string method)
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

    int N;
    double dx, dy, gamma, zeta, tau_nu, eta, theta;
    tie(dx, dy, N, gamma, zeta, tau_nu, eta, theta) = args;
    
    state qprime(N);

    smatrix rho(N);
    smatrix vx(N);
    smatrix vy(N);
    smatrix Momx(N);
    smatrix Momy(N);
    smatrix cs(N);

    
    auto C = [&](double t, state y) {return scheme(t, y, dx, dy, N, gamma, zeta, tau_nu, eta, theta);};

    while (t < tEnd) {

        rho   = q.get(0);
        Momx  = q.get(1);
        Momy  = q.get(2);
        cs    = getSpeedOfSound(rho, gamma);

        for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            vx.set( Momx.get(i,j) / rho.get(i,j), i,j);
            vy.set( Momy.get(i,j) / rho.get(i,j), i,j);
        }
    }        
        // condition to ensure that the time steps are small enough so that
        // waves do not interfere with each other
        smatrix propagation_speed = local_propagation_speed(rho, vx, vy, eta, zeta, tau_nu, cs);
        double courant_number;

        courant_number = dx / propagation_speed.max();            

        double dt = std::min(dtmax, 0.4 * (courant_number));

        if (method == "Heuns") {
            qprime = heuns(q, C, dt, t);
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


