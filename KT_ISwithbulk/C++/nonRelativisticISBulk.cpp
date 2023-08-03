#include <cmath>
#include <vector>
#include <tuple>
#include <list>
#include <iostream>
#include <string>
#include "KTmethods2dBulk.h"
#include "nonRelativisticISBulk.h"

stateb KTschemeNonRelativisticIS(double t,  stateb& IC, double dx, double dy, int N, double gamma, double zeta, double tau_nu, double theta = 1) {
     
    /* Finite Volume simulation */

    // Generate Initial Conditions

    /* Initial conditions for rho */ 
    std::vector<std::vector<double>> rho(N,std::vector<double>(N,0.0));
    std::vector<std::vector<double>> P(N,std::vector<double>(N,0.0));

    /* Initial conditions for v */
    std::vector<std::vector<double>> Momx(N, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> Momy(N, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> vx(N, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> vy(N, std::vector<double>(N,0.0));

    /* Pi initial condition */
    std::vector<std::vector<double>> Pi(N,std::vector<double>(N,0.0));
  

    rho  = IC.get(0);
    Momx = IC.get(1);
    Momy = IC.get(2);
    Pi   = IC.get(3);
    

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            vx[i][j] = Momx[i][j] / rho[i][j];
            vy[i][j] = Momy[i][j] / rho[i][j];

             /* Pressure from equation of stateb */
            P[i][j] = pow(rho[i][j], gamma);

            // set boundary conditions
            rho[0][j] = rho[1][j];
            rho[i][0] = rho[i][1];
            rho[N-1][j] = rho[N-2][j];
            rho[i][N-1] = rho[i][N-2];

            vx[0][j] = -vx[1][j];
            vx[i][0] = 0;
            vx[N-1][j] = -vx[N-2][j];
            vx[i][N-1] = 0;

            vy[0][j] = 0;
            vy[i][0] = -vy[i][1];
            vy[N-1][j] = 0;
            vy[i][N-1] = -vy[i][N-2];

            P[0][j] = P[1][j];
            P[i][0] = P[i][1];
            P[N-1][j] = P[N-2][j];
            P[i][N-1] = P[i][N-2];

            Pi[0][j] = Pi[1][j];
            Pi[i][0] = Pi[i][1];
            Pi[N-1][j] = Pi[N-2][j];
            Pi[i][N-1] = Pi[i][N-2];

           
        }
    }        



    // calculate speed of sound
    // getSpeedOfSound(rho, gamma)
    std::vector<std::vector<double>> cs = getSpeedOfSound(rho, gamma);

    // calculate gradients
    std::vector<std::vector<double>> rho_dx = getGradient(rho, dx, 0, theta);
    std::vector<std::vector<double>> vx_dx = getGradient(vx, dx, 0, theta);
    std::vector<std::vector<double>> vy_dx = getGradient(vy, dx, 0, theta);
    std::vector<std::vector<double>> P_dx = getGradient(P, dx, 0, theta);
    std::vector<std::vector<double>> Pi_dx = getGradient(Pi, dx, 0, theta);

    std::vector<std::vector<double>> rho_dy = getGradient(rho, dy, 1, theta);
    std::vector<std::vector<double>> vx_dy = getGradient(vx, dy, 1, theta);
    std::vector<std::vector<double>> vy_dy = getGradient(vy, dy, 1, theta);
    std::vector<std::vector<double>> P_dy = getGradient(P, dy, 1, theta);
    std::vector<std::vector<double>> Pi_dy = getGradient(Pi, dy, 1, theta);
    

    // extrapolate fluxes

    std::vector<std::vector<double>> rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR;
    std::tie(rhoM_XL, rhoP_XL, rhoM_XR, rhoP_XR) = extrapolateInSpaceToFace(rho, rho_dx, dx, 0);
    std::vector<std::vector<double>> vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR;   
    std::tie(vxM_XL,  vxP_XL,  vxM_XR,  vxP_XR)  = extrapolateInSpaceToFace(vx,  vx_dx,   dx, 0);
    std::vector<std::vector<double>> vyM_XL,  vyP_XL,  vyM_XR,  vyP_XR;   
    std::tie(vyM_XL,  vyP_XL,  vyM_XR,  vyP_XR)  = extrapolateInSpaceToFace(vy,  vy_dx,   dx, 0);
    std::vector<std::vector<double>> PM_XL,   PP_XL,   PM_XR,   PP_XR;   
    std::tie(PM_XL,   PP_XL,   PM_XR,   PP_XR)   = extrapolateInSpaceToFace(P,   P_dx,   dx, 0);
    std::vector<std::vector<double>> PiM_XL,  PiP_XL,  PiM_XR,  PiP_XR;   
    std::tie(PiM_XL,  PiP_XL,  PiM_XR,  PiP_XR)  = extrapolateInSpaceToFace(Pi,  Pi_dx,  dx, 0);

    std::vector<std::vector<double>> rhoM_YL, rhoP_YL, rhoM_YR, rhoP_YR;
    std::tie(rhoM_YL, rhoP_YL, rhoM_YR, rhoP_YR) = extrapolateInSpaceToFace(rho, rho_dy, dy, 1);
    std::vector<std::vector<double>> vxM_YL,  vxP_YL,  vxM_YR,  vxP_YR;   
    std::tie(vxM_YL,  vxP_YL,  vxM_YR,  vxP_YR)  = extrapolateInSpaceToFace(vx,  vx_dy,  dy, 1);
    std::vector<std::vector<double>> vyM_YL,  vyP_YL,  vyM_YR,  vyP_YR;   
    std::tie(vyM_YL,  vyP_YL,  vyM_YR,  vyP_YR)  = extrapolateInSpaceToFace(vy,  vy_dy,  dy, 1);
    std::vector<std::vector<double>> PM_YL,   PP_YL,   PM_YR,   PP_YR;   
    std::tie(PM_YL,   PP_YL,   PM_YR,   PP_YR)   = extrapolateInSpaceToFace(P,   P_dy,   dy, 1);
    std::vector<std::vector<double>> PiM_YL,  PiP_YL,  PiM_YR,  PiP_YR;   
    std::tie(PiM_YL,  PiP_YL,  PiM_YR,  PiP_YR)  = extrapolateInSpaceToFace(Pi,  Pi_dy,  dy, 1);

    // Compute fluxes

    std::vector<std::vector<double>> flux_Mass_XR, flux_Momx_XR, flux_Momy_XR, flux_Pi_vxR;
    std::tie(flux_Mass_XR, flux_Momx_XR, flux_Momy_XR, flux_Pi_vxR) = getXFlux(rhoP_XR, rhoM_XR, vxP_XR, vxM_XR,
                                            vyP_XR, vyM_XR, PiP_XR, PiM_XR, PP_XR, PM_XR, gamma, zeta, tau_nu);

    std::vector<std::vector<double>> flux_Mass_XL, flux_Momx_XL, flux_Momy_XL, flux_Pi_vxL;
    std::tie(flux_Mass_XL, flux_Momx_XL, flux_Momy_XL, flux_Pi_vxL) = getXFlux(rhoP_XL, rhoM_XL, vxP_XL, vxM_XL,
                                            vyP_XL, vyM_XL, PiP_XL, PiM_XL, PP_XL, PM_XL, gamma, zeta, tau_nu);

    std::vector<std::vector<double>> flux_Mass_YR, flux_Momx_YR, flux_Momy_YR, flux_Pi_vyR;
    std::tie(flux_Mass_YR, flux_Momx_YR, flux_Momy_YR, flux_Pi_vyR) = getYFlux(rhoP_YR, rhoM_YR, vxP_YR, vxM_YR,
                                            vyP_YR, vyM_YR, PiP_YR, PiM_YR, PP_YR, PM_YR, gamma, zeta, tau_nu);

    std::vector<std::vector<double>> flux_Mass_YL, flux_Momx_YL, flux_Momy_YL, flux_Pi_vyL;
    std::tie(flux_Mass_YL, flux_Momx_YL, flux_Momy_YL, flux_Pi_vyL) = getYFlux(rhoP_YL, rhoM_YL, vxP_YL, vxM_YL,
                                            vyP_YL, vyM_YL, PiP_YL, PiM_YL, PP_YL, PM_YL, gamma, zeta, tau_nu);
    
    
    
    // Update conservative variables

    std::vector<std::vector<double>> J(N,std::vector<double>(N,0.0));
    std::vector<std::vector<double>>  Jxx(N,std::vector<double>(N,0.0));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Jxx[i][j] = -Pi[i][j]/tau_nu;

        }
    }

    std::vector<std::vector<double>> timederivative_rho, timederivative_Momx, timederivative_Momy;
    std::vector<std::vector<double>> timederivative_Pi;

    timederivative_rho  = applyFluxes( flux_Mass_XR,   flux_Mass_XL, flux_Mass_YR,   flux_Mass_YL, dx, dy, J);
    timederivative_Momx = applyFluxes( flux_Momx_XR,   flux_Momx_XL, flux_Momx_YR,   flux_Momx_YL, dx, dy, J); 
    timederivative_Momy = applyFluxes( flux_Momy_XR,   flux_Momy_XL, flux_Momy_YR,   flux_Momy_YL, dx, dy, J);
    timederivative_Pi   = applyFluxes( flux_Pi_vxR,    flux_Pi_vxL,  flux_Pi_vyR,    flux_Pi_vyL,  dx, dy, Jxx);



    

   
    
    return {timederivative_rho,timederivative_Momx,timederivative_Momy,timederivative_Pi};
    
}


std::list<stateb> integrator(stateb (*scheme)(double, stateb&, double, double, int, double, double, double, double), std::tuple<double,double> time, std::list<stateb> Q, double dtmax,  std::tuple<double, double, int, double, double, double, double> args, std::string method)
{
    /*
    This is an integrator that evolves a

    scheme     is the method to get dy/dt e.g. KTscheme
    time       is the current time
    q0         is the current stateb
    dtmax      is the upperbound of dt set by the user
    BC         is a function that enforces the boundary conditions
    method     is the method used in the integrator
    args       are additional arguments for scheme
    */

    std::cout.precision(2);

    double t = std::get<0>(time);
    double tEnd = std::get<1>(time);
    int outCount = 1;

    stateb q = Q.back();
    stateb qprime;

    int N;
    double dx, dy, gamma, zeta, tau_nu, theta;
    std::tie(dx, dy, N, gamma, zeta, tau_nu, theta) = args;

    std::vector<std::vector<double>> rho(N,std::vector<double>(N,0.0));
    std::vector<std::vector<double>> vx(N, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> vy(N, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> Momx(N, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> Momy(N, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> cs(N, std::vector<double>(N,0.0));

    auto C = [&](double t, stateb y) {return scheme(t, y, dx, dy, N, gamma, zeta, tau_nu, theta);};

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
        std::vector<std::vector<double>> propagation_speed = local_propagation_speed(rho, vx, vy, zeta, tau_nu, cs);
        double courant_number;

        
        courant_number = dx / max_value(propagation_speed);
            


        double dt = std::min(dtmax, 0.4 * (courant_number));

        

        if (method == "Heuns") {
            qprime = heuns(q, C, dt, t);
        } if (method == "RK4") {
            qprime = rK4(q, C, dt, t);
        } 

        
        
        //cout << dt << endl;
        t = t + dt;
        q = qprime;

        if (t > outCount*dtmax) {
            Q.push_back(qprime);
            std::cout << t << '/' << tEnd << std::endl;
            ++outCount;
        }
    }

    return Q;
}


