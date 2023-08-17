#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <functional>
#include "KTmethods2d.h"


double max_value(const std::vector<std::vector<double>>& value) {

    // Returns the maximum value in a vector<vector<double>>

    if (value.empty() || value[0].empty()) {
        // Handle the case when the input matrix is empty
        throw std::runtime_error("Input matrix is empty.");
    }

    double max = value[0][0];

    for (const auto& row : value) {
        for (const double& element : row) {
            max = std::max(max, element);
        }
    }

    return max;
}

double sign(double value){   

    // Returns the mathematical sign of a double value 

    return (value > 0) - (value < 0);
}



std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>> getConserved(std::vector<std::vector<double>>& rho,
 std::vector<std::vector<double>>& vx, std::vector<std::vector<double>>& vy, double vol) {
   
    // Takes density, x-velocity, y-velocity, volume as inputs
    // Returns a tuple with the Mass, x-Momentum, y-Momentum


    std::vector<std::vector<double>> Mass(rho.size(), std::vector<double>(rho[0].size()));
    std::vector<std::vector<double>> Momx(rho.size(), std::vector<double>(rho[0].size()));
    std::vector<std::vector<double>> Momy(rho.size(), std::vector<double>(rho[0].size()));

    for (unsigned int i = 0; i < rho.size(); i++) {
        for (unsigned int j = 0; j < rho[0].size(); j++) {
            Mass[i][j] = rho[i][j] * vol;
            Momx[i][j] = rho[i][j] * vx[i][j];
            Momy[i][j] = rho[i][j] * vy[i][j];
        }
    }

    return make_tuple(Mass, Momx, Momy);
}

std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>> getPrimitive(std::vector<std::vector<double>> &Mass,
 std::vector<std::vector<double>>& Momx, std::vector<std::vector<double>>& Momy, double gamma, double vol) {
    
    // Takes Mass, x-Momentum, y-Momentum, volume as inputs
    // Returns a tuple with the density, x-velocity, y-velocity, and pressure using the equation of state


    std::vector<std::vector<double>> rho(Mass.size(), std::vector<double>(Mass[0].size()));
    std::vector<std::vector<double>> vx(Mass.size(), std::vector<double>(Mass[0].size()));
    std::vector<std::vector<double>> vy(Mass.size(), std::vector<double>(Mass[0].size()));
    std::vector<std::vector<double>> P(Mass.size(), std::vector<double>(Mass[0].size()));

    for (unsigned int i = 0; i < Mass.size(); i++) {
        for (unsigned int j = 0; j < Mass[0].size(); j++) {
            rho[i][j] = Mass[i][j] / vol;
            vx[i][j] = Momx[i][j] / rho[i][j];
            vy[i][j] = Momy[i][j] / rho[i][j];
            P[i][j] = pow(rho[i][j], gamma);

            

        }
    }

    return make_tuple(rho, vx, vy, P);
}

std::vector<std::vector<double>> getSpeedOfSound(std::vector<std::vector<double>>& rho, double gamma) {

    // get the speed of sound based on the underlying EOS

    unsigned int rows,cols;
    rows = rho.size();
    cols = rho[0].size();

    std::vector<std::vector<double>> cs(rows, std::vector<double>(cols, 0.0));

    for (unsigned int i=0; i < rows; i++){
        for (unsigned int j=0; j < cols; j++){
            cs[i][j] = sqrt(gamma * pow(rho[i][j], gamma - 1)); 
        }
    }

    return cs;
}

double minmod2(double x, double y) {
    

    double sign1 = sign(x);
    double sign2 = sign(y);
    double abs1 = std::abs(x);
    double abs2 = std::abs(y);

    double mmod = (sign1 + sign2) * std::min(abs1, abs2) / 2;

    return mmod;
       
}

double minmod3(double x, double y, double z) {
    return minmod2(x,minmod2(y,z));
}

double minmod(double x, double y, double z) {

    // minmod function to ensure TVD and therefore the stability of the scheme

    return std::max(0.0, std::min({x, y, z}));
}

std::vector<std::vector<double>> getGradient (std::vector<std::vector<double>>& f, double dx, int axis, double theta = 1) {


    // f     is a variable
    // dx    is the cell size
    // axis  is the axis over which the gradient is returned
    // theta is the flux limiter variable


    // Returns the gradient of the f variable 



    int n = f.size();
    int m = f[0].size();

    std::vector<std::vector<double>> df_dx(n, std::vector<double>(m, 0.0));
    

    if (axis == 0) {
        for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < m - 1; j++) {
                double term1 = theta * (f[i][j] - f[i-1][j]) / dx;
                double term2 = (f[i+1][j] - f[i-1][j]) / (2 * dx);
                double term3 = theta * (f[i+1][j] - f[i][j]) / dx;
                df_dx[i][j] = minmod(term1, term2, term3);
            }
        }
    }
    else if (axis == 1) {
        for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < m - 1; j++) {
                double term1 = theta * (f[i][j] - f[i][j-1]) / dx;
                double term2 = (f[i][j+1] - f[i][j-1]) / (2 * dx);
                double term3 = theta * (f[i][j+1] - f[i][j]) / dx;
                df_dx[i][j] = minmod(term1, term2, term3);
            }
        }
    }

    return df_dx;
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>,
 std::vector<std::vector<double>>, std::vector<std::vector<double>>> extrapolateInSpaceToFace (std::vector<std::vector<double>>& q,
  std::vector<std::vector<double>>& q_dx, double dx, int axis) {
    
    // q     is the variable to be extrapolated
    // q_dx  is gradient in x-axis
    // dx    is  cell size
    // axis  is the axis which the q is extrapolated 

    // Returns the tuple of linear half-steps in space 

    int n = q.size();
    int m = q[0].size();

    std::vector<std::vector<double>> qP_XL(q.size(), std::vector<double>(q[0].size(), 0.0));
    std::vector<std::vector<double>> qP_XR(q.size(), std::vector<double>(q[0].size(), 0.0));
    std::vector<std::vector<double>> qM_XR(q.size(), std::vector<double>(q[0].size(), 0.0));
    std::vector<std::vector<double>> qM_XL(q.size(), std::vector<double>(q[0].size(), 0.0));

    if (axis == 0) {
        for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < m - 1; j++) {
                qP_XL[i][j] = q[i][j] - q_dx[i][j] * dx / 2;
                qP_XR[i][j] = q[i+1][j] - q_dx[i+1][j] * dx / 2;
                qM_XR[i][j] = q[i][j] + q_dx[i][j] * dx / 2;
                qM_XL[i][j] = q[i-1][j] + q_dx[i-1][j] * dx / 2;
            }
        }
    }
    else if (axis == 1) {
        for (int i = 1; i < n-1; i++) {
            for (int j = 1; j < m-1; j++) {
                qP_XL[i][j] = q[i][j] - q_dx[i][j] * dx / 2;
                qP_XR[i][j] = q[i][j+1] - q_dx[i][j+1] * dx / 2;
                qM_XR[i][j] = q[i][j] + q_dx[i][j] * dx / 2;
                qM_XL[i][j] = q[i][j-1] + q_dx[i][j-1] * dx / 2;
            }
        }
    }

    return make_tuple(qM_XL, qP_XL, qM_XR, qP_XR);
}

std::vector<std::vector<double>> local_propagation_speed (std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& vx, std::vector<std::vector<double>>& vy, double eta, double zeta, double tau_nu, std::vector<std::vector<double>>& cs) {
   
    // Returns the maximum local propagation speed which is the spectral radius of the hyperbolic PDE

    int rows = rho.size();
    int cols = rho[0].size();

    std::vector<std::vector<double>> C1(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> C2(rows, std::vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (rho[i][j] != 0) {
                C1[i][j] = sqrt(eta * tau_nu / rho[i][j]);
                C2[i][j] = cs[i][j] * sqrt( (zeta + 4.0 / 3.0 * tau_nu) / (rho[i][j] * tau_nu));
            }
        }
    }

    std::vector<std::vector<double>> maxC(rows, std::vector<double>(cols, 0.0));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            maxC[i][j] = std::max(C1[i][j], C2[i][j]);
        }
    }

    return maxC;
}


std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>,
 std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> getXFlux(std::vector<std::vector<double>>& rho_P,
    std::vector<std::vector<double>>& rho_M, std::vector<std::vector<double>>& vx_P, std::vector<std::vector<double>>& vx_M, 
    std::vector<std::vector<double>>& vy_P, std::vector<std::vector<double>>& vy_M, std::vector<std::vector<double>>& Pixx_P, 
   std::vector<std::vector<double>>& Pixx_M, std::vector<std::vector<double>>& Pixy_P, std::vector<std::vector<double>>& Pixy_M, 
   std::vector<std::vector<double>>& Piyx_P, std::vector<std::vector<double>>& Piyx_M, std::vector<std::vector<double>>& Piyy_P, 
   std::vector<std::vector<double>>& Piyy_M, std::vector<std::vector<double>>& P_P, std::vector<std::vector<double>>& P_M, double gamma, 
  double eta, double zeta, double tau_nu) {

    int rows = rho_P.size();
    int cols = rho_P[0].size();

    std::vector<std::vector<double>> flux_Mass(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Momx(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Momy(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Pixx_vx(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Pixy_vx(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Piyx_vx(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Piyy_vx(rows, std::vector<double>(cols, 0.0));

    double vx_av;
    double vy_av;
    double momx_av;
    double Pixx_av;
    double Piyx_av;
    double Pixx_vx_av;
    double Pixy_vx_av;
    double Piyx_vx_av;
    double Piyy_vx_av;
    double P_av;

    std::vector<std::vector<double>> cs_P = getSpeedOfSound(rho_P, gamma);
    std::vector<std::vector<double>> C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P);
    std::vector<std::vector<double>> cs_M = getSpeedOfSound(rho_M, gamma);
    std::vector<std::vector<double>> C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M);

    double C;

    double B = eta / tau_nu;
    double A = zeta / tau_nu;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            vx_av = 0.5 * (vx_P[i][j] + vx_M[i][j]);
            vy_av = 0.5 * (vy_P[i][j] + vy_M[i][j]);
            momx_av = 0.5 * (rho_P[i][j] * vx_P[i][j] + rho_M[i][j] * vx_M[i][j]);
            Pixx_av = 0.5 * (Pixx_P[i][j] + Pixx_M[i][j]);
            Piyx_av = 0.5 * (Piyx_P[i][j] + Piyx_M[i][j]);
            Pixx_vx_av = 0.5 * (Pixx_P[i][j] * vx_P[i][j] + Pixx_M[i][j] * vx_M[i][j]);
            Pixy_vx_av = 0.5 * (Pixy_P[i][j] * vx_P[i][j] + Piyx_M[i][j] * vx_M[i][j]);
            Piyx_vx_av = 0.5 * (Piyx_P[i][j] * vx_P[i][j] + Pixy_M[i][j] * vx_M[i][j]);
            Piyy_vx_av = 0.5 * (Piyy_P[i][j] * vx_P[i][j] + Piyy_M[i][j] * vx_M[i][j]);
            P_av = 0.5 * (P_P[i][j] + P_M[i][j]);

            

            flux_Mass[i][j] = momx_av;
            flux_Momx[i][j] = 0.5 * (rho_P[i][j] * vx_P[i][j] * vx_P[i][j] + rho_M[i][j] * vx_M[i][j] * vx_M[i][j]) + ((P_av) + (Pixx_av)) / gamma;
            flux_Momy[i][j] = 0.5 * (rho_P[i][j] * (vx_P[i][j] * vy_P[i][j]) + rho_M[i][j] * (vx_M[i][j] * vy_M[i][j])) + (Piyx_av) / gamma;
            flux_Pixx_vx[i][j] = Pixx_vx_av + 2 * B * (vx_av) + (A - 2.0 / 3.0 * B) * vx_av;
            flux_Pixy_vx[i][j] = Pixy_vx_av + B * vy_av;
            flux_Piyx_vx[i][j] = Piyx_vx_av + B * vy_av;
            flux_Piyy_vx[i][j] = Piyy_vx_av + (A - 2.0 / 3.0 * B) * vx_av;

            C = std::max(C_M[i][j], C_P[i][j]);

            flux_Mass[i][j] -= C * 0.5 * (rho_P[i][j] - rho_M[i][j]);
            flux_Momx[i][j] -= C * 0.5 * (rho_P[i][j] * vx_P[i][j] - rho_M[i][j] * vx_M[i][j]);
            flux_Momy[i][j] -= C * 0.5 * (rho_P[i][j] * vy_P[i][j] - rho_M[i][j] * vy_M[i][j]);
            flux_Pixx_vx[i][j] -= C * 0.5 * (Pixx_P[i][j] - Pixx_M[i][j]);
            flux_Pixy_vx[i][j] -= C * 0.5 * (Pixy_P[i][j] - Pixy_M[i][j]);
            flux_Piyx_vx[i][j] -= C * 0.5 * (Piyx_P[i][j] - Piyx_M[i][j]);
            flux_Piyy_vx[i][j] -= C * 0.5 * (Piyy_P[i][j] - Piyy_M[i][j]);

        }
    }

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vx, flux_Pixy_vx, flux_Piyx_vx, flux_Piyy_vx);
  }



std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>,
 std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> getYFlux(std::vector<std::vector<double>>& rho_P, std::vector<std::vector<double>>& rho_M,
             std::vector<std::vector<double>>& vx_P, std::vector<std::vector<double>>& vx_M,
             std::vector<std::vector<double>>& vy_P, std::vector<std::vector<double>>& vy_M,
             std::vector<std::vector<double>>& Pixx_P, std::vector<std::vector<double>>& Pixx_M,
             std::vector<std::vector<double>>& Pixy_P, std::vector<std::vector<double>>& Pixy_M,
             std::vector<std::vector<double>>& Piyx_P, std::vector<std::vector<double>>& Piyx_M,
             std::vector<std::vector<double>>& Piyy_P, std::vector<std::vector<double>>& Piyy_M,
             std::vector<std::vector<double>>& P_P, std::vector<std::vector<double>>& P_M,
             double gamma, double eta,
             double zeta, double tau_nu){

    int rows = rho_P.size();
    int cols = rho_P[0].size();

    std::vector<std::vector<double>> flux_Mass(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Momx(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Momy(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Pixx_vy(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Pixy_vy(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Piyx_vy(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Piyy_vy(rows, std::vector<double>(cols, 0.0));

    double vx_av;
    double vy_av;
    double momy_av;
    double Piyy_av;
    double Pixy_av;
    double Pixx_vy_av;
    double Pixy_vy_av;
    double Piyx_vy_av;
    double Piyy_vy_av;
    double P_av;

    std::vector<std::vector<double>> cs_P = getSpeedOfSound(rho_P, gamma);
    std::vector<std::vector<double>> C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P);

    std::vector<std::vector<double>> cs_M = getSpeedOfSound(rho_M, gamma);
    std::vector<std::vector<double>> C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M);

    double C;

    double B = eta / tau_nu;
    double A = zeta/ tau_nu;


    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            vx_av = 0.5 * (vx_P[i][j] + vx_M[i][j]);
            vy_av = 0.5 * (vy_P[i][j] + vy_M[i][j]);
            momy_av = 0.5 * (rho_P[i][j] * vy_P[i][j] + rho_M[i][j] * vy_M[i][j]);
            Piyy_av = 0.5 * (Piyy_P[i][j] + Piyy_M[i][j]);
            Pixy_av = 0.5 * (Pixy_P[i][j] + Pixy_M[i][j]);
            Pixx_vy_av = 0.5 * (Pixx_P[i][j] * vy_P[i][j] + Pixx_M[i][j] * vy_M[i][j]);
            Pixy_vy_av = 0.5 * (Pixy_P[i][j] * vy_P[i][j] + Piyx_M[i][j] * vy_M[i][j]);
            Piyx_vy_av = 0.5 * (Piyx_P[i][j] * vy_P[i][j] + Pixy_M[i][j] * vy_M[i][j]);
            Piyy_vy_av = 0.5 * (Piyy_P[i][j] * vy_P[i][j] + Piyy_M[i][j] * vy_M[i][j]);
            P_av = 0.5 * (P_P[i][j] + P_M[i][j]);


            flux_Mass[i][j] = momy_av;
            flux_Momx[i][j] = 0.5 * (rho_P[i][j] * vx_P[i][j] * vy_P[i][j] + rho_M[i][j] * vx_M[i][j] * vy_M[i][j]) + Pixy_av / gamma;
            flux_Momy[i][j] = 0.5 * (rho_P[i][j] * vy_P[i][j] * vy_P[i][j] + rho_M[i][j] * vy_M[i][j] * vy_M[i][j]) + ((P_av)  + (Piyy_av)) / gamma;
            flux_Pixx_vy[i][j] = Pixx_vy_av + (A - 2.0 / 3.0 * B) * (vy_av);
            flux_Pixy_vy[i][j] = Pixy_vy_av + B * (vx_av);
            flux_Piyx_vy[i][j] = Piyx_vy_av + B * (vx_av);
            flux_Piyy_vy[i][j] = Piyy_vy_av + 2 * B * (vy_av) + (A - 2.0 / 3.0 * B) * (vy_av);

            C = std::max(C_M[i][j], C_P[i][j]);

            flux_Mass[i][j] -= C * 0.5 * (rho_P[i][j] - rho_M[i][j]);
            flux_Momx[i][j] -= C * 0.5 * (rho_P[i][j] * vx_P[i][j] - rho_M[i][j] * vx_M[i][j]);
            flux_Momy[i][j] -= C * 0.5 * (rho_P[i][j] * vy_P[i][j] - rho_M[i][j] * vy_M[i][j]);
            flux_Pixx_vy[i][j] -= C * 0.5 * (Pixx_P[i][j] - Pixx_M[i][j]);
            flux_Pixy_vy[i][j] -= C * 0.5 * (Pixy_P[i][j] - Pixy_M[i][j]);
            flux_Piyx_vy[i][j] -= C * 0.5 * (Piyx_P[i][j] - Piyx_M[i][j]);
            flux_Piyy_vy[i][j] -= C * 0.5 * (Piyy_P[i][j] - Piyy_M[i][j]);

        }
    }

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vy, flux_Pixy_vy, flux_Piyx_vy, flux_Piyy_vy);
        
}
// Apply fluxes to conserved variables

std::vector<std::vector<double>> applyFluxes(std::vector<std::vector<double>>& flux_H1_X, std::vector<std::vector<double>>& flux_H2_X,
  std::vector<std::vector<double>>& flux_H1_Y, std::vector<std::vector<double>>& flux_H2_Y,
   double dx, double dy, std::vector<std::vector<double>>& J){

    int rows = flux_H1_X.size();
    int cols = flux_H1_X[0].size();


    std::vector<std::vector<double>> C(rows, std::vector<double>(cols, 0.0));
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            C[i][j] -= (flux_H1_X[i][j] - flux_H2_X[i][j]) / dx;
            C[i][j] -= (flux_H1_Y[i][j] - flux_H2_Y[i][j]) / dy;
            C[i][j] += J[i][j];
        }
    }

    return C;
}

// Heun's method
state heuns (state& q, std::function<state(double,state)> f, double dt, double t) {
    
    int rows = (q.get(0)).size();
    int cols = (q.get(0))[0].size();

    double k1,k2;
    std::vector<std::vector<double>> c1(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> c2(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> y(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> yprime(rows, std::vector<double>(cols, 0.0));


    state qprime;
    state C1,C2;

    C1 = f(t,q);

    for (int n = 0; n < 7; ++n){
        c1 = C1.get(n);
        y = q.get(n);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                k1 = dt * c1[i][j];
                yprime[i][j] = y[i][j] + k1;
            }
        }
        qprime.set(n , yprime);
    }
    
    C2 = f(t+dt,qprime);

    for (int n = 0; n < 7; ++n){
        c1 = C1.get(n);
        c2 = C2.get(n);
        y = q.get(n);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                k1 = dt * c1[i][j];
                k2 = dt * c2[i][j];
                yprime[i][j] = y[i][j] + 0.5*(k1 + k2);
            }
        }
        qprime.set(n , yprime);
    }
    

    return qprime;
}


state rK4 (state& q, std::function<state(double,state)> f, double dt, double t) {
    
    int rows = (q.get(0)).size();
    int cols = (q.get(0))[0].size();

    double k1,k2,k3,k4;
    std::vector<std::vector<double>> c1(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> c2(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> c3(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> c4(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> y(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> yprime(rows, std::vector<double>(cols, 0.0));


    state qprime;
    state C1,C2,C3,C4;

    C1 = f(t,q);

    for (int n = 0; n < 7; ++n){
        c1 = C1.get(n);
        y = q.get(n);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                k1 = dt * c1[i][j];
                yprime[i][j] = y[i][j] + 0.5 * k1;
            }
        }
        qprime.set(n , yprime);
    }
    
    C2 = f(t+dt/2,qprime);

    for (int n = 0; n < 7; ++n){
        c2 = C2.get(n);
        y = q.get(n);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                k2 = dt * c2[i][j];
                yprime[i][j] = y[i][j] + 0.5 * k2;
            }
        }
        qprime.set(n , yprime);
    }

    C3 = f(t+dt/2,qprime);

    for (int n = 0; n < 7; ++n){
        c3 = C3.get(n);
        y = q.get(n);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                k3 = dt * c3[i][j];
                yprime[i][j] = y[i][j] + k3;
            }
        }
        qprime.set(n , yprime);
    }

    C4 = f(t+dt,qprime);

    for (int n = 0; n < 7; ++n){
        c1 = C1.get(n);
        c2 = C2.get(n);
        c3 = C3.get(n);
        c4 = C4.get(n);
        y = q.get(n);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                k1 = dt * c1[i][j];
                k2 = dt * c2[i][j];
                k3 = dt * c3[i][j];
                k4 = dt * c4[i][j];
                yprime[i][j] = y[i][j] + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
            }
        }
        qprime.set(n , yprime);
    }

    return qprime;
}
