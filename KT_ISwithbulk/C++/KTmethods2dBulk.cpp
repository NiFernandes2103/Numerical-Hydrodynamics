#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <functional>
#include "KTmethods2dBulk.h"


double max_value(const std::vector<std::vector<double>>& value) {
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
    return (value > 0) - (value < 0);
}



std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>> getConserved(std::vector<std::vector<double>>& rho,
 std::vector<std::vector<double>>& vx, std::vector<std::vector<double>>& vy, double vol) {
   
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

double minmod(double x, double y, double z) {
    return std::max(0.0, std::min({x, y, z}));
}

std::vector<std::vector<double>> getGradient (std::vector<std::vector<double>>& f, double dx, int axis, double theta = 1) {
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

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> extrapolateInSpaceToFace (std::vector<std::vector<double>>& q, std::vector<std::vector<double>>& q_dx, double dx, int axis) {
    
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

std::vector<std::vector<double>> local_propagation_speed (std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& vx, std::vector<std::vector<double>>& vy, double zeta, double tau_nu, std::vector<std::vector<double>>& cs) {
   
    int rows = rho.size();
    int cols = rho[0].size();

    std::vector<std::vector<double>> C1(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> C2(rows, std::vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (rho[i][j] != 0) {
                C1[i][j] = std::abs(vx[i][j]) +  sqrt(cs[i][j]*cs[i][j] + zeta/tau_nu);
                C2[i][j] = std::abs(vy[i][j]) +  sqrt(cs[i][j]*cs[i][j] + zeta/tau_nu);
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


std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>,
 std::vector<std::vector<double>>> getXFlux(std::vector<std::vector<double>>& rho_P,
    std::vector<std::vector<double>>& rho_M, std::vector<std::vector<double>>& vx_P, std::vector<std::vector<double>>& vx_M, 
    std::vector<std::vector<double>>& vy_P, std::vector<std::vector<double>>& vy_M, std::vector<std::vector<double>>& Pi_P, 
    std::vector<std::vector<double>>& Pi_M, std::vector<std::vector<double>>& P_P, std::vector<std::vector<double>>& P_M, double gamma,
    double zeta, double tau_nu) {

    int rows = rho_P.size();
    int cols = rho_P[0].size();

    std::vector<std::vector<double>> flux_Mass(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Momx(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Momy(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Pi_vx(rows, std::vector<double>(cols, 0.0));

    double vx_av;
    double vy_av;
    double momx_av;
    double Pi_av;
    double Pi_vx_av;
    double P_av;

    std::vector<std::vector<double>> cs_P = getSpeedOfSound(rho_P, gamma);
    std::vector<std::vector<double>> C_P = local_propagation_speed(rho_P, vx_P, vy_P, zeta, tau_nu, cs_P);
    std::vector<std::vector<double>> cs_M = getSpeedOfSound(rho_M, gamma);
    std::vector<std::vector<double>> C_M = local_propagation_speed(rho_M, vx_M, vy_M, zeta, tau_nu, cs_M);

    double C;

    double B = zeta / tau_nu;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            vx_av = 0.5 * (vx_P[i][j] + vx_M[i][j]);
            vy_av = 0.5 * (vy_P[i][j] + vy_M[i][j]);
            momx_av = 0.5 * (rho_P[i][j] * vx_P[i][j] + rho_M[i][j] * vx_M[i][j]);
            Pi_av = 0.5 * (Pi_P[i][j] + Pi_M[i][j]);
            Pi_vx_av = 0.5 * (Pi_P[i][j] * vx_P[i][j] + Pi_M[i][j] * vx_M[i][j]);
            P_av = 0.5 * (P_P[i][j] + P_M[i][j]);

            

            flux_Mass[i][j] = momx_av;
            flux_Momx[i][j] = 0.5 * (rho_P[i][j] * vx_P[i][j] * vx_P[i][j] + rho_M[i][j] * vx_M[i][j] * vx_M[i][j]) + ((P_av) + (Pi_av)) / gamma;
            flux_Momy[i][j] = 0.5 * (rho_P[i][j] * (vx_P[i][j] * vy_P[i][j]) + rho_M[i][j] * (vx_M[i][j] * vy_M[i][j]));
            flux_Pi_vx[i][j] = Pi_vx_av + B * (vx_av);

            C = std::max(C_M[i][j], C_P[i][j]);

            flux_Mass[i][j] -= C * 0.5 * (rho_P[i][j] - rho_M[i][j]);
            flux_Momx[i][j] -= C * 0.5 * (rho_P[i][j] * vx_P[i][j] - rho_M[i][j] * vx_M[i][j]);
            flux_Momy[i][j] -= C * 0.5 * (rho_P[i][j] * vy_P[i][j] - rho_M[i][j] * vy_M[i][j]);
            flux_Pi_vx[i][j] -= C * 0.5 * (Pi_P[i][j] - Pi_M[i][j]);

        }
    }

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pi_vx);
  }



std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> getYFlux(std::vector<std::vector<double>>& rho_P, std::vector<std::vector<double>>& rho_M,
             std::vector<std::vector<double>>& vx_P, std::vector<std::vector<double>>& vx_M,
             std::vector<std::vector<double>>& vy_P, std::vector<std::vector<double>>& vy_M,
             std::vector<std::vector<double>>& Pi_P, std::vector<std::vector<double>>& Pi_M,
             std::vector<std::vector<double>>& P_P, std::vector<std::vector<double>>& P_M,
             double gamma, double zeta, double tau_nu){

    int rows = rho_P.size();
    int cols = rho_P[0].size();

    std::vector<std::vector<double>> flux_Mass(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Momx(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Momy(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> flux_Pi_vy(rows, std::vector<double>(cols, 0.0));
   

    double vx_av;
    double vy_av;
    double momy_av;
    double Pi_av;
    double Pi_vy_av;
    double P_av;

    std::vector<std::vector<double>> cs_P = getSpeedOfSound(rho_P, gamma);
    std::vector<std::vector<double>> C_P = local_propagation_speed(rho_P, vx_P, vy_P, zeta, tau_nu, cs_P);

    std::vector<std::vector<double>> cs_M = getSpeedOfSound(rho_M, gamma);
    std::vector<std::vector<double>> C_M = local_propagation_speed(rho_M, vx_M, vy_M, zeta, tau_nu, cs_M);

    double C;

    double B = zeta / tau_nu;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            vx_av = 0.5 * (vx_P[i][j] + vx_M[i][j]);
            vy_av = 0.5 * (vy_P[i][j] + vy_M[i][j]);
            momy_av = 0.5 * (rho_P[i][j] * vy_P[i][j] + rho_M[i][j] * vy_M[i][j]);
            Pi_av = 0.5 * (Pi_P[i][j] + Pi_M[i][j]);
            Pi_vy_av = 0.5 * (Pi_P[i][j] * vy_P[i][j] + Pi_M[i][j] * vy_M[i][j]);
            P_av = 0.5 * (P_P[i][j] + P_M[i][j]);


            flux_Mass[i][j] = momy_av;
            flux_Momx[i][j] = 0.5 * (rho_P[i][j] * vx_P[i][j] * vy_P[i][j] + rho_M[i][j] * vx_M[i][j] * vy_M[i][j]);
            flux_Momy[i][j] = 0.5 * (rho_P[i][j] * vy_P[i][j] * vy_P[i][j] + rho_M[i][j] * vy_M[i][j] * vy_M[i][j]) + ((P_av)  + (Pi_av)) / gamma;
            flux_Pi_vy[i][j] = Pi_vy_av + (B) * (vy_av);

            C = std::max(C_M[i][j], C_P[i][j]);

            flux_Mass[i][j] -= C * 0.5 * (rho_P[i][j] - rho_M[i][j]);
            flux_Momx[i][j] -= C * 0.5 * (rho_P[i][j] * vx_P[i][j] - rho_M[i][j] * vx_M[i][j]);
            flux_Momy[i][j] -= C * 0.5 * (rho_P[i][j] * vy_P[i][j] - rho_M[i][j] * vy_M[i][j]);
            flux_Pi_vy[i][j] -= C * 0.5 * (Pi_P[i][j] - Pi_M[i][j]);

        }
    }

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pi_vy);
        
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
stateb heuns (stateb& q, std::function<stateb(double,stateb)> f, double dt, double t) {
    
    int rows = (q.get(0)).size();
    int cols = (q.get(0))[0].size();

    double k1,k2;
    std::vector<std::vector<double>> c1(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> c2(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> y(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> yprime(rows, std::vector<double>(cols, 0.0));


    stateb qprime;
    stateb C1,C2;

    C1 = f(t,q);

    for (int n = 0; n < 4; ++n){
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

    for (int n = 0; n < 4; ++n){
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


stateb rK4 (stateb& q, std::function<stateb(double,stateb)> f, double dt, double t) {
    
    int rows = (q.get(0)).size();
    int cols = (q.get(0))[0].size();

    double k1,k2,k3,k4;
    std::vector<std::vector<double>> c1(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> c2(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> c3(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> c4(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> y(rows, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> yprime(rows, std::vector<double>(cols, 0.0));


    stateb qprime;
    stateb C1,C2,C3,C4;

    C1 = f(t,q);

    for (int n = 0; n < 4; ++n){
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

    for (int n = 0; n < 4; ++n){
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

    for (int n = 0; n < 4; ++n){
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

    for (int n = 0; n < 4; ++n){
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
