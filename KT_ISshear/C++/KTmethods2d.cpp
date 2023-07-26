#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <functional>
#include "KTmethods2d.h"
using namespace std;



double max_value(double** value, int rows, int cols)
{  
    
    double max = value[0][0];

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            if (value[i][j] > max){
                max = value[i][j];
            }
        }
    }
    return max;

}

double sign(double value){   
    
    if (value > 0.0)
    {
        return 1.0;
    } else if (value < 0.0)
    {
        return -1.0;
    } else 
    {
        return 0.0;
    }
    
    
    
}



tuple<double**,double**,double**> getConserved(double**& rho,
 double**& vx, double**& vy, double vol, int rows, int cols) {
   
    double** Mass = new double*[rows];
    double** Momx = new double*[rows];
    double** Momy = new double*[rows];

    for (int i = 0; i < rows; i++) {
        Mass[i] = new double[cols];
        Momx[i] = new double[cols];
        Momy[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            
            Mass[i][j] = rho[i][j] * vol;
            Momx[i][j] = rho[i][j] * vx[i][j];
            Momy[i][j] = rho[i][j] * vy[i][j];
        }
    }

    return make_tuple(Mass, Momx, Momy);
}

tuple<double**,double**,double**,double**> getPrimitive(double** &Mass,
 double**& Momx, double**& Momy, double gamma, double vol, int rows, int cols) {
    
    double** rho = new double*[rows];
    double** vx = new double*[rows];
    double** vy = new double*[rows];
    double** P = new double*[rows];

    for (int i = 0; i < rows; i++) {
        rho[i] = new double[cols];
        vx[i] = new double[cols];
        vy[i] = new double[cols];
        P[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            
            rho[i][j] = Mass[i][j] / vol;
            vx[i][j] = Momx[i][j] / rho[i][j];
            vy[i][j] = Momy[i][j] / rho[i][j];
            P[i][j] = pow(rho[i][j], gamma);

        }
    }

    return make_tuple(rho, vx, vy, P);
}

double** getSpeedOfSound(double**& rho, double gamma, int rows, int cols) {
    double** cs = new double*[rows];
    for (int i=0; i < rows; i++){
        cs[i] = new double[cols];
        for (int j=0; j < cols; j++){
            cs[i][j] = sqrt(gamma * pow(rho[i][j], gamma - 1));
        }
    }

    return cs;
}

double minmod(double x, double y) {
    

    double sign1 = sign(x);
    double sign2 = sign(y);
    double abs1 = std::abs(x);
    double abs2 = std::abs(y);

    double mmod = (sign1 + sign2) * min(abs1, abs2) / 2;

    return mmod;
       
}

double minmod3(double x, double y, double z) {
    return minmod(x,minmod(y,z));
}

double** getGradient (double**& f, double dx, int axis, int rows, int cols, double theta = 1) {
  
    double** df_dx = new double*[rows];

    df_dx[0] = new double[cols];
    df_dx[rows-1] = new double[cols];
    
    if (axis == 0) {
        for (int i = 1; i < rows - 1; i++) {
            df_dx[i] = new double[cols];

            df_dx[i][0] = 0;
            df_dx[i][cols-1] = 0;
            df_dx[0][i] = 0;
            df_dx[rows-1][i] = 0;

            for (int j = 1; j < cols - 1; j++) {
                double term1 = theta * (f[i][j] - f[i-1][j]) / dx;
                double term2 = (f[i+1][j] - f[i-1][j]) / (2 * dx);
                double term3 = theta * (f[i+1][j] - f[i][j]) / dx;
                df_dx[i][j] = minmod3(term1, term2, term3);
            }
        }
    }
    else if (axis == 1) {
        for (int i = 1; i < rows - 1; i++) {
            df_dx[i] = new double[cols];

            df_dx[i][0] = 0;
            df_dx[i][cols-1] = 0;
            df_dx[0][i] = 0;
            df_dx[rows-1][i] = 0;

            for (int j = 1; j < cols - 1; j++) {
                double term1 = theta * (f[i][j] - f[i][j-1]) / dx;
                double term2 = (f[i][j+1] - f[i][j-1]) / (2 * dx);
                double term3 = theta * (f[i][j+1] - f[i][j]) / dx;
                df_dx[i][j] = minmod3(term1, term2, term3);    
                
            }
        }
    }

    return df_dx;
}

tuple<double**, double**, double**, double**> extrapolateInSpaceToFace (double**& q, double**& q_dx, double dx, int axis, int rows, int cols) {
    

    double** qP_XL = new double*[rows];
    double** qP_XR = new double*[rows];
    double** qM_XR = new double*[rows];
    double** qM_XL = new double*[rows];
    qP_XL[0] = new double[cols];
    qP_XR[0] = new double[cols];
    qM_XR[0] = new double[cols];
    qM_XL[0] = new double[cols];
    qP_XL[rows-1] = new double[cols];
    qP_XR[rows-1] = new double[cols];
    qM_XR[rows-1] = new double[cols];
    qM_XL[rows-1] = new double[cols];

    if (axis == 0) {
        for (int i = 1; i < rows - 1; i++) {
            qP_XL[i] = new double[cols];
            qP_XR[i] = new double[cols];
            qM_XR[i] = new double[cols];
            qM_XL[i] = new double[cols];

            qP_XL[i][0] = 0;
            qP_XR[i][cols-1] = 0;
            qM_XR[0][i] = 0;
            qM_XL[rows-1][i] = 0;

            for (int j = 1; j < cols - 1; j++) {
                qP_XL[i][j] = q[i][j] - q_dx[i][j] * dx / 2;
                qP_XR[i][j] = q[i+1][j] - q_dx[i+1][j] * dx / 2;
                qM_XR[i][j] = q[i][j] + q_dx[i][j] * dx / 2;
                qM_XL[i][j] = q[i-1][j] + q_dx[i-1][j] * dx / 2;
            }
        }
    }
    else if (axis == 1) {
        for (int i = 1; i < rows-1; i++) {
            qP_XL[i] = new double[cols];
            qP_XR[i] = new double[cols];
            qM_XR[i] = new double[cols];
            qM_XL[i] = new double[cols];

            qP_XL[i][0] = 0;
            qP_XR[i][cols-1] = 0;
            qM_XR[0][i] = 0;
            qM_XL[rows-1][i] = 0;

            for (int j = 1; j < cols-1; j++) {
                qP_XL[i][j] = q[i][j] - q_dx[i][j] * dx / 2;
                qP_XR[i][j] = q[i][j+1] - q_dx[i][j+1] * dx / 2;
                qM_XR[i][j] = q[i][j] + q_dx[i][j] * dx / 2;
                qM_XL[i][j] = q[i][j-1] + q_dx[i][j-1] * dx / 2;
            }
        }
    }

    return make_tuple(qM_XL, qP_XL, qM_XR, qP_XR);
}

double** local_propagation_speed (double**& rho, double**& vx, double**& vy, double eta, double zeta, double tau_nu, double**& cs, int rows, int cols) {
   
    double** maxC = new double*[rows];
    for (int i = 0; i < rows; i++) {
        maxC[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
                double C1, C2;
                C1 = sqrt(eta * tau_nu / rho[i][j]);
                C2 = cs[i][j] * sqrt( (zeta + 4.0 / 3.0 * tau_nu) / (rho[i][j] * tau_nu));
                maxC[i][j] = std::max(C1, C2);
        }
    }

    return maxC;
}


tuple<double**, double**, double**, double**,
 double**, double**, double**> getXFlux(double**& rho_P,
    double**& rho_M, double**& vx_P, double**& vx_M, 
    double**& vy_P, double**& vy_M, double**& Pixx_P, 
   double**& Pixx_M, double**& Pixy_P, double**& Pixy_M, 
   double**& Piyx_P, double**& Piyx_M, double**& Piyy_P, 
   double**& Piyy_M, double**& P_P, double**& P_M, double gamma, 
  double eta, double zeta, double tau_nu, int rows, int cols) {


    double** flux_Mass = new double*[rows];
    double** flux_Momx = new double*[rows];
    double** flux_Momy = new double*[rows];
    double** flux_Pixx_vx = new double*[rows];
    double** flux_Pixy_vx = new double*[rows];
    double** flux_Piyx_vx = new double*[rows];
    double** flux_Piyy_vx = new double*[rows];

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

    double** cs_P = getSpeedOfSound(rho_P, gamma, rows, cols);
    double** C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P, rows, cols);
    double** cs_M = getSpeedOfSound(rho_M, gamma, rows, cols);
    double** C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M, rows, cols);

    double C;

    double B = eta / tau_nu;
    double A = zeta / tau_nu;

    for (int i = 0; i < rows; i++) {
        flux_Mass[i] = new double[cols];
        flux_Momx[i] = new double[cols];
        flux_Momy[i] = new double[cols];
        flux_Pixx_vx[i] = new double[cols];
        flux_Pixy_vx[i] = new double[cols];
        flux_Piyx_vx[i] = new double[cols];
        flux_Piyy_vx[i] = new double[cols];

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

            C = max(C_M[i][j], C_P[i][j]);

            flux_Mass[i][j] -= C * 0.5 * (rho_P[i][j] - rho_M[i][j]);
            flux_Momx[i][j] -= C * 0.5 * (rho_P[i][j] * vx_P[i][j] - rho_M[i][j] * vx_M[i][j]);
            flux_Momy[i][j] -= C * 0.5 * (rho_P[i][j] * vy_P[i][j] - rho_M[i][j] * vy_M[i][j]);
            flux_Pixx_vx[i][j] -= C * 0.5 * (Pixx_P[i][j] - Pixx_M[i][j]);
            flux_Pixy_vx[i][j] -= C * 0.5 * (Pixy_P[i][j] - Pixy_M[i][j]);
            flux_Piyx_vx[i][j] -= C * 0.5 * (Piyx_P[i][j] - Piyx_M[i][j]);
            flux_Piyy_vx[i][j] -= C * 0.5 * (Piyy_P[i][j] - Piyy_M[i][j]);

        }
    }


    delete[] cs_M; 
    delete[] cs_P;
    delete[] C_M;
    delete[] C_P;

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vx, flux_Pixy_vx, flux_Piyx_vx, flux_Piyy_vx);
  }



tuple<double**, double**, double**, double**,
 double**, double**, double**> getYFlux(double**& rho_P, double**& rho_M,
             double**& vx_P, double**& vx_M,
             double**& vy_P, double**& vy_M,
             double**& Pixx_P, double**& Pixx_M,
             double**& Pixy_P, double**& Pixy_M,
             double**& Piyx_P, double**& Piyx_M,
             double**& Piyy_P, double**& Piyy_M,
             double**& P_P, double**& P_M,
             double gamma, double eta,
             double zeta, double tau_nu,  int rows, int cols){

    

    double** flux_Mass = new double*[rows];
    double** flux_Momx = new double*[rows];
    double** flux_Momy = new double*[rows];
    double** flux_Pixx_vy = new double*[rows];
    double** flux_Pixy_vy = new double*[rows];
    double** flux_Piyx_vy = new double*[rows];
    double** flux_Piyy_vy = new double*[rows];

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

    double** cs_P = getSpeedOfSound(rho_P, gamma, rows, cols);
    double** C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P, rows, cols);

    double** cs_M = getSpeedOfSound(rho_M, gamma, rows, cols);
    double** C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M, rows, cols);

    double C;
    double B = eta / tau_nu;
    double A = zeta/ tau_nu;


    for (int i = 0; i < rows; i++) {
        flux_Mass[i] = new double[cols];
        flux_Momx[i] = new double[cols];
        flux_Momy[i] = new double[cols];
        flux_Pixx_vy[i] = new double[cols];
        flux_Pixy_vy[i] = new double[cols];
        flux_Piyx_vy[i] = new double[cols];
        flux_Piyy_vy[i] = new double[cols];
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

            C = max(C_M[i][j], C_P[i][j]);

            flux_Mass[i][j] -= C * 0.5 * (rho_P[i][j] - rho_M[i][j]);
            flux_Momx[i][j] -= C * 0.5 * (rho_P[i][j] * vx_P[i][j] - rho_M[i][j] * vx_M[i][j]);
            flux_Momy[i][j] -= C * 0.5 * (rho_P[i][j] * vy_P[i][j] - rho_M[i][j] * vy_M[i][j]);
            flux_Pixx_vy[i][j] -= C * 0.5 * (Pixx_P[i][j] - Pixx_M[i][j]);
            flux_Pixy_vy[i][j] -= C * 0.5 * (Pixy_P[i][j] - Pixy_M[i][j]);
            flux_Piyx_vy[i][j] -= C * 0.5 * (Piyx_P[i][j] - Piyx_M[i][j]);
            flux_Piyy_vy[i][j] -= C * 0.5 * (Piyy_P[i][j] - Piyy_M[i][j]);

        }
    }

    delete[] cs_M;
    delete[] cs_P;
    delete[] C_M;
    delete[] C_P;

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vy, flux_Pixy_vy, flux_Piyx_vy, flux_Piyy_vy);
        
}
// Apply fluxes to conserved variables

double** applyFluxes(double**& flux_H1_X, double**& flux_H2_X,
  double**& flux_H1_Y, double**& flux_H2_Y,
   double dx, double dy, double**& J, int rows, int cols){


    double** C = new double*[rows];
    
    for (int i = 0; i < rows; i++) {
        C[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            C[i][j] = (flux_H1_X[i][j] - flux_H2_X[i][j]) / dx;
            C[i][j] -= (flux_H1_Y[i][j] - flux_H2_Y[i][j]) / dy;
            C[i][j] += J[i][j];
        }
    }

    return C;
}

// Heun's method
state heuns (state& q, function<state(double,state)> f, double dt, double t, int rows, int cols) {
    

    double k1,k2;
    double** c1 = new double*[rows];
    double** c2 = new double*[rows];
    double** y = new double*[rows];
    double** yprime = new double*[rows];


    state qprime;
    state C1,C2;

    C1 = f(t,q);

    for (int n = 0; n < 7; ++n){
        c1 = C1.get(n);
        y = q.get(n);
        for (int i = 0; i < rows; i++) {
            yprime[i] = new double[cols];
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
            delete[] yprime[i];
        }
        qprime.set(n , yprime);
    }
    
    delete[] c1;
    delete[] c2; 
    delete[] y;
    delete[] yprime;

    return qprime;
}


state rK4 (state& q, function<state(double,state)> f, double dt, double t, int rows, int cols) {
  

    double k1,k2,k3,k4;
    double** c1 = new double*[rows];
    double** c2 = new double*[rows];
    double** c3 = new double*[rows];
    double** c4 = new double*[rows];
    double** y = new double*[rows];
    double** yprime = new double*[rows];


    state qprime;
    state C1,C2,C3,C4;

    C1 = f(t,q);

    for (int n = 0; n < 7; ++n){
        c1 = C1.get(n);
        y = q.get(n);
        for (int i = 0; i < rows; i++) {
            yprime[i] = new double[cols];
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
                yprime[i][j] = y[i][j] + 0.5*k2;
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
            delete[] yprime[i];

        }
        qprime.set(n , yprime);
    }

    delete[] c1;
    delete[] c2;
    delete[] c3;
    delete[] c4;
    delete[] y;
    delete[] yprime;

    return qprime;
}

