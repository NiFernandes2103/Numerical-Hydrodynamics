#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <map>
#include <functional>
#include "KTmethods2d.h"
using namespace std;

double max_value(vector<vector<double>> value)
{  
    unsigned int rows,cols;
    rows = value.size();
    cols = value[0].size();

    double max = value[0][0];

    for (unsigned int i = 0; i < rows; i++){
        for (unsigned int j = 0; j < cols; j++){
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



tuple<vector<vector<double>>,vector<vector<double>>,vector<vector<double>>> getConserved(vector<vector<double>>& rho,
 vector<vector<double>>& vx, vector<vector<double>>& vy, double vol) {
   
    vector<vector<double>> Mass(rho.size(), vector<double>(rho[0].size()));
    vector<vector<double>> Momx(rho.size(), vector<double>(rho[0].size()));
    vector<vector<double>> Momy(rho.size(), vector<double>(rho[0].size()));

    for (unsigned int i = 0; i < rho.size(); i++) {
        for (unsigned int j = 0; j < rho[0].size(); j++) {
            Mass[i][j] = rho[i][j] * vol;
            Momx[i][j] = rho[i][j] * vx[i][j];
            Momy[i][j] = rho[i][j] * vy[i][j];
        }
    }

    return make_tuple(Mass, Momx, Momy);
}

tuple<vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,vector<vector<double>>> getPrimitive(vector<vector<double>> &Mass,
 vector<vector<double>>& Momx, vector<vector<double>>& Momy, double gamma, double vol) {
    
    vector<vector<double>> rho(Mass.size(), vector<double>(Mass[0].size()));
    vector<vector<double>> vx(Mass.size(), vector<double>(Mass[0].size()));
    vector<vector<double>> vy(Mass.size(), vector<double>(Mass[0].size()));
    vector<vector<double>> P(Mass.size(), vector<double>(Mass[0].size()));

    for (unsigned int i = 0; i < Mass.size(); i++) {
        for (unsigned int j = 0; j < Mass[0].size(); j++) {
            rho[i][j] = Mass[i][j] / vol;
            vx[i][j] = Momx[i][j] / rho[i][j];
            vy[i][j] = Momy[i][j] / rho[i][j];
            P[i][j] = pow(rho[i][j], gamma);

            
            // set boundary conditions
            rho[0][j] = rho[1][j];
            rho[i][0] = rho[i][1];
            rho[-1][j] = rho[-2][j];
            rho[i][-1] = rho[i][-2];

            vx[0][j] = -vx[1][j];
            vx[i][0] = 0;
            vx[-1][j] = -vx[-2][j];
            vx[i][-1] = 0;

            vy[0][j] = 0;
            vy[i][0] = -vy[i][1];
            vy[-1][j] = 0;
            vy[i][-1] = -vy[i][-2];

            P[0][j] = P[1][j];
            P[i][0] = P[i][1];
            P[-1][j] = P[-2][j];
            P[i][-1] = P[i][-2];

        }
    }

    return make_tuple(rho, vx, vy, P);
}

vector<vector<double>> getSpeedOfSound(vector<vector<double>>& rho, double gamma) {

    unsigned int rows,cols;
    rows = rho.size();
    cols = rho[0].size();

    vector<vector<double>> cs(rows, vector<double>(cols, 0.0));

    for (unsigned int i=0; i < rows; i++){
        for (unsigned int j=0; j < cols; j++){
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

/*
vector<vector<double>> getGradient (vector<vector<double>>& f, double dx, int axis, double theta = 1) {
    vector<vector<double>> df_dx(f.size(), vector<double>(f[0].size(), 0.0));
    int n = f[axis].size();

    vector<int> K(n);
    for (int i = 0; i < n; i++) {
        K[i] = i;
    }

    vector<int> Kp1(K);
    rotate(Kp1.begin(), Kp1.begin() + 1, Kp1.end());

    vector<int> Km1(K);
    rotate(Km1.rbegin(), Km1.rbegin() + 1, Km1.rend());


    if (axis == 0) {
        for (int i = 0; i < f.size(); i++) {
            for (int j = 0; j < f[i].size(); j++) {
                double term1 = theta * (f[i][j] - f[Km1[i]][j]) / dx;
                double term2 = (f[Kp1[i]][j] - f[Km1[i]][j]) / (2 * dx);
                double term3 = theta * (f[Kp1[i]][j] - f[i][j]) / dx;
                df_dx[i][j] = minmod3(term1, term2, term3);
            }
        }
    }
    else if (axis == 1) {
        for (int i = 0; i < f.size(); i++) {
            for (int j = 0; j < f[i].size(); j++) {
                double term1 = theta * (f[i][j] - f[i][Km1[j]]) / dx;
                double term2 = (f[i][Kp1[j]] - f[i][Km1[j]]) / (2 * dx);
                double term3 = theta * (f[i][Kp1[j]] - f[i][j]) / dx;
                df_dx[i][j] = minmod3(term1, term2, term3);
            }
        }
    }


    return df_dx;
}

tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> extrapolateInSpaceToFace (vector<vector<double>>& q, vector<vector<double>>& q_dx, double dx, int axis) {

    int n = q.size();

    vector<int> K(n);
    for (int i = 0; i < n; i++) {
        K[i] = i;
    }

    vector<int> Kp1(K);
    rotate(Kp1.begin(), Kp1.begin() + 1, Kp1.end());

    vector<int> Km1(K);
    rotate(Km1.rbegin(), Km1.rbegin() + 1, Km1.rend());

    vector<vector<double>> qP_XL(q.size(), vector<double>(q[0].size(), 0.0));
    vector<vector<double>> qP_XR(q.size(), vector<double>(q[0].size(), 0.0));
    vector<vector<double>> qM_XR(q.size(), vector<double>(q[0].size(), 0.0));
    vector<vector<double>> qM_XL(q.size(), vector<double>(q[0].size(), 0.0));

    if (axis == 0) {
        for (int i = 0; i < q.size(); i++) {
            for (int j = 0; j < q[i].size(); j++) {
                qP_XL[i][j] = q[i][j] - q_dx[i][j] * dx / 2;
                qP_XR[i][j] = q[Kp1[i]][j] - q_dx[Kp1[i]][j] * dx / 2;
                qM_XR[i][j] = q[i][j] + q_dx[i][j] * dx / 2;
                qM_XL[i][j] = q[Km1[i]][j] + q_dx[Km1[i]][j] * dx / 2;
            }
        }
    }
    else if (axis == 1) {
        for (int i = 0; i < q.size(); i++) {
            for (int j = 0; j < q[i].size(); j++) {
                qP_XL[i][j] = q[i][j] - q_dx[i][j] * dx / 2;
                qP_XR[i][j] = q[i][Kp1[j]] - q_dx[i][Kp1[j]] * dx / 2;
                qM_XR[i][j] = q[i][j] + q_dx[i][j] * dx / 2;
                qM_XL[i][j] = q[i][Km1[j]] + q_dx[i][Km1[j]] * dx / 2;
            }
        }
    }

    return make_tuple(qM_XL, qP_XL, qM_XR, qP_XR);
}
*/
vector<vector<double>> getGradient (vector<vector<double>>& f, double dx, int axis, double theta = 1) {
    int n = f.size();
    int m = f[0].size();

    vector<vector<double>> df_dx(n, vector<double>(m, 0.0));
    

    if (axis == 0) {
        for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < m - 1; j++) {
                double term1 = theta * (f[i][j] - f[i-1][j]) / dx;
                double term2 = (f[i+1][j] - f[i-1][j]) / (2 * dx);
                double term3 = theta * (f[i+1][j] - f[i][j]) / dx;
                df_dx[i][j] = minmod3(term1, term2, term3);
            }
        }
    }
    else if (axis == 1) {
        for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < m - 1; j++) {
                double term1 = theta * (f[i][j] - f[i][j-1]) / dx;
                double term2 = (f[i][j+1] - f[i][j-1]) / (2 * dx);
                double term3 = theta * (f[i][j+1] - f[i][j]) / dx;
                df_dx[i][j] = minmod3(term1, term2, term3);
            }
        }
    }

    return df_dx;
}

tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> extrapolateInSpaceToFace (vector<vector<double>>& q, vector<vector<double>>& q_dx, double dx, int axis) {
    
    int n = q.size();
    int m = q[0].size();

    vector<vector<double>> qP_XL(q.size(), vector<double>(q[0].size(), 0.0));
    vector<vector<double>> qP_XR(q.size(), vector<double>(q[0].size(), 0.0));
    vector<vector<double>> qM_XR(q.size(), vector<double>(q[0].size(), 0.0));
    vector<vector<double>> qM_XL(q.size(), vector<double>(q[0].size(), 0.0));

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

vector<vector<double>> local_propagation_speed (vector<vector<double>>& rho, vector<vector<double>>& vx, vector<vector<double>>& vy, double eta, double zeta, double tau_nu, vector<vector<double>>& cs) {
   
    int rows = rho.size();
    int cols = rho[0].size();

    vector<vector<double>> C1(rows, vector<double>(cols, 0.0));
    vector<vector<double>> C2(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (rho[i][j] != 0) {
                C1[i][j] = sqrt(eta * tau_nu / rho[i][j]);
                C2[i][j] = cs[i][j] * sqrt( (zeta + 4.0 / 3.0 * tau_nu) / (rho[i][j] * tau_nu));
            }
        }
    }

    vector<vector<double>> maxC(rows, vector<double>(cols, 0.0));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            maxC[i][j] = std::max(C1[i][j], C2[i][j]);
        }
    }

    return maxC;
}


tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>,
 vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> getXFlux(vector<vector<double>>& rho_P,
    vector<vector<double>>& rho_M, vector<vector<double>>& vx_P, vector<vector<double>>& vx_M, 
    vector<vector<double>>& vy_P, vector<vector<double>>& vy_M, vector<vector<double>>& Pixx_P, 
   vector<vector<double>>& Pixx_M, vector<vector<double>>& Pixy_P, vector<vector<double>>& Pixy_M, 
   vector<vector<double>>& Piyx_P, vector<vector<double>>& Piyx_M, vector<vector<double>>& Piyy_P, 
   vector<vector<double>>& Piyy_M, vector<vector<double>>& P_P, vector<vector<double>>& P_M, double gamma, 
  double eta, double zeta, double tau_nu) {

    int rows = rho_P.size();
    int cols = rho_P[0].size();

    vector<vector<double>> flux_Mass(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Momx(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Momy(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Pixx_vx(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Pixy_vx(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Piyx_vx(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Piyy_vx(rows, vector<double>(cols, 0.0));

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

    vector<vector<double>> cs_P = getSpeedOfSound(rho_P, gamma);
    vector<vector<double>> C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P);
    vector<vector<double>> cs_M = getSpeedOfSound(rho_M, gamma);
    vector<vector<double>> C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M);

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

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vx, flux_Pixy_vx, flux_Piyx_vx, flux_Piyy_vx);
  }



tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>,
 vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> getYFlux(vector<vector<double>>& rho_P, vector<vector<double>>& rho_M,
             vector<vector<double>>& vx_P, vector<vector<double>>& vx_M,
             vector<vector<double>>& vy_P, vector<vector<double>>& vy_M,
             vector<vector<double>>& Pixx_P, vector<vector<double>>& Pixx_M,
             vector<vector<double>>& Pixy_P, vector<vector<double>>& Pixy_M,
             vector<vector<double>>& Piyx_P, vector<vector<double>>& Piyx_M,
             vector<vector<double>>& Piyy_P, vector<vector<double>>& Piyy_M,
             vector<vector<double>>& P_P, vector<vector<double>>& P_M,
             double gamma, double eta,
             double zeta, double tau_nu){

    int rows = rho_P.size();
    int cols = rho_P[0].size();

    vector<vector<double>> flux_Mass(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Momx(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Momy(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Pixx_vy(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Pixy_vy(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Piyx_vy(rows, vector<double>(cols, 0.0));
    vector<vector<double>> flux_Piyy_vy(rows, vector<double>(cols, 0.0));

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

    vector<vector<double>> cs_P = getSpeedOfSound(rho_P, gamma);
    vector<vector<double>> C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P);

    vector<vector<double>> cs_M = getSpeedOfSound(rho_M, gamma);
    vector<vector<double>> C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M);

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
            flux_Pixx_vy[i][j] = Pixx_vy_av + (A - 2.0 / 3.0 * B) * (vy_P[i][j] + vy_M[i][j]) * 0.5;
            flux_Pixy_vy[i][j] = Pixy_vy_av + B * (vx_P[i][j] + vx_M[i][j]) * 0.5;
            flux_Piyx_vy[i][j] = Piyx_vy_av + B * (vx_P[i][j] + vx_M[i][j]) * 0.5;
            flux_Piyy_vy[i][j] = Piyy_vy_av + B * (vy_P[i][j] + vy_M[i][j]) + (A - 2.0 / 3.0 * B) * (vy_P[i][j] + vy_M[i][j]) * 0.5;

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

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vy, flux_Pixy_vy, flux_Piyx_vy, flux_Piyy_vy);
        
}
// Apply fluxes to conserved variables

vector<vector<double>> applyFluxes(vector<vector<double>>& flux_H1_X, vector<vector<double>>& flux_H2_X,
  vector<vector<double>>& flux_H1_Y, vector<vector<double>>& flux_H2_Y,
   double dx, double dy, vector<vector<double>>& J){

    int rows = flux_H1_X.size();
    int cols = flux_H1_X[0].size();


    vector<vector<double>> C(rows, vector<double>(cols, 0.0));
    
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
state heuns (state& q, function<state(double,state)> f, double dt, double t) {
    
    int rows = (q.get(0)).size();
    int cols = (q.get(0))[0].size();

    double k1,k2;
    vector<vector<double>> c1(rows, vector<double>(cols, 0.0));
    vector<vector<double>> c2(rows, vector<double>(cols, 0.0));
    vector<vector<double>> y(rows, vector<double>(cols, 0.0));
    vector<vector<double>> yprime(rows, vector<double>(cols, 0.0));


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


state rK4 (state& q, function<state(double,state)> f, double dt, double t) {
    
    int rows = (q.get(0)).size();
    int cols = (q.get(0))[0].size();

    double k1,k2,k3,k4;
    vector<vector<double>> c1(rows, vector<double>(cols, 0.0));
    vector<vector<double>> c2(rows, vector<double>(cols, 0.0));
    vector<vector<double>> c3(rows, vector<double>(cols, 0.0));
    vector<vector<double>> c4(rows, vector<double>(cols, 0.0));
    vector<vector<double>> y(rows, vector<double>(cols, 0.0));
    vector<vector<double>> yprime(rows, vector<double>(cols, 0.0));


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
        }
        qprime.set(n , yprime);
    }

    return qprime;
}
