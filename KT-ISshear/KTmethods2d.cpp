#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <map>
//#include <Ktmethods2d.h>
using namespace std;

class State {       
  public:             
    vector<vector<double>> rho;        
    vector<vector<double>> Momx;
    vector<vector<double>> Momy;
    vector<vector<double>> Pixx;
    vector<vector<double>> Pixy;
    vector<vector<double>> Piyx;
    vector<vector<double>> Piyy;

    State() : rho(), Momx(), Momy(), Pixx(), Pixy(), Piyx(), Piyy() {}

    State(vector<vector<double>> v1,
     vector<vector<double>> v2,
     vector<vector<double>> v3,
     vector<vector<double>> v4,
     vector<vector<double>> v5,
     vector<vector<double>> v6,
     vector<vector<double>> v7) : rho(v1), Momx(v2), Momy(v3), Pixx(v4), Pixy(v5), Piyx(v6), Piyy(v7) {}

    vector<vector<double>> get(int n) {
        if (n == 0){
            return rho;
        }
        if (n == 1){
            return Momx;
        }
        if (n == 2){
            return Momy;
        }
        if (n == 3){
            return Pixx;
        }
        if (n == 4){
            return Pixy;
        }
        if (n == 5){
            return Piyx;
        }
        if (n == 6){
            return Piyy;
        }
        else { 
            return vector<vector<double>>();
        } 
        
    }

    void set(int n , vector<vector<double>> v) {

        if (n == 0){
            rho = v;
        }
        if (n == 1){
            Momx = v;
        }
        if (n == 2){
            Momy = v;
        }
        if (n == 3){
            Pixx = v;
        }
        if (n == 4){
            Pixy = v;
        }
        if (n == 5){
            Piyx = v;
        }
        if (n == 6){
            Piyy = v;
        }
    }

    

};






double max(vector<vector<double>> value)
{  
    int rows,cols;
    rows = value.size();
    cols = value[0].size();

    double max = value[0][0];

    for (int i = 0; i < value.size(); i++){
        for (int j = 0; j < value[0].size(); j++){
            if (value[i][j] > max){
                max = value[i][j];
            }
        }
    }
    return max;

}

vector<vector<double>> sign(vector<vector<double>> value)
{   
    int rows,cols;
    rows = value.size();
    cols = value[0].size();

    vector<vector<double>> sign(rows, vector<double>(cols, 0.0));
    int rows,cols;
    rows = value.size();
    cols = value[0].size();

    vector<vector<double>> sign(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < value.size(); i++){
        for (int j = 0; j < value[0].size(); j++){
            if (value[i][j] > 0.0){
                return 1;
            } else if (value[i][j] < 0.0){
                return -1;
            } else {
                return 0;
            }
        }
    }
}

tuple<vector<vector<double>>,vector<vector<double>>,vector<vector<double>>> getConserved(vector<vector<double>>& rho,
 vector<vector<double>>& vx, vector<vector<double>>& vy, double vol) {
   
    vector<vector<double>> Mass(rho.size(), vector<double>(rho[0].size()));
    vector<vector<double>> Momx(rho.size(), vector<double>(rho[0].size()));
    vector<vector<double>> Momy(rho.size(), vector<double>(rho[0].size()));

    for (int i = 0; i < rho.size(); i++) {
        for (int j = 0; j < rho[0].size(); j++) {
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

    for (int i = 0; i < Mass.size(); i++) {
        for (int j = 0; j < Mass[0].size(); j++) {
            rho[i][j] = Mass[i][j] / vol;
            vx[i][j] = (rho[i][j] != 0) ? Momx[i][j] / rho[i][j] : 0;
            vy[i][j] = (rho[i][j] != 0) ? Momy[i][j] / rho[i][j] : 0;
            P[i][j] = pow(abs(rho[i][j]), gamma);
        }
    }

    return make_tuple(rho, vx, vy, P);
}

vector<vector<double>> getSpeedOfSound(vector<vector<double>>& rho, double gamma) {

    int rows,cols;
    rows = rho.size();
    cols = rho[0].size();

    vector<vector<double>> cs(rows, vector<double>(cols, 0.0));

    for (int i=0; i < rows; i++){
        for (int j=0; j < cols; j++){
            cs[i][j] = sqrt(gamma * pow(abs(rho[i][j]), gamma - 1));
        }
    }

    return cs;
}

vector<vector<double>> minmod2(vector<vector<double>>& x, vector<vector<double>>& y) {
    return (sign(x) + sign(y)) * min(abs(x), abs(y)) / 2;
}

vector<vector<double>> minmod3(vector<vector<double>>& x, vector<vector<double>>& y, vector<vector<double>>& z) {
    return minmod2(x, minmod2(y, z));
}

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
                double term1 = theta * (f[i][j] - f[Km1[j]][j]) / dx;
                double term2 = (f[Kp1[j]][j] - f[Km1[j]][j]) / (2 * dx);
                double term3 = theta * (f[Kp1[j]][j] - f[i][j]) / dx;
                df_dx[i][j] = min({term1, term2, term3});
            }
        }
    }
    else if (axis == 1) {
        for (int i = 0; i < f.size(); i++) {
            for (int j = 0; j < f[i].size(); j++) {
                double term1 = theta * (f[i][j] - f[i][Km1[j]]) / dx;
                double term2 = (f[i][Kp1[j]] - f[i][Km1[j]]) / (2 * dx);
                double term3 = theta * (f[i][Kp1[j]] - f[i][j]) / dx;
                df_dx[i][j] = min({term1, term2, term3});
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
                qP_XR[i][j] = q[Kp1[j]][j] - q_dx[Kp1[j]][j] * dx / 2;
                qM_XR[i][j] = q[i][j] + q_dx[i][j] * dx / 2;
                qM_XL[i][j] = q[Km1[j]][j] + q_dx[Km1[j]][j] * dx / 2;
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


vector<vector<double>> local_propagation_speed (vector<vector<double>>& rho, vector<vector<double>>& vx, vector<vector<double>>& vy, double eta, double zeta, double tau_nu, vector<vector<double>>& cs) {
   
    int rows = rho.size();
    int cols = rho[0].size();

    vector<vector<double>> C1(rows, vector<double>(cols, 0.0));
    vector<vector<double>> C2(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (rho[i][j] != 0) {
                C1[i][j] = sqrt(eta * tau_nu / rho[i][j]);
                C2[i][j] = sqrt(cs[i][j] * cs[i][j] * (zeta + 4.0 / 3.0 * tau_nu) / (rho[i][j] * tau_nu));
            }
        }
    }

    vector<vector<double>> maxC(rows, vector<double>(cols, 0.0));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            maxC[i][j] = max(C1[i][j], C2[i][j]);
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

    vector<vector<double>> rho_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> momx_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Pixx_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Piyx_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Pixx_vx_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Pixy_vx_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Piyx_vx_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Piyy_vx_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> P_av(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            rho_av[i][j] = 0.5 * (rho_P[i][j] + rho_M[i][j]);
            momx_av[i][j] = 0.5 * (rho_P[i][j] * vx_P[i][j] + rho_M[i][j] * vx_M[i][j]);
            Pixx_av[i][j] = 0.5 * (Pixx_P[i][j] + Pixx_M[i][j]);
            Piyx_av[i][j] = 0.5 * (Piyx_P[i][j] + Piyx_M[i][j]);
            Pixx_vx_av[i][j] = 0.5 * (Pixx_P[i][j] * vx_P[i][j] + Pixx_M[i][j] * vx_M[i][j]);
            Pixy_vx_av[i][j] = 0.5 * (Pixy_P[i][j] * vx_P[i][j] + Piyx_M[i][j] * vx_M[i][j]);
            Piyx_vx_av[i][j] = 0.5 * (Piyx_P[i][j] * vx_P[i][j] + Pixy_M[i][j] * vx_M[i][j]);
            Piyy_vx_av[i][j] = 0.5 * (Piyy_P[i][j] * vx_P[i][j] + Piyy_M[i][j] * vx_M[i][j]);
            P_av[i][j] = 0.5 * (P_P[i][j] + P_M[i][j]);

            double B = eta / tau_nu;
            double A = zeta / tau_nu;

            flux_Mass[i][j] = momx_av[i][j];
            flux_Momx[i][j] = 0.5 * (rho_P[i][j] * pow(vx_P[i][j], 2) + rho_M[i][j] * pow(vx_M[i][j], 2)) + (P_av[i][j]) + (Pixx_av[i][j]) / gamma;
            flux_Momy[i][j] = 0.5 * (rho_P[i][j] * (vx_P[i][j] * vy_P[i][j]) + rho_M[i][j] * (vx_M[i][j] * vy_M[i][j])) + (Piyx_av[i][j]) / gamma;
            flux_Pixx_vx[i][j] = Pixx_vx_av[i][j] + B * (vx_P[i][j] + vx_M[i][j]) + (A - 2.0 / 3.0 * B) * (vx_P[i][j] + vx_M[i][j]) * 0.5;
            flux_Pixy_vx[i][j] = Pixy_vx_av[i][j] + B * (vy_P[i][j] + vy_M[i][j]) * 0.5;
            flux_Piyx_vx[i][j] = Piyx_vx_av[i][j] + B * (vy_P[i][j] + vy_M[i][j]) * 0.5;
            flux_Piyy_vx[i][j] = Piyy_vx_av[i][j] + (A - 2.0 / 3.0 * B) * (vx_P[i][j] + vx_M[i][j]) * 0.5;
            }
            }

    vector<vector<double>> cs_P = getSpeedOfSound(rho_P, gamma);
    vector<vector<double>> C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P);
    vector<vector<double>> cs_M = getSpeedOfSound(rho_M, gamma);
    vector<vector<double>> C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M);

    vector<vector<double>> C(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            C[i][j] = max(C_M[i][j], C_P[i][j]);
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            flux_Mass[i][j] -= C[i][j] * 0.5 * (rho_P[i][j] - rho_M[i][j]);
            flux_Momx[i][j] -= C[i][j] * 0.5 * (rho_P[i][j] * vx_P[i][j] - rho_M[i][j] * vx_M[i][j]);
            flux_Momy[i][j] -= C[i][j] * 0.5 * (rho_P[i][j] * vy_P[i][j] - rho_M[i][j] * vy_M[i][j]);
            flux_Pixx_vx[i][j] -= C[i][j] * 0.5 * (Pixx_P[i][j] - Pixx_M[i][j]);
            flux_Pixy_vx[i][j] -= C[i][j] * 0.5 * (Pixy_P[i][j] - Pixy_M[i][j]);
            flux_Piyx_vx[i][j] -= C[i][j] * 0.5 * (Piyx_P[i][j] - Piyx_M[i][j]);
            flux_Piyy_vx[i][j] -= C[i][j] * 0.5 * (Piyy_P[i][j] - Piyy_M[i][j]);
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

    vector<vector<double>> rho_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> momy_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Piyy_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Pixy_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Pixx_vy_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Pixy_vy_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Piyx_vy_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> Piyy_vy_av(rows, vector<double>(cols, 0.0));
    vector<vector<double>> P_av(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            rho_av[i][j] = 0.5 * (rho_P[i][j] + rho_M[i][j]);
            momy_av[i][j] = 0.5 * (rho_P[i][j] * vy_P[i][j] + rho_M[i][j] * vy_M[i][j]);
            Piyy_av[i][j] = 0.5 * (Piyy_P[i][j] + Piyy_M[i][j]);
            Pixy_av[i][j] = 0.5 * (Pixy_P[i][j] + Pixy_M[i][j]);
            Pixx_vy_av[i][j] = 0.5 * (Pixx_P[i][j] * vy_P[i][j] + Pixx_M[i][j] * vy_M[i][j]);
            Pixy_vy_av[i][j] = 0.5 * (Pixy_P[i][j] * vy_P[i][j] + Piyx_M[i][j] * vy_M[i][j]);
            Piyx_vy_av[i][j] = 0.5 * (Piyx_P[i][j] * vy_P[i][j] + Pixy_M[i][j] * vy_M[i][j]);
            Piyy_vy_av[i][j] = 0.5 * (Piyy_P[i][j] * vy_P[i][j] + Piyy_M[i][j] * vy_M[i][j]);
            P_av[i][j] = 0.5 * (P_P[i][j] + P_M[i][j]);
        }
    }

    vector<vector<double>> cs_P = getSpeedOfSound(rho_P, gamma);
    vector<vector<double>> C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P);

    vector<vector<double>> cs_M = getSpeedOfSound(rho_M, gamma);
    vector<vector<double>> C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M);

    vector<vector<double>> C(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            C[i][j] = max(C_M[i][j], C_P[i][j]);
        }
    }

    double B = eta / tau_nu;
    double A = zeta/ tau_nu;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            flux_Mass[i][j] = momy_av[i][j];
            flux_Momx[i][j] = 0.5 * (rho_P[i][j] * vx_P[i][j] * vy_P[i][j] + rho_M[i][j] * vx_M[i][j] * vy_M[i][j]) + Pixy_av[i][j] / gamma;
            flux_Momy[i][j] = 0.5 * (rho_P[i][j] * vy_P[i][j] * vy_P[i][j] + rho_M[i][j] * vy_M[i][j] * vy_M[i][j]) + (P_av[i][j]) + (Piyy_av[i][j]) / gamma;
            flux_Pixx_vy[i][j] = Pixx_vy_av[i][j] + (A - 2.0 / 3.0 * B) * (vx_P[i][j] + vx_M[i][j]) * 0.5;
            flux_Pixy_vy[i][j] = Pixy_vy_av[i][j] + B * (vy_P[i][j] + vy_M[i][j]) * 0.5;
            flux_Piyx_vy[i][j] = Piyx_vy_av[i][j] + B * (vy_P[i][j] + vy_M[i][j]) * 0.5;
            flux_Piyy_vy[i][j] = Piyy_vy_av[i][j] + B * (vx_P[i][j] + vx_M[i][j]) + (A - 2.0 / 3.0 * B) * (vx_P[i][j] + vx_M[i][j]) * 0.5;
            flux_Mass[i][j] -= C[i][j] * 0.5 * (rho_P[i][j] - rho_M[i][j]);
            flux_Momx[i][j] -= C[i][j] * 0.5 * (rho_P[i][j] * vx_P[i][j] - rho_M[i][j] * vx_M[i][j]);
            flux_Momy[i][j] -= C[i][j] * 0.5 * (rho_P[i][j] * vy_P[i][j] - rho_M[i][j] * vy_M[i][j]);
            flux_Pixx_vy[i][j] -= C[i][j] * 0.5 * (Pixx_P[i][j] - Pixx_M[i][j]);
            flux_Pixy_vy[i][j] -= C[i][j] * 0.5 * (Pixy_P[i][j] - Pixy_M[i][j]);
            flux_Piyx_vy[i][j] -= C[i][j] * 0.5 * (Piyx_P[i][j] - Piyx_M[i][j]);
            flux_Piyy_vy[i][j] -= C[i][j] * 0.5 * (Piyy_P[i][j] - Piyy_M[i][j]);
            }
        }
        
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
State Heuns (
State& q, 
function<State(double,State)> f, double dt, double t) {
    
    int rows = (q.get(0)).size();
    int cols = (q.get(0))[0].size();
    int rows = (q.get(0)).size();
    int cols = (q.get(0))[0].size();

    vector<vector<double>> k1(rows, vector<double>(cols, 0.0));
    vector<vector<double>> k2(rows, vector<double>(cols, 0.0));
    vector<vector<double>> c(rows, vector<double>(cols, 0.0));
    vector<vector<double>> y(rows, vector<double>(cols, 0.0));
    vector<vector<double>> yprime(rows, vector<double>(cols, 0.0));


    State qprime;
    State C;

    C = f(t,q);
    for (int n = 0; n < tuple_size<decltype(C)>::value; ++n){
        c = C.get(n);
        y = q.get(n);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                k1[i][j] = dt * c[i][j];
                yprime[i][j] = y[i][j] + k1[i][j];
            }
        }
        tuple_cat( qprime, make_tuple(yprime) );
    }

    C = f(t,qprime);


    return qprime;
}
