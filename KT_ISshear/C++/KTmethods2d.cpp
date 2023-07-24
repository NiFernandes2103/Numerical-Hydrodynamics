#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <map>
#include <functional>
#include "matrix.h"
#include "KTmethods2d.h"
using namespace std;

double max_value(smatrix value)
{  
    int N;
    N = value.N;


    double max = value.get(0,0);

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (value.get(i,j) > max){
                max = value.get(i,j);
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

double absolute(double value) {

    return sign(value)*value;
}


tuple<smatrix,smatrix,smatrix> getConserved(smatrix& rho,
 smatrix& vx, smatrix& vy, double vol) {

    int N = rho.N;
   
    smatrix Mass(N);
    smatrix Momx(N);
    smatrix Momy(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Mass.set(rho.get(i,j) * vol,i,j);
            Momx.set(rho.get(i,j) * vx.get(i,j),i,j);
            Momy.set(rho.get(i,j) * vy.get(i,j),i,j);
        }
    }

    return make_tuple(Mass, Momx, Momy);
}

tuple<smatrix,smatrix,smatrix,smatrix> getPrimitive(smatrix &Mass,
 smatrix& Momx, smatrix& Momy, double gamma, double vol) {

    int N = Mass.N;
    
    smatrix rho(N);
    smatrix vx(N);
    smatrix vy(N);
    smatrix P(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            rho.set(Mass.get(i,j) / vol,i,j);
            vx.set(Momx.get(i,j) / rho.get(i,j),i,j);
            vy.set(Momy.get(i,j) / rho.get(i,j),i,j);
            P.set(pow(rho.get(i,j), gamma),i,j);
        }
    }

    return make_tuple(rho, vx, vy, P);
}

smatrix getSpeedOfSound(smatrix& rho, double gamma) {

    int N = rho.N;

    smatrix cs(N);

    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){
            cs.set(sqrt(gamma * pow(rho.get(i,j), gamma - 1)),i,j);
        }
    }

    return cs;
}

double minmod(double a, double b) {
    

    double sign1 = sign(a);
    double sign2 = sign(b);
    double abs1 = absolute(a);
    double abs2 = absolute(b);

    double mmod = (sign1 + sign2) * min(abs1, abs2) / 2;

    return mmod;
       
}

double minmod3(double a, double b, double c) {
    return minmod(a,minmod(b,c));
}

/*

smatrix getGradient (smatrix& f, double dx, int axis, double theta = 1) {
    smatrix df_dx(f.size(), vector<double>(f[0].size(), 0.0));
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
                double term1 = theta * (f.get(i,j) - f[Km1[i]][j]) / dx;
                double term2 = (f[Kp1[i]][j] - f[Km1[i]][j]) / (2 * dx);
                double term3 = theta * (f[Kp1[i]][j] - f.get(i,j)) / dx;
                df_dx.set(,i,j) minmod3(term1, term2, term3);
            }
        }
    }
    else if (axis == 1) {
        for (int i = 0; i < f.size(); i++) {
            for (int j = 0; j < f[i].size(); j++) {
                double term1 = theta * (f.get(i,j) - f[i][Km1[j]]) / dx;
                double term2 = (f[i][Kp1[j]] - f[i][Km1[j]]) / (2 * dx);
                double term3 = theta * (f[i][Kp1[j]] - f.get(i,j)) / dx;
                df_dx.set(,i,j) minmod3(term1, term2, term3);
            }
        }
    }


    return df_dx;
}
tuple<smatrix, smatrix, smatrix, smatrix> extrapolateInSpaceToFace (smatrix& q, smatrix& q_dx, double dx, int axis) {

    int n = q.size();

    vector<int> K(n);
    for (int i = 0; i < n; i++) {
        K[i] = i;
    }

    vector<int> Kp1(K);
    rotate(Kp1.begin(), Kp1.begin() + 1, Kp1.end());

    vector<int> Km1(K);
    rotate(Km1.rbegin(), Km1.rbegin() + 1, Km1.rend());

    smatrix qP_XL(N);
    smatrix qP_XR(N);
    smatrix qM_XR(N);
    smatrix qM_XL(N);

    if (axis == 0) {
        for (int i = 0; i < q.size(); i++) {
            for (int j = 0; j < q[i].size(); j++) {
                qP_XL.set(,i,j) q.get(i,j) - q_dx.get(i,j) * dx / 2;
                qP_XR.set(,i,j) q[Kp1[i]][j] - q_dx[Kp1[i]][j] * dx / 2;
                qM_XR.set(,i,j) q.get(i,j) + q_dx.get(i,j) * dx / 2;
                qM_XL.set(,i,j) q[Km1[i]][j] + q_dx[Km1[i]][j] * dx / 2;
            }
        }
    }
    else if (axis == 1) {
        for (int i = 0; i < q.size(); i++) {
            for (int j = 0; j < q[i].size(); j++) {
                qP_XL.set(,i,j) q.get(i,j) - q_dx.get(i,j) * dx / 2;
                qP_XR.set(,i,j) q[i][Kp1[j]] - q_dx[i][Kp1[j]] * dx / 2;
                qM_XR.set(,i,j) q.get(i,j) + q_dx.get(i,j) * dx / 2;
                qM_XL.set(,i,j) q[i][Km1[j]] + q_dx[i][Km1[j]] * dx / 2;
            }
        }
    }

    return make_tuple(qM_XL, qP_XL, qM_XR, qP_XR);
}
*/

smatrix getGradient (smatrix& f, double dx, int axis, double theta = 1) {
    
    int N = f.N;
   
    smatrix df_dx(N);
    

    if (axis == 0) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                double term1 = theta * (f.get(i,j) - f.get(i-1,j)) / dx;
                double term2 = (f.get(i+1,j) - f.get(i-1,j)) / (2 * dx);
                double term3 = theta * (f.get(i+1,j) - f.get(i,j)) / dx;
                df_dx.set(minmod3(term1, term2, term3),i,j);
            }
        }
    }
    else if (axis == 1) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                double term1 = theta * (f.get(i,j) - f.get(i,j-1)) / dx;
                double term2 = (f.get(i, j+1) - f.get(i,j-1)) / (2 * dx);
                double term3 = theta * (f.get(i, j+1) - f.get(i,j)) / dx;
                df_dx.set(minmod3(term1, term2, term3),i,j);
            }
        }
    }

    return df_dx;
}

tuple<smatrix, smatrix, smatrix, smatrix> extrapolateInSpaceToFace (smatrix& q, smatrix& q_dx, double dx, int axis) {
    
    int N = q.N;

    smatrix qP_XL(N);
    smatrix qP_XR(N);
    smatrix qM_XR(N);
    smatrix qM_XL(N);

    if (axis == 0) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                qP_XL.set(q.get(i,j) - q_dx.get(i,j) * dx / 2,i,j);
                qP_XR.set(q.get(i+1,j) - q_dx.get(i+1,j) * dx / 2,i,j);
                qM_XR.set(q.get(i,j) + q_dx.get(i,j) * dx / 2,i,j);
                qM_XL.set(q.get(i-1,j) + q_dx.get(i-1,j) * dx / 2,i,j);
            }
        }
    }
    else if (axis == 1) {
        for (int i = 1; i < N-1; i++) {
            for (int j = 1; j < N-1; j++) {
                qP_XL.set(q.get(i,j) - q_dx.get(i,j) * dx / 2,i,j);
                qP_XR.set(q.get(i, j+1) - q_dx.get(i, j+1) * dx / 2,i,j);
                qM_XR.set(q.get(i,j) + q_dx.get(i,j) * dx / 2,i,j);
                qM_XL.set( q.get(i,j-1) + q_dx.get(i,j-1) * dx / 2,i,j);
            }
        }
    }

    return make_tuple(qM_XL, qP_XL, qM_XR, qP_XR);
}

smatrix local_propagation_speed (smatrix& rho, smatrix& vx, smatrix& vy, double eta, double zeta, double tau_nu, smatrix& cs) {
   
    int N = rho.N;
    

    smatrix C1(N);
    smatrix C2(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (rho.get(i,j) != 0) {
                C1.set(sqrt(eta * tau_nu / rho.get(i,j)),i,j);
                C2.set(sqrt(cs.get(i,j) * cs.get(i,j) * (zeta + 4.0 / 3.0 * tau_nu) / (rho.get(i,j) * tau_nu)),i,j);
            }
        }
    }


    smatrix maxC(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            maxC.set(max(C1.get(i,j), C2.get(i,j)),i,j);
        }
    }

    return maxC;
}


tuple<smatrix, smatrix, smatrix, smatrix,
 smatrix, smatrix, smatrix> getXFlux(smatrix& rho_P,
    smatrix& rho_M, smatrix& vx_P, smatrix& vx_M, 
    smatrix& vy_P, smatrix& vy_M, smatrix& Pixx_P, 
   smatrix& Pixx_M, smatrix& Pixy_P, smatrix& Pixy_M, 
   smatrix& Piyx_P, smatrix& Piyx_M, smatrix& Piyy_P, 
   smatrix& Piyy_M, smatrix& P_P, smatrix& P_M, double gamma, 
  double eta, double zeta, double tau_nu) {

    int N = rho_P.N;

    smatrix flux_Mass(N);
    smatrix flux_Momx(N);
    smatrix flux_Momy(N);
    smatrix flux_Pixx_vx(N);
    smatrix flux_Pixy_vx(N);
    smatrix flux_Piyx_vx(N);
    smatrix flux_Piyy_vx(N);

    double momx_av;
    double Pixx_av;
    double Piyx_av;
    double Pixx_vx_av;
    double Pixy_vx_av;
    double Piyx_vx_av;
    double Piyy_vx_av;
    double P_av;

    smatrix cs_P = getSpeedOfSound(rho_P, gamma);
    smatrix C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P);
    smatrix cs_M = getSpeedOfSound(rho_M, gamma);
    smatrix C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M);

    double C;

    double B = eta / tau_nu;
    double A = zeta / tau_nu;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            momx_av = 0.5 * (rho_P.get(i,j) * vx_P.get(i,j) + rho_M.get(i,j) * vx_M.get(i,j));
            Pixx_av = 0.5 * (Pixx_P.get(i,j) + Pixx_M.get(i,j));
            Piyx_av = 0.5 * (Piyx_P.get(i,j) + Piyx_M.get(i,j));
            Pixx_vx_av = 0.5 * (Pixx_P.get(i,j) * vx_P.get(i,j) + Pixx_M.get(i,j) * vx_M.get(i,j));
            Pixy_vx_av = 0.5 * (Pixy_P.get(i,j) * vx_P.get(i,j) + Piyx_M.get(i,j) * vx_M.get(i,j));
            Piyx_vx_av = 0.5 * (Piyx_P.get(i,j) * vx_P.get(i,j) + Pixy_M.get(i,j) * vx_M.get(i,j));
            Piyy_vx_av = 0.5 * (Piyy_P.get(i,j) * vx_P.get(i,j) + Piyy_M.get(i,j) * vx_M.get(i,j));
            P_av = 0.5 * (P_P.get(i,j) + P_M.get(i,j));

            C = std::max(C_M.get(i,j), C_P.get(i,j));


            flux_Mass.set(momx_av - C * 0.5 * (rho_P.get(i,j) - rho_M.get(i,j)),i,j);
            flux_Momx.set(0.5 * (rho_P.get(i,j) * pow(vx_P.get(i,j), 2) + rho_M.get(i,j) * pow(vx_M.get(i,j), 2)) + ((P_av) + (Pixx_av)) / gamma - C * 0.5 * (rho_P.get(i,j) * vx_P.get(i,j) - rho_M.get(i,j) * vx_M.get(i,j)),i,j);
            flux_Momy.set(0.5 * (rho_P.get(i,j) * (vx_P.get(i,j) * vy_P.get(i,j)) + rho_M.get(i,j) * (vx_M.get(i,j) * vy_M.get(i,j))) + (Piyx_av) / gamma - C * 0.5 * (rho_P.get(i,j) * vy_P.get(i,j) - rho_M.get(i,j) * vy_M.get(i,j)),i,j);
            flux_Pixx_vx.set(Pixx_vx_av + B * (vx_P.get(i,j) + vx_M.get(i,j)) + (A - 2.0 / 3.0 * B) * (vx_P.get(i,j) + vx_M.get(i,j)) * 0.5 - C * 0.5 * (Pixx_P.get(i,j) - Pixx_M.get(i,j)),i,j);
            flux_Pixy_vx.set(Pixy_vx_av + B * (vy_P.get(i,j) + vy_M.get(i,j)) * 0.5 - C * 0.5 * (Pixy_P.get(i,j) - Pixy_M.get(i,j)),i,j);
            flux_Piyx_vx.set(Piyx_vx_av + B * (vy_P.get(i,j) + vy_M.get(i,j)) * 0.5 - C * 0.5 * (Piyx_P.get(i,j) - Piyx_M.get(i,j)),i,j);
            flux_Piyy_vx.set(Piyy_vx_av + (A - 2.0 / 3.0 * B) * (vx_P.get(i,j) + vx_M.get(i,j)) * 0.5 - C * 0.5 * (Piyy_P.get(i,j) - Piyy_M.get(i,j)),i,j);

        }
    }

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vx, flux_Pixy_vx, flux_Piyx_vx, flux_Piyy_vx);
  }

tuple<smatrix, smatrix, smatrix, smatrix,
 smatrix, smatrix, smatrix> getYFlux(smatrix& rho_P, smatrix& rho_M,
             smatrix& vx_P, smatrix& vx_M,
             smatrix& vy_P, smatrix& vy_M,
             smatrix& Pixx_P, smatrix& Pixx_M,
             smatrix& Pixy_P, smatrix& Pixy_M,
             smatrix& Piyx_P, smatrix& Piyx_M,
             smatrix& Piyy_P, smatrix& Piyy_M,
             smatrix& P_P, smatrix& P_M,
             double gamma, double eta,
             double zeta, double tau_nu){

    int N = rho_P.N;

    smatrix flux_Mass(N);
    smatrix flux_Momx(N);
    smatrix flux_Momy(N);
    smatrix flux_Pixx_vy(N);
    smatrix flux_Pixy_vy(N);
    smatrix flux_Piyx_vy(N);
    smatrix flux_Piyy_vy(N);

    double momy_av;
    double Piyy_av;
    double Pixy_av;
    double Pixx_vy_av;
    double Pixy_vy_av;
    double Piyx_vy_av;
    double Piyy_vy_av;
    double P_av;

    smatrix cs_P = getSpeedOfSound(rho_P, gamma);
    smatrix C_P = local_propagation_speed(rho_P, vx_P, vy_P, eta, zeta, tau_nu, cs_P);

    smatrix cs_M = getSpeedOfSound(rho_M, gamma);
    smatrix C_M = local_propagation_speed(rho_M, vx_M, vy_M, eta, zeta, tau_nu, cs_M);

    double C;

    double B = eta / tau_nu;
    double A = zeta/ tau_nu;


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {

            momy_av = 0.5 * (rho_P.get(i,j) * vy_P.get(i,j) + rho_M.get(i,j) * vy_M.get(i,j));
            Piyy_av = 0.5 * (Piyy_P.get(i,j) + Piyy_M.get(i,j));
            Pixy_av = 0.5 * (Pixy_P.get(i,j) + Pixy_M.get(i,j));
            Pixx_vy_av = 0.5 * (Pixx_P.get(i,j) * vy_P.get(i,j) + Pixx_M.get(i,j) * vy_M.get(i,j));
            Pixy_vy_av = 0.5 * (Pixy_P.get(i,j) * vy_P.get(i,j) + Piyx_M.get(i,j) * vy_M.get(i,j));
            Piyx_vy_av = 0.5 * (Piyx_P.get(i,j) * vy_P.get(i,j) + Pixy_M.get(i,j) * vy_M.get(i,j));
            Piyy_vy_av = 0.5 * (Piyy_P.get(i,j) * vy_P.get(i,j) + Piyy_M.get(i,j) * vy_M.get(i,j));
            P_av = 0.5 * (P_P.get(i,j) + P_M.get(i,j));

            C = std::max(C_M.get(i,j), C_P.get(i,j));

            flux_Mass.set(momy_av - C * 0.5 * (rho_P.get(i,j) - rho_M.get(i,j)),i,j);
            flux_Momx.set(0.5 * (rho_P.get(i,j) * vx_P.get(i,j) * vy_P.get(i,j) + rho_M.get(i,j) * vx_M.get(i,j) * vy_M.get(i,j)) + Pixy_av / gamma - C * 0.5 * (rho_P.get(i,j) * vx_P.get(i,j) - rho_M.get(i,j) * vx_M.get(i,j)),i,j);
            flux_Momy.set(0.5 * (rho_P.get(i,j) * vy_P.get(i,j) * vy_P.get(i,j) + rho_M.get(i,j) * vy_M.get(i,j) * vy_M.get(i,j)) + ((P_av)  + (Piyy_av)) / gamma - C * 0.5 * (rho_P.get(i,j) * vy_P.get(i,j) - rho_M.get(i,j) * vy_M.get(i,j)),i,j);
            flux_Pixx_vy.set(Pixx_vy_av + (A - 2.0 / 3.0 * B) * (vy_P.get(i,j) + vy_M.get(i,j)) * 0.5 - C * 0.5 * (Pixx_P.get(i,j) - Pixx_M.get(i,j)) ,i,j);
            flux_Pixy_vy.set(Pixy_vy_av + B * (vx_P.get(i,j) + vx_M.get(i,j)) * 0.5 - C * 0.5 * (Pixy_P.get(i,j) - Pixy_M.get(i,j)),i,j);
            flux_Piyx_vy.set(Piyx_vy_av + B * (vx_P.get(i,j) + vx_M.get(i,j)) * 0.5 - C * 0.5 * (Piyx_P.get(i,j) - Piyx_M.get(i,j)),i,j);
            flux_Piyy_vy.set(Piyy_vy_av + B * (vy_P.get(i,j) + vy_M.get(i,j)) + (A - 2.0 / 3.0 * B) * (vy_P.get(i,j) + vy_M.get(i,j)) * 0.5 - C * 0.5 * (Piyy_P.get(i,j) - Piyy_M.get(i,j)),i,j);


        }
    }

    return make_tuple(flux_Mass, flux_Momx, flux_Momy, flux_Pixx_vy, flux_Pixy_vy, flux_Piyx_vy, flux_Piyy_vy);
        
}
// Apply fluxes to conserved variables

smatrix applyFluxes(smatrix& flux_H1_X, smatrix& flux_H2_X,
  smatrix& flux_H1_Y, smatrix& flux_H2_Y,
   double dx, double dy, smatrix& J){

    int N = flux_H1_X.N;


    smatrix C(N);
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C.set(C.get(i,j) - (flux_H1_X.get(i,j) - flux_H2_X.get(i,j)) / dx,i,j);
            C.set(C.get(i,j) - (flux_H1_Y.get(i,j) - flux_H2_Y.get(i,j)) / dy,i,j);
            C.set(C.get(i,j) + J.get(i,j), i,j);
        }
    }

    return C;
}

// Heun's method
state heuns (state& q, function<state(double,state)> f, double dt, double t) {
    
    int N = q.get(0).N;

    smatrix k1(N);
    smatrix k2(N);
    smatrix c1(N);
    smatrix c2(N);
    smatrix y(N);
    smatrix yprime(N);


    state qprime(N);
    state C1(N),C2(N);

    C1 = f(t,q);

    for (int n = 0; n < 7; ++n){
        c1 = C1.get(n);
        y = q.get(n);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                k1.set(dt * c1.get(i,j),i,j);
                yprime.set(y.get(i,j) + k1.get(i,j),i,j);
            }
        }
        qprime.set(n , yprime);
    }
    
    C2 = f(t+dt,qprime);

    for (int n = 0; n < 7; ++n){
        c1 = C1.get(n);
        c2 = C2.get(n);
        y = q.get(n);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                k2.set(dt * (c2.get(i,j)),i,j);
                yprime.set(y.get(i,j) + 0.5*(k1.get(i,j) + k2.get(i,j)),i,j);
            }
        }
        qprime.set(n , yprime);
    }
    

    return qprime;
}
