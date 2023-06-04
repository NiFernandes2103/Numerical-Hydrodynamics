
#ifndef KTmethods2d.h  
#define KTmethods2d.h
#include <vector>



int sign(vector<vector<double>> value);

tuple<vector<vector<double>>,vector<vector<double>>,vector<vector<double>>> getConserved(vector<vector<double>>& rho,
 vector<vector<double>>& vx, vector<vector<double>>& vy, double vol);

vector<vector<double>> getSpeedOfSound(vector<vector<double>>& rho, double& gamma); 

vector<vector<double>> minmod2(vector<vector<double>>& x, vector<vector<double>>& y);

vector<vector<double>> minmod3(vector<vector<double>>& x, vector<vector<double>>& y, vector<vector<double>>& z);

vector<vector<double>> getGradient (vector<vector<double>>& f, double dx, int axis, double theta = 1);

tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> extrapolateInSpaceToFace (vector<vector<double>>& q, vector<vector<double>>& q_dx, double dx, int axis);

vector<vector<double>> local_propagation_speed (vector<vector<double>>& rho, vector<vector<double>>& vx, vector<vector<double>>& vy, double eta, double zeta, double tau_nu, vector<vector<double>>& cs);

tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>,
 vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> getXFlux(vector<vector<double>>& rho_P,
    vector<vector<double>>& rho_M, vector<vector<double>>& vx_P, vector<vector<double>>& vx_M, 
    vector<vector<double>>& vy_P, vector<vector<double>>& vy_M, vector<vector<double>>& Pixx_P, 
   vector<vector<double>>& Pixx_M, vector<vector<double>>& Pixy_P, vector<vector<double>>& Pixy_M, 
   vector<vector<double>>& Piyx_P, vector<vector<double>>& Piyx_M, vector<vector<double>>& Piyy_P, 
   vector<vector<double>>& Piyy_M, vector<vector<double>>& P_P, vector<vector<double>>& P_M, double gamma, 
  double eta, double zeta, double tau_nu);

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
             double zeta, double tau_nu); 
#endif