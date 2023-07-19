// file: KTmethods2d.h

#ifndef _KTmethods2d_h
#define _KTmethods2d_h

#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <map>
#include <functional>
#pragma once
using namespace std;

class state {       
  private:             
    vector<vector<double>> rho;        
    vector<vector<double>> Momx;
    vector<vector<double>> Momy;
    vector<vector<double>> Pixx;
    vector<vector<double>> Pixy;
    vector<vector<double>> Piyx;
    vector<vector<double>> Piyy;

public:
    state() : rho(), Momx(), Momy(), Pixx(), Pixy(), Piyx(), Piyy() {}

    state(vector<vector<double>> v1,
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




double max_value(vector<vector<double>> value);

double sign(double value);

tuple<vector<vector<double>>,vector<vector<double>>,vector<vector<double>>> getConserved(vector<vector<double>>& rho,
 vector<vector<double>>& vx, vector<vector<double>>& vy, double vol);

tuple<vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,vector<vector<double>>> getPrimitive(vector<vector<double>> &Mass,
 vector<vector<double>>& Momx, vector<vector<double>>& Momy, double gamma, double vol);

vector<vector<double>> getSpeedOfSound(vector<vector<double>>& rho, double gamma);

double minmod(double a, double b);

double minmod3(double a, double b, double c);

vector<vector<double>> getGradient (vector<vector<double>>& f, double dx, int axis, double theta);

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

// Apply fluxes to conserved variables
vector<vector<double>> applyFluxes(vector<vector<double>>& flux_H1_X, vector<vector<double>>& flux_H2_X,
  vector<vector<double>>& flux_H1_Y, vector<vector<double>>& flux_H2_Y,
   double dx, double dy, vector<vector<double>>& J);

// Heun's method
state heuns (state& q, function<state(double,state)> f, double dt, double t);



# endif