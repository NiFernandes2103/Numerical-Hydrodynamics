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
#include "matrix.h"
#pragma once
using namespace std;

class state {       
  private:             
    smatrix rho;        
    smatrix Momx;
    smatrix Momy;
    smatrix Pixx;
    smatrix Pixy;
    smatrix Piyx;
    smatrix Piyy;

public:
    state() : rho(), Momx(), Momy(), Pixx(), Pixy(), Piyx(), Piyy() {}

    state(smatrix v1,
     smatrix v2,
     smatrix v3,
     smatrix v4,
     smatrix v5,
     smatrix v6,
     smatrix v7) : rho(v1), Momx(v2), Momy(v3), Pixx(v4), Pixy(v5), Piyx(v6), Piyy(v7) {}

    smatrix get(int n) {
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
            return smatrix();
        } 
        
    }

    void set(int n , smatrix v) {

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




double max_value(smatrix value);

double sign(double value);

tuple<smatrix,smatrix,smatrix> getConserved(smatrix& rho,
 smatrix& vx, smatrix& vy, double vol);

tuple<smatrix,smatrix,smatrix,smatrix> getPrimitive(smatrix &Mass,
 smatrix& Momx, smatrix& Momy, double gamma, double vol);

smatrix getSpeedOfSound(smatrix& rho, double gamma);

double minmod(double a, double b);

double minmod3(double a, double b, double c);

smatrix getGradient (smatrix& f, double dx, int axis, double theta);

tuple<smatrix, smatrix, smatrix, smatrix> extrapolateInSpaceToFace (smatrix& q, smatrix& q_dx, double dx, int axis);

smatrix local_propagation_speed (smatrix& rho, smatrix& vx, smatrix& vy, double eta, double zeta, double tau_nu, smatrix& cs);


tuple<smatrix, smatrix, smatrix, smatrix,
 smatrix, smatrix, smatrix> getXFlux(smatrix& rho_P,
    smatrix& rho_M, smatrix& vx_P, smatrix& vx_M, 
    smatrix& vy_P, smatrix& vy_M, smatrix& Pixx_P, 
   smatrix& Pixx_M, smatrix& Pixy_P, smatrix& Pixy_M, 
   smatrix& Piyx_P, smatrix& Piyx_M, smatrix& Piyy_P, 
   smatrix& Piyy_M, smatrix& P_P, smatrix& P_M, double gamma, 
  double eta, double zeta, double tau_nu);

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
             double zeta, double tau_nu);

// Apply fluxes to conserved variables
smatrix applyFluxes(smatrix& flux_H1_X, smatrix& flux_H2_X,
  smatrix& flux_H1_Y, smatrix& flux_H2_Y,
   double dx, double dy, smatrix& J);

// Heun's method
state heuns (state& q, function<state(double,state)> f, double dt, double t);



# endif