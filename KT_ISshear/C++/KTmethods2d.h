// file: KTmethods2d.h

#ifndef _KTmethods2d_h
#define _KTmethods2d_h

#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <functional>
#pragma once
using namespace std;

class state {       
  private:             
    double** rho;        
    double** Momx;
    double** Momy;
    double** Pixx;
    double** Pixy;
    double** Piyx;
    double** Piyy;

public:
    state() : rho(), Momx(), Momy(), Pixx(), Pixy(), Piyx(), Piyy() {}

    state(double** v1,
     double** v2,
     double** v3,
     double** v4,
     double** v5,
     double** v6,
     double** v7) : rho(v1), Momx(v2), Momy(v3), Pixx(v4), Pixy(v5), Piyx(v6), Piyy(v7) {}

    double** get(int n) {
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
            return nullptr;
        }
        
    }

    void set(int n , double** v) {

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


double max_value(double** value, int rows, int cols);

double** sign(double** value,  int rows, int cols);

tuple<double**,double**,double**> getConserved(double**& rho,
 double**& vx, double**& vy, double vol, int rows, int cols);

tuple<double**,double**,double**,double**> getPrimitive(double** &Mass,
 double**& Momx, double**& Momy, double gamma, double vol, int rows, int cols);

double** getSpeedOfSound(double**& rho, double gamma, int rows, int cols);

double minmod(double x, double y);

double minmod3(double x, double y, double z);

double** getGradient (double**& f, double dx, int axis, int rows, int cols, double theta);

tuple<double**, double**, double**, double**> extrapolateInSpaceToFace (double**& q, double**& q_dx, double dx, int axis, int rows, int cols);

double** local_propagation_speed (double**& rho, double**& vx, double**& vy, double eta, double zeta, double tau_nu, double**& cs, int rows, int cols);


tuple<double**, double**, double**, double**,
 double**, double**, double**> getXFlux(double**& rho_P,
    double**& rho_M, double**& vx_P, double**& vx_M, 
    double**& vy_P, double**& vy_M, double**& Pixx_P, 
   double**& Pixx_M, double**& Pixy_P, double**& Pixy_M, 
   double**& Piyx_P, double**& Piyx_M, double**& Piyy_P, 
   double**& Piyy_M, double**& P_P, double**& P_M, double gamma, 
  double eta, double zeta, double tau_nu, int rows, int cols);

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
             double zeta, double tau_nu, int rows, int cols);

// Apply fluxes to conserved variables
double** applyFluxes(double**& flux_H1_X, double**& flux_H2_X,
  double**& flux_H1_Y, double**& flux_H2_Y,
   double dx, double dy, double**& J, int rows, int cols);

// Heun's method
state heuns (state& q, function<state(double,state)> f, double dt, double t, int rows, int cols);

state rK4 (state& q, function<state(double,state)> f, double dt, double t, int rows, int cols);

# endif