// file: KTmethods2d.h

#ifndef _KTmethods2d_h
#define _KTmethods2d_h

#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <functional>
#pragma once

class state {       
  private:             
    std::vector<std::vector<double>> rho;        
    std::vector<std::vector<double>> Momx;
    std::vector<std::vector<double>> Momy;
    std::vector<std::vector<double>> Pixx;
    std::vector<std::vector<double>> Pixy;
    std::vector<std::vector<double>> Piyx;
    std::vector<std::vector<double>> Piyy;

public:
    state() : rho(), Momx(), Momy(), Pixx(), Pixy(), Piyx(), Piyy() {}

    state(std::vector<std::vector<double>> v1,
     std::vector<std::vector<double>> v2,
     std::vector<std::vector<double>> v3,
     std::vector<std::vector<double>> v4,
     std::vector<std::vector<double>> v5,
     std::vector<std::vector<double>> v6,
     std::vector<std::vector<double>> v7) : rho(v1), Momx(v2), Momy(v3), Pixx(v4), Pixy(v5), Piyx(v6), Piyy(v7) {}

    std::vector<std::vector<double>> get(int n) {
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
            return std::vector<std::vector<double>>();
        } 
        
    }

    void set(int n , std::vector<std::vector<double>> v) {

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


double max_value(const std::vector<std::vector<double>>& value);

std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>> getConserved(std::vector<std::vector<double>>& rho,
 std::vector<std::vector<double>>& vx, std::vector<std::vector<double>>& vy, double vol);

std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>> getPrimitive(std::vector<std::vector<double>> &Mass,
 std::vector<std::vector<double>>& Momx, std::vector<std::vector<double>>& Momy, double gamma, double vol);

std::vector<std::vector<double>> getSpeedOfSound(std::vector<std::vector<double>>& rho, double gamma);

double minmod2(double x, double y);

double minmod3(double x, double y, double z);

double minmod(double x, double y, double z);

std::vector<std::vector<double>> getGradient (std::vector<std::vector<double>>& f, double dx, int axis, double theta);

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> extrapolateInSpaceToFace (std::vector<std::vector<double>>& q, std::vector<std::vector<double>>& q_dx, double dx, int axis);

std::vector<std::vector<double>> local_propagation_speed (std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& vx, std::vector<std::vector<double>>& vy, double eta, double zeta, double tau_nu, std::vector<std::vector<double>>& cs);

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>,
 std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> getXFlux(std::vector<std::vector<double>>& rho_P,
    std::vector<std::vector<double>>& rho_M, std::vector<std::vector<double>>& vx_P, std::vector<std::vector<double>>& vx_M, 
    std::vector<std::vector<double>>& vy_P, std::vector<std::vector<double>>& vy_M, std::vector<std::vector<double>>& Pixx_P, 
   std::vector<std::vector<double>>& Pixx_M, std::vector<std::vector<double>>& Pixy_P, std::vector<std::vector<double>>& Pixy_M, 
   std::vector<std::vector<double>>& Piyx_P, std::vector<std::vector<double>>& Piyx_M, std::vector<std::vector<double>>& Piyy_P, 
   std::vector<std::vector<double>>& Piyy_M, std::vector<std::vector<double>>& P_P, std::vector<std::vector<double>>& P_M, double gamma, 
  double eta, double zeta, double tau_nu);

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
             double zeta, double tau_nu);

// Apply fluxes to conserved variables
std::vector<std::vector<double>> applyFluxes(std::vector<std::vector<double>>& flux_H1_X, std::vector<std::vector<double>>& flux_H2_X,
  std::vector<std::vector<double>>& flux_H1_Y, std::vector<std::vector<double>>& flux_H2_Y,
   double dx, double dy, std::vector<std::vector<double>>& J);

// Heun's method
state heuns (state& q, std::function<state(double,state)> f, double dt, double t);

state rK4 (state& q, std::function<state(double,state)> f, double dt, double t);

# endif