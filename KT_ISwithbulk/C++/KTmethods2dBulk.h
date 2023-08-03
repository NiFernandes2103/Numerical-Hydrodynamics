// file: KTmethods2d.h

#ifndef _KTmethods2dBulk_h
#define _KTmethods2dBulk_h

#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <functional>
#pragma once

class stateb {       
  private:             
    std::vector<std::vector<double>> rho;        
    std::vector<std::vector<double>> Momx;
    std::vector<std::vector<double>> Momy;
    std::vector<std::vector<double>> Pi;


public:
    stateb() : rho(), Momx(), Momy(), Pi() {}

    stateb(std::vector<std::vector<double>> v1,
     std::vector<std::vector<double>> v2,
     std::vector<std::vector<double>> v3,
     std::vector<std::vector<double>> v4) : rho(v1), Momx(v2), Momy(v3), Pi(v4) {}

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
            return Pi;
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
            Pi = v;
        }
    }

};


double max_value(const std::vector<std::vector<double>>& value);

std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>> getConserved(std::vector<std::vector<double>>& rho,
 std::vector<std::vector<double>>& vx, std::vector<std::vector<double>>& vy, double vol);

std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>,std::vector<std::vector<double>>> getPrimitive(std::vector<std::vector<double>> &Mass,
 std::vector<std::vector<double>>& Momx, std::vector<std::vector<double>>& Momy, double gamma, double vol);

std::vector<std::vector<double>> getSpeedOfSound(std::vector<std::vector<double>>& rho, double gamma);

double minmod(double x, double y, double z);

std::vector<std::vector<double>> getGradient (std::vector<std::vector<double>>& f, double dx, int axis, double theta);

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> extrapolateInSpaceToFace (std::vector<std::vector<double>>& q, std::vector<std::vector<double>>& q_dx, double dx, int axis);

std::vector<std::vector<double>> local_propagation_speed (std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& vx, std::vector<std::vector<double>>& vy, double zeta,  double tau_nu, std::vector<std::vector<double>>& cs);

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, 
std::vector<std::vector<double>>, std::vector<std::vector<double>>> getXFlux(std::vector<std::vector<double>>& rho_P,
    std::vector<std::vector<double>>& rho_M, std::vector<std::vector<double>>& vx_P, std::vector<std::vector<double>>& vx_M, 
    std::vector<std::vector<double>>& vy_P, std::vector<std::vector<double>>& vy_M, std::vector<std::vector<double>>& Pi_P, 
   std::vector<std::vector<double>>& Pi_M, std::vector<std::vector<double>>& P_P, std::vector<std::vector<double>>& P_M, double gamma, 
  double zeta, double tau_nu);

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>,
 std::vector<std::vector<double>>, std::vector<std::vector<double>>> getYFlux(std::vector<std::vector<double>>& rho_P, std::vector<std::vector<double>>& rho_M,
             std::vector<std::vector<double>>& vx_P, std::vector<std::vector<double>>& vx_M,
             std::vector<std::vector<double>>& vy_P, std::vector<std::vector<double>>& vy_M,
             std::vector<std::vector<double>>& Pi_P, std::vector<std::vector<double>>& Pi_M,
             std::vector<std::vector<double>>& P_P, std::vector<std::vector<double>>& P_M,
             double gamma, double zeta, double tau_nu);

// Apply fluxes to conserved variables
std::vector<std::vector<double>> applyFluxes(std::vector<std::vector<double>>& flux_H1_X, std::vector<std::vector<double>>& flux_H2_X,
  std::vector<std::vector<double>>& flux_H1_Y, std::vector<std::vector<double>>& flux_H2_Y,
   double dx, double dy, std::vector<std::vector<double>>& J);

// Heun's method
stateb heuns (stateb& q, std::function<stateb(double,stateb)> f, double dt, double t);

stateb rK4 (stateb& q, std::function<stateb(double,stateb)> f, double dt, double t);

# endif