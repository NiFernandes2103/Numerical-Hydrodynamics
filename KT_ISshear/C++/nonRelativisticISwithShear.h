// Purpose: Header file for nonRelativisticISwithShear.cpp

#pragma once
#include <cmath>
#include <vector>
#include <tuple>
#include <list>
#include <iostream>
#include "KTmethods2d.h"


state KTschemeNonRelativisticIS(double t,  state& IC, double dx, double dy, int N, double gamma, double zeta, double tau_nu, double eta, double theta);

std::list<state> integrator(state (*scheme)(double, state&, double, double, int, double, double, double, double, double), std::tuple<double,double> time, std::list<state> Q, double dtmax,  std::tuple<double, double, int, double, double, double, double, double> args, std::string method = "Heuns");