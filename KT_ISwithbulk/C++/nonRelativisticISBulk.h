// Purpose: Header file for nonRelativisticISwithShear.cpp

#pragma once
#include <cmath>
#include <vector>
#include <tuple>
#include <list>
#include <iostream>
#include "KTmethods2dBulk.h"


stateb KTschemeNonRelativisticIS(double t,  stateb& IC, double dx, double dy, int N, double gamma, double zeta, double tau_nu, double theta);

std::list<stateb> integrator(stateb (*scheme)(double, stateb&, double, double, int, double, double, double, double), std::tuple<double,double> time, std::list<stateb> Q, double dtmax,  std::tuple<double, double, int, double, double, double, double> args, std::string method = "Heuns");