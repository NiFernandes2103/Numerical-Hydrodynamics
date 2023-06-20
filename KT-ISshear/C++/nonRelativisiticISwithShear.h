
#include <cmath>
#include <vector>
#include <tuple>
#include <list>
#include <iostream>
#include <map> 
#include "KTmethods2d.h"
using namespace std;

state KTschemeNonRelativisticIS(double t,  state& IC, double dx, double dy, int N, double gamma, double zeta, double tau_nu, double eta, double theta);

map<double,state> integrator(state (*scheme)(double, state&, double, double, int, double, double, double, double, double), tuple<double,double> time, map<double, state> Q, double dtmax,  tuple<double, double, int, double, double, double, double, double> args, string method = "Heuns");