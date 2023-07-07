// Purpose: Header file for fileFunc.cpp

#ifndef _fileFunc_h
#define _fileFunc_h

#include "KTmethods2d.h"
#include <iostream>
#include <map>
#include <fstream>
using namespace std;

void parameters_csv(double t, double tEnd, double tOut,int N, double boxsize, double gamma, double zeta, double eta, double tau_nu, double theta, string filename);

void create(state initial, string filename);

void write(map<double, state> solution, string filename);

void write_each(map<double, state> solution, string filename, int n);

#endif