// Purpose: Header file for fileFunc.cpp

#ifndef _file_handler_h
#define _file_handler_h

#include "KTmethods2dBulk.h"
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <fstream>


void parameters_csv(double t, double tEnd, double tOut,int N, double boxsize, double a, double b, double gamma, double zeta, double tau_nu, double theta, std::string filename);

void create(state initial, std::string filename);

void write(std::list<state> solution, std::string filename);

void write_each(std::list<state> solution, std::string filename, int n);

#endif