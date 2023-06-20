// Purpose: Header file for fileFunc.cpp

#ifndef _fileFunc_h
#define _fileFunc_h

#include "KTmethods2d.h"
#include <iostream>
#include <map>
#include <fstream>
using namespace std;

void create(state initial, string filename);

void write(map<double, state> solution, string filename);

#endif