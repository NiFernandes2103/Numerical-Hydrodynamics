// file: matrix.h

#ifndef _matrix_h
#define _matrix_h

# include <tuple>
#pragma once
using namespace std;


struct smatrix{

    int N;
    double* m;

    smatrix() : N(0), m(nullptr) {}

    smatrix(int n) : N(n), m(new double[n * n]{0}) {}

    double get(int i, int j) {
        return m[i + j*N];
    }

    void set(double var, int i, int j) {
        m[i + j*N] = var; 
    }

    ~smatrix() {
        delete[] m;
    }



};


#endif