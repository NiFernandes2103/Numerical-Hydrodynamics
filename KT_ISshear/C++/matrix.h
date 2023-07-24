// file: matrix.h

#ifndef _matrix_h
#define _matrix_h

# include <tuple>
#pragma once
using namespace std;


struct smatrix{

    double* m;
    int N;

    smatrix(smatrix &old) : m(old.m), N(old.N) {}

    smatrix(int n) : m(new double[n*n]{0}), N(n) {}

    double get(int i, int j) {
        return m[i*N + j];
    }

    void set(double var, int i, int j) {
        m[i*N + j] = var; 
    }

    ~smatrix() {
        delete[] m;
    }



};


#endif