// file: matrix.h

#ifndef _matrix_h
#define _matrix_h

# include <algorithm>
#pragma once


struct smatrix{

    double* m;
    int N;

    smatrix(smatrix &old) : N(old.N) {
        m = new double[N * N];
        std::copy(old.m, old.m + N * N, m);   
     }

    // Copy Constructor
    smatrix(const smatrix& old) : N(old.N) {
        m = new double[N * N];
        std::copy(old.m, old.m + N * N, m);
    }

    // Copy Assignment Operator
    smatrix& operator=(const smatrix& old) {
        if (this != &old) {
            delete[] m;
            N = old.N;
            m = new double[N * N];
            std::copy(old.m, old.m + N * N, m);
        }
        return *this;
    }

    smatrix(int n) : m(new double[n*n]{0}), N(n) {}

    double get(int i, int j) {
        return m[i*N + j];
    }

    void set(double var, int i, int j) {
        m[i*N + j] = var; 
    }

    double max() {
        return *std::max_element(m, m + N*N);
    }

    

    ~smatrix() {
        delete[] m;
    }



};


#endif