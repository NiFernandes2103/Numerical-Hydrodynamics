// file: matrix.h

#ifndef _matrix_h
#define _matrix_h

# include <tuple>
#pragma once
using namespace std;


struct smatrix{

    int N;
    int M;
    double* arr = new double[N*M];

    smatrix() : N(0),M(0),arr({0}) {}

    smatrix(int n , int m) : N(n),M(m),arr({0}) {}

    smatrix(int n , int m, double* arr) : N(n),M(m),arr(arr) {}

    tuple<int, int> shape() {
        return make_tuple(N,M);
    }

    double get(int i, int j) {
        return arr[i + N*j];
    }

    void set(double var, int i, int j) {
        arr[i + N*j] = var; 
    }



};


#endif