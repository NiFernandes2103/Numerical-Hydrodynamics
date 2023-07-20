#include "KTmethods2d.h"
#include "matrix.h"
#include "fileFunc.h"
#include <iostream>
#include <map>
#include <fstream>
using namespace std;


void parameters_csv(double t, double tEnd, double tOut,int N, double boxsize, double a, double b, double gamma, double zeta, double eta, double tau_nu, double theta, string filename ) {
    
    // file pointer
    fstream fout;
  
    // opens an existing csv file or creates a new file.
    fout.open(filename, ios::out | ios::trunc);


    fout << "t" << ",";
    fout << "tEnd" << ",";
    fout << "tOut" << ",";
    fout << "N" << ",";
    fout << "boxsize" << ",";
    fout << "a" << ",";
    fout << "b" << ",";
    fout << "gamma" << ",";
    fout << "zeta" << ",";
    fout << "eta" << ",";
    fout << "tau_nu" << ",";
    fout << "theta" << endl;
    
    fout << t << ",";
    fout << tEnd << ",";
    fout << tOut << ",";
    fout << N << ",";
    fout << boxsize << ",";
    fout << a << ",";
    fout << b << ",";
    fout << gamma << ",";
    fout << zeta << ",";
    fout << eta << ",";
    fout << tau_nu << ",";
    fout << theta << endl;

    fout.close();
    
}


void create(state initial,string filename) {
    
    // file pointer
    fstream fout;
  
    // opens an existing csv file or creates a new file.
    fout.open(filename, ios::out | ios::trunc);


    smatrix rhoIC  = initial.get(0);        
    smatrix MomxIC = initial.get(1);
    smatrix MomyIC = initial.get(2);
    smatrix PixxIC = initial.get(3);
    smatrix PixyIC = initial.get(4);
    smatrix PiyxIC = initial.get(5);
    smatrix PiyyIC = initial.get(6);

    int A = rhoIC.N;


    // Read the input
    for (int i = 0; i < A; i++) {
        for (int j = 0; j < A; j++) {
            fout << rhoIC.get(i,j)  <<',';
            fout << MomxIC.get(i,j) <<',';
            fout << MomyIC.get(i,j) <<',';
            fout << PixxIC.get(i,j) <<',';
            fout << PixyIC.get(i,j) <<',';
            fout << PiyxIC.get(i,j) <<',';
            fout << PiyyIC.get(i,j) <<',';
        }
    }

    fout.close();
    
}


void write(map<double, state> solution, string filename)
{
    // file pointer
    fstream fout;

  
    // opens an existing csv file or creates a new file.
    fout.open(filename, ios::out | ios::trunc);


    // Iterate over the list using the iterator
    // Read the input
    for (map<double,state>::iterator it = solution.begin(); it != solution.end(); ++it) {
        smatrix rho  = (*it).second.get(0);        
        smatrix Momx = (*it).second.get(1);
        smatrix Momy = (*it).second.get(2);
        smatrix Pixx = (*it).second.get(3);
        smatrix Pixy = (*it).second.get(4);
        smatrix Piyx = (*it).second.get(5);
        smatrix Piyy = (*it).second.get(6);

        int N = rho.N;

        // write the csv file
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fout << rho.get(i,j)  <<',';
                fout << Momx.get(i,j) <<',';
                fout << Momy.get(i,j) <<',';
                fout << Pixx.get(i,j) <<',';
                fout << Pixy.get(i,j) <<',';
                fout << Piyx.get(i,j) <<',';
                fout << Piyy.get(i,j) <<',';
            }
        }
        
        fout << "\n";

        
    }

    fout.close();

}

void write_each(map<double, state> solution, string filename, int n)
{
    // file pointer
    fstream fout;

  
    // opens an existing csv file or creates a new file.
    fout.open(filename, ios::out | ios::trunc);


    // Iterate over the list using the iterator
    // Read the input
    for (map<double,state>::iterator it = solution.begin(); it != solution.end(); ++it) {
        smatrix v  = (*it).second.get(n);        

        int N = v.N;

        // write the csv file
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fout << v.get(i,j)  <<',';
            }
        }
        
        fout << "\n";
        
    }

    fout.close();

   
    
}