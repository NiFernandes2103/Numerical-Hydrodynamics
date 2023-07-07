#include "KTmethods2d.h"
#include "fileFunc.h"
#include <iostream>
#include <map>
#include <fstream>
using namespace std;


void parameters_csv(double t, double tEnd, double tOut,int N, double boxsize, double gamma, double zeta, double eta, double tau_nu, double theta, string filename ) {
    
    // file pointer
    fstream fout;
  
    // opens an existing csv file or creates a new file.
    fout.open(filename, ios::out | ios::trunc);


    fout << "t" << ",";
    fout << "tEnd" << ",";
    fout << "tOut" << ",";
    fout << "N" << ",";
    fout << "boxsize" << ",";
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


    vector<vector<double>> rhoIC  = initial.get(0);        
    vector<vector<double>> MomxIC = initial.get(1);
    vector<vector<double>> MomyIC = initial.get(2);
    vector<vector<double>> PixxIC = initial.get(3);
    vector<vector<double>> PixyIC = initial.get(4);
    vector<vector<double>> PiyxIC = initial.get(5);
    vector<vector<double>> PiyyIC = initial.get(6);

    int A = rhoIC.size();

    // Read the input
    for (int i = 0; i < A; i++) {
        for (int j = 0; j < A; j++) {
            fout << rhoIC[i][j]  <<',';
            fout << MomxIC[i][j] <<',';
            fout << MomyIC[i][j] <<',';
            fout << PixxIC[i][j] <<',';
            fout << PixyIC[i][j] <<',';
            fout << PiyxIC[i][j] <<',';
            fout << PiyyIC[i][j] <<',';
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
        vector<vector<double>> rho  = (*it).second.get(0);        
        vector<vector<double>> Momx = (*it).second.get(1);
        vector<vector<double>> Momy = (*it).second.get(2);
        vector<vector<double>> Pixx = (*it).second.get(3);
        vector<vector<double>> Pixy = (*it).second.get(4);
        vector<vector<double>> Piyx = (*it).second.get(5);
        vector<vector<double>> Piyy = (*it).second.get(6);

        int N = rho.size();

        // write the csv file
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fout << rho[i][j]  <<',';
                fout << Momx[i][j] <<',';
                fout << Momy[i][j] <<',';
                fout << Pixx[i][j] <<',';
                fout << Pixy[i][j] <<',';
                fout << Piyx[i][j] <<',';
                fout << Piyy[i][j] <<',';
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
        vector<vector<double>> v  = (*it).second.get(n);        

        int N = v.size();

        // write the csv file
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fout << v[i][j]  <<',';
            }
        }
        
        fout << "\n";
        
    }

    fout.close();



    
    
}