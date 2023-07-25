#include "KTmethods2d.h"
#include "fileFunc.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

void parameters_csv(double t, double tEnd, double tOut, int N, double boxsize, double a, double b, double gamma, double zeta, double eta, double tau_nu, double theta, std::string filename ) {
    
    // file pointer
    std::fstream fout;
  
    // opens an existing csv file or creates a new file.
    fout.open(filename, std::ios::out | std::ios::trunc);


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
    fout << "theta" << std::endl;
    
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
    fout << theta << std::endl;

    fout.close();
    
}


void create(state initial,std::string filename) {
    
    // file pointer
    std::fstream fout;
  
    // opens an existing csv file or creates a new file.
    fout.open(filename, std::ios::out | std::ios::trunc);


    std::vector<std::vector<double>> rhoIC  = initial.get(0);        
    std::vector<std::vector<double>> MomxIC = initial.get(1);
    std::vector<std::vector<double>> MomyIC = initial.get(2);
    std::vector<std::vector<double>> PixxIC = initial.get(3);
    std::vector<std::vector<double>> PixyIC = initial.get(4);
    std::vector<std::vector<double>> PiyxIC = initial.get(5);
    std::vector<std::vector<double>> PiyyIC = initial.get(6);

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


void write(std::list<state> solution, std::string filename)
{
    // file pointer
    std::fstream fout;

  
    // opens an existing csv file or creates a new file.
    fout.open(filename, std::ios::out | std::ios::trunc);


    // Iterate over the list using the iterator
    // Read the input
    for (std::list<state>::iterator it = solution.begin(); it != solution.end(); ++it) {
        std::vector<std::vector<double>> rho  = (*it).get(0);        
        std::vector<std::vector<double>> Momx = (*it).get(1);
        std::vector<std::vector<double>> Momy = (*it).get(2);
        std::vector<std::vector<double>> Pixx = (*it).get(3);
        std::vector<std::vector<double>> Pixy = (*it).get(4);
        std::vector<std::vector<double>> Piyx = (*it).get(5);
        std::vector<std::vector<double>> Piyy = (*it).get(6);

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

void write_each(std::list<state> solution, std::string filename, int n)
{
    // file pointer
    std::fstream fout;

  
    // opens an existing csv file or creates a new file.
    fout.open(filename, std::ios::out | std::ios::trunc);


    // Iterate over the list using the iterator
    // Read the input
    for (std::list<state>::iterator it = solution.begin(); it != solution.end(); ++it) {
        std::vector<std::vector<double>> v  = (*it).get(n);        

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