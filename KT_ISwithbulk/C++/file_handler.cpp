#include "KTmethods2dBulk.h"
#include "file_handler.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

void parameters_csv(double t, double tEnd, double tOut, int N, double boxsize, double a, double b, double gamma, double zeta, double tau_nu, double theta, std::string filename ) {
    
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
    fout << tau_nu << ",";
    fout << theta << std::endl;

    fout.close();
    
}


void create(stateb initial,std::string filename) {
    
    // file pointer
    std::fstream fout;
  
    // opens an existing csv file or creates a new file.
    fout.open(filename, std::ios::out | std::ios::trunc);


    std::vector<std::vector<double>> rhoIC  = initial.get(0);        
    std::vector<std::vector<double>> MomxIC = initial.get(1);
    std::vector<std::vector<double>> MomyIC = initial.get(2);
    std::vector<std::vector<double>> PiIC = initial.get(3);

    int A = rhoIC.size();

    // Read the input
    for (int i = 0; i < A; i++) {
        for (int j = 0; j < A; j++) {
            fout << rhoIC[i][j]  <<',';
            fout << MomxIC[i][j] <<',';
            fout << MomyIC[i][j] <<',';
            fout << PiIC[i][j] <<',';

        }
    }

    fout.close();
    
}


void write(std::list<stateb> solution, std::string filename)
{
    // file pointer
    std::fstream fout;

  
    // opens an existing csv file or creates a new file.
    fout.open(filename, std::ios::out | std::ios::trunc);


    // Iterate over the list using the iterator
    // Read the input
    for (std::list<stateb>::iterator it = solution.begin(); it != solution.end(); ++it) {
        std::vector<std::vector<double>> rho  = (*it).get(0);        
        std::vector<std::vector<double>> Momx = (*it).get(1);
        std::vector<std::vector<double>> Momy = (*it).get(2);
        std::vector<std::vector<double>> Pi   = (*it).get(3);

        int N = rho.size();

        // write the csv file
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fout << rho[i][j]  <<',';
                fout << Momx[i][j] <<',';
                fout << Momy[i][j] <<',';
                fout << Pi[i][j] <<',';
  
            }
        }
        
        fout << "\n";

        
    }

    fout.close();

}

void write_each(std::list<stateb> solution, std::string filename, int n)
{
    // file pointer
    std::fstream fout;

  
    // opens an existing csv file or creates a new file.
    fout.open(filename, std::ios::out | std::ios::trunc);


    // Iterate over the list using the iterator
    // Read the input
    for (std::list<stateb>::iterator it = solution.begin(); it != solution.end(); ++it) {
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