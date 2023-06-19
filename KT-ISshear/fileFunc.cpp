#include "KTmethods2d.cpp"
#include <iostream>
#include <fstream>

void create(Solution &solution)
{
    // file pointer
    fstream fout;
  
    // opens an existing csv file or creates a new file.
    fout.open("NonRelatISwithShear.csv", ios::out | ios::app);

    std::list<State>::iterator it;

    // Iterate over the list using the iterator
    for (it = solution.begin(); it != solution.end(); ++it) {
        vector<vector<double>> rho  = *it.get(1);        
        vector<vector<double>> Momx = *it.get(2);
        vector<vector<double>> Momy = *it.get(3);
        vector<vector<double>> Pixx = *it.get(4);
        vector<vector<double>> Pixy = *it.get(5);
        vector<vector<double>> Piyx = *it.get(6);
        vector<vector<double>> Piyy = *it.get(7);

        int N = rho.size();

        // Read the input
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
        
    }

    
    
}