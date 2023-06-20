#include "KTmethods2d.cpp"
#include <iostream>
#include <map>
#include <fstream>

void create(map<double, State> &solution)
{
    // file pointer
    fstream fout;
  
    // opens an existing csv file or creates a new file.
    fout.open("NonRelatISwithShear.csv", ios::out | ios::app);


    // Iterate over the list using the iterator
    for (map<double,State>::iterator it = solution.begin(); it != solution.end(); ++it) {
        vector<vector<double>> rho  = (*it).second.get(1);        
        vector<vector<double>> Momx = (*it).second.get(2);
        vector<vector<double>> Momy = (*it).second.get(3);
        vector<vector<double>> Pixx = (*it).second.get(4);
        vector<vector<double>> Pixy = (*it).second.get(5);
        vector<vector<double>> Piyx = (*it).second.get(6);
        vector<vector<double>> Piyy = (*it).second.get(7);

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