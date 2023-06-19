
#include KTmethods.h

void create(Solution &solution)
{
    // file pointer
    fstream fout;
  
    // opens an existing csv file or creates a new file.
    fout.open("NonRelatISwithShear.csv", ios::out | ios::app);
  
    
    vector<vector<double>> rho;        
    vector<vector<double>> Momx;
    vector<vector<double>> Momy;
    vector<vector<double>> Pixx;
    vector<vector<double>> Pixy;
    vector<vector<double>> Piyx;
    vector<vector<double>> Piyy;
  
    // Read the input
    for (i = 0; i < 5; i++) {
  
        
  
        // Insert the data to file
        fout << roll << ", "
             << name << ", "
             << math << ", "
             << phy << ", "
             << chem << ", "
             << bio
             << "\n";
    }
}