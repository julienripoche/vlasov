#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"

using namespace std;

int main()
{
    //Initialize positions and momenta values
    vector<vector<double> > r;
    vector<vector<double> > p;
    read(r, "coords.gnu");
    read(p, "momenta.gnu");

    //Useful variables
    int N = 5;
    int n_box = 30;
    double L = 20.;
    double l0 = L/n_box;
    double sigma = 1.;

    //Initialize rho map
    vector<double> rho(n_box*n_box*n_box,0);
    rho_map(rho, r, sigma, l0, n_box, N);

    //Write density profile in "density.gnu" 
    ofstream densityFile("density.gnu");
    for(int i=0 ; i<n_box ; i++)
    {
        for(int j=0 ; j<n_box ; j++)
        {
            densityFile << i << " " << j << " " << rho[key(i,j,n_box/2,n_box)] << endl;
        }
    }

    return 0;
}
